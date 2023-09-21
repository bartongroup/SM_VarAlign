import itertools
import logging
import os
from copy import deepcopy

import pandas as pd
import tqdm
import vcf  ## Requires pysam functionality
from vcf.utils import trim_common_suffix

log = logging.getLogger(__name__)


class Reader(vcf.Reader):
    """
    Extension of PyVCF vcf.Reader class that parses VCF data to a Pandas DataFrame.

    VCF compatibility requirements:
        - VEP annotated
        - VEP ALLELE_NUM annotation
    """

    def __init__(self, *args, **kwargs):
        super(Reader, self).__init__(*args, **kwargs)

        # Override default ascii encoding for better compatibility
        if 'encoding' not in kwargs.keys():
            self.encoding = 'utf_8'

        # Log some key metadata
        log.info('Loading %s...', self.filename)
        log.info('File format: {}'.format(self.metadata.get('fileformat')))
        log.info('Reference: {}'.format(self.metadata.get('reference')))

        # Parse VEP annotation format
        CSQ_Format = self.infos['CSQ'].desc.split(' Format: ')[1].split('|')
        log.info('CSQ Format: {}'.format(CSQ_Format))
        self.CSQ_Format = CSQ_Format

        # Parse INFO header
        info_header = pd.DataFrame([v for k, v in list(self.infos.items())], dtype='str')
        # PyVCF parses the VCF codes to int or None as follows:
        # field_counts = {
        #     '.': None,  # Unknown number of values
        #     'A': -1,  # Equal to the number of alternate alleles in a given record
        #     'G': -2,  # Equal to the number of genotypes in a given record
        #     'R': -3,  # Equal to the number of alleles including reference in a given record
        # }
        log.info('gnomAD header:\n%s', info_header.head().to_string())
        log.info('gnomAD header types:\n%s', info_header.dtypes.to_string())
        info_flag_num = '0.0'
        info_value_num = '1.0'
        info_allele_num = '-1.0'
        info_flag_fields = info_header.query('num == @info_flag_num').id.tolist()
        log.info('gnomAD info value header:\n%s', info_header.query('num == @info_value_num').to_string())
        info_value_fields = info_header.query('num == @info_value_num').id.tolist()
        info_allele_fields = info_header.query('num == @info_allele_num').id.tolist()
        info_other_fields = set(info_header[info_header['num'].isnull()].id)
        info_other_fields.discard('CSQ')
        log.info('INFO flags: {}'.format(str(info_flag_fields)))
        log.info('INFO value: {}'.format(str(info_value_fields)))
        log.info('INFO per allele: {}'.format(str(info_allele_fields)))
        log.info('Other INFO fields: {}'.format(str(info_other_fields)))
        self.info_flag_fields = info_flag_fields
        self.info_value_fields = info_value_fields
        self.info_allele_fields = info_allele_fields
        self.info_other_fields = info_other_fields
        self.info_header = info_header

        # Annotations that need special handling during variant allele expansion
        standard_num_values = [info_flag_num, info_value_num, info_allele_num]
        special_handling = {'INFO': list(info_header.query('num not in @standard_num_values').id)}
        msg = 'The following INFO fields cannot be resolved to a variant and will be stripped when melting variants: {}'
        log.info(msg.format(str(special_handling['INFO'])))
        self.special_handling = special_handling

        # Parse FORMAT header
        format_header = pd.DataFrame([v for k, v in list(self.formats.items())], dtype='str')
        try:
            log.info('FORMAT fields: {}'.format(str(list(format_header.id))))
        except AttributeError:
            log.info('No FORMAT header.')
        self.format_header = format_header

        # Parse FILTER header
        filter_header = pd.DataFrame([v for k, v in list(self.filters.items())], dtype='str')
        try:
            log.info('FILTER flags: {}'.format(str(list(filter_header.id))))
        except AttributeError:
            log.info('No FILTER header.')
        self.filter_header = filter_header

    def get_vep_annotation(self, record, fields=None):
        """
        Retrieve and format the VEP annotation from a VCF record.
        """
        if fields is None:
            fields = self.CSQ_Format
        parsed = pd.DataFrame([x.split('|') for x in record.INFO['CSQ']],
                              columns=self.CSQ_Format)
        if len(fields) < len(self.CSQ_Format):
            parsed = parsed[fields]

        return parsed

    @staticmethod
    def get_vep_raw(record):
        """
        Retrive VEP annotation from VCF record with minimal parsing.
        """
        return [x.split('|') for x in record.INFO['CSQ']]

    def tabulate_variant_effects(self, variants):
        """
        Efficiently parse the VEP CSQ field from a list of VCF records into a table.

        :param variants: List of VCF records.
        :return: DataFrame of VEP annotations indexed to variant list.
        """
        veps = [self.get_vep_raw(x) for x in variants]
        variant_indices = [i for i, x in enumerate(veps) for _ in x]  # index each vep entry to the variant record
        vep_table = pd.DataFrame(list(itertools.chain.from_iterable(veps)),
                                 columns=self.CSQ_Format, index=variant_indices)
        return vep_table

    def tabulate_variant_info(self, variants, split=False):
        """
        Efficiently parse INFO from a list of VCF records into a table.

        :param variants: List of VCF records.
        :param split: Whether to return seperate tables for site and allele annotations.
        :return: DataFrame of INFO annotations indexed to variant list.
        """
        info_table = pd.DataFrame([x.INFO for x in variants])
        if split:
            return self._split_info_table(info_table)

        else:
            return info_table

    def _split_info_table(self, info_table):
        """
        Split an INFO table into separate tables for site and allele level annotations.
        :param info_table:
        :return:
        """
        available = info_table.columns
        allele_fields = [x for x in self.info_allele_fields if x in available]
        site_fields = [x for x in self.info_value_fields + self.info_flag_fields if x in available]
        return info_table[site_fields], info_table[allele_fields]

    def split_variant(self, variant, alleles=[], exclude=None, value_only=False):
        """
        Split a multiallelic variant.

        :param variant:
        :param alleles: Alleles to be extracted, either by ALT allele or ALT allele index (default: all alleles).
        :param exclude: INFO fields to be dropped from the new variant records.
        :return:
        """
        if exclude is None:
            exclude = self.special_handling['INFO']

        def _extract_allele(variant, allele_index, exclude):
            new_variant = deepcopy(variant)
            # ALT field
            new_variant.ALT = [new_variant.ALT[allele_index], ]
            # INFO field
            info = new_variant.INFO
            [info.pop(x, 0) for x in exclude]  # Exclude ignored info fields
            if value_only:
                new_variant.INFO = {k: v[allele_index] if isinstance(v, list) else v for k, v in list(info.items())}
            else:
                new_variant.INFO = {k: [v[allele_index], ] if isinstance(v, list) else v for k, v in list(info.items())}
            return new_variant

        # Work out alleles to process, either all or selection by ALT allele base or index
        if alleles == ():
            alleles = list(range(len(variant.ALT)))
        else:
            if isinstance(alleles, int):
                pass
            if isinstance(alleles[0], str):  # Will work with single string or list
                try:
                    alleles = [variant.ALT.index(x) for x in alleles]
                except ValueError:
                    # Happens when VEP allele doesn't match variant
                    # Could occur because a multiallelic site contains an indel, MNP, and/or a SNP
                    trimmed_alts = [trim_common_suffix(str(x), variant.REF)[0] for x in variant.ALT]
                    alleles = [trimmed_alts.index(x) for x in alleles]

        # Extract desired alleles
        allelic_records = []
        for i in alleles:
            new_variant = _extract_allele(variant, i, exclude)
            allelic_records.append(new_variant)

        return allelic_records

    def vcf_row_to_table(self, variants, source_ids=None, include_other_info=False):
        """
        Parse PyVCF Variant to a DataFrame.


        :param variants:
        :return:
        """

        log.info('Processing {} variants...'.format(len(variants)))

        # Build source_id table
        if source_ids:
            source_id_series = pd.Series(source_ids)
            source_id_series.name = ('External', 'SOURCE_IDS')
            source_id_series.index.name = 'SITE'

        # Build record table
        log.info('Constructing site record table...')
        records = []
        for site, variant in enumerate(variants):
            for allele, alt in enumerate(variant.ALT):
                records.append([site, allele, variant.ID, variant.CHROM, variant.POS, variant.REF, alt, variant.QUAL,
                                variant.FILTER])
        row_record = pd.DataFrame(records, columns=['SITE', 'ALLELE_NUM', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'QUAL',
                                                    'FILTER'])
        row_record.set_index(['SITE', 'ALLELE_NUM'], inplace=True)

        # VEP table
        log.info('Tabulating VEP records...')
        vep_table = self.tabulate_variant_effects(variants)
        vep_table['ALLELE_NUM'] = vep_table['ALLELE_NUM'].astype(int)
        vep_table['ALLELE_NUM'] = vep_table['ALLELE_NUM'] - 1
        vep_table.index.name = 'SITE'
        vep_table.set_index('ALLELE_NUM', inplace=True, append=True)  # NB. Exclude 'Feature' to join on index later

        # INFO tables
        log.info('Processing INFO fields...')
        site_info, allele_info = self.tabulate_variant_info(variants, split=True)
        site_info.index.name = 'SITE'
        # Split allele INFO fields
        tables = []
        if self.info_allele_fields:
            for site_num, variant in enumerate(variants):
                info_dict = variant.INFO
                # FIXME: Need to handle absence of these types of fields
                allele_info_records = {k: info_dict[k] for k in self.info_allele_fields if k in info_dict}
                fields, both_allele_values = list(zip(*list(allele_info_records.items())))
                for allele_num, single_allele_values in enumerate(zip(*both_allele_values)):
                    allele_entry = dict(list(zip(fields, single_allele_values)))
                    allele_entry['ALLELE_NUM'] = allele_num
                    allele_entry['SITE'] = site_num
                    tables.append(allele_entry)
            split_allele_info = pd.DataFrame(tables).set_index(['SITE', 'ALLELE_NUM'])
            split_allele_info.columns = pd.MultiIndex.from_tuples([('Allele_INFO', x)
                                                                   for x in split_allele_info.columns],
                                                                  names=['Type', 'Field'])
        else:
            split_allele_info = pd.DataFrame()
        # Unbound info
        if include_other_info and self.info_other_fields:
            log.info('Processing "Other" INFO fields...')
            other_info = self.tabulate_variant_info(variants, split=False)
            available_info_other_fields = [x for x in self.info_other_fields if x in other_info.columns]
            other_info = other_info[available_info_other_fields]
            other_info.index.name = 'SITE'
            other_info.columns = pd.MultiIndex.from_tuples([('Other_INFO', x) for x in other_info.columns],
                                                           names=['Type', 'Field'])
        else:
            other_info = pd.DataFrame()

        # Create MultiIndex for sub-table columns
        row_record.columns = pd.MultiIndex.from_tuples([('Row', x) for x in row_record.columns],
                                                       names=['Type', 'Field'])
        vep_table.columns = pd.MultiIndex.from_tuples([('VEP', x) for x in vep_table.columns],
                                                      names=['Type', 'Field'])
        site_info.columns = pd.MultiIndex.from_tuples([('Site_INFO', x) for x in site_info.columns],
                                                      names=['Type', 'Field'])

        # Merge all the sub-tables
        log.info('Joining sub-tables...')
        # 1. Merged at site level
        merged_variant_table = row_record.join(site_info)
        if source_ids:
            merged_variant_table = merged_variant_table.join(source_id_series)
        if not other_info.empty:
            merged_variant_table = merged_variant_table.join(other_info)
        # 2. Merged at allele level
        if not split_allele_info.empty:
            merged_variant_table = merged_variant_table.join(split_allele_info)
        merged_variant_table = merged_variant_table.join(vep_table)
        # 3. Add Feature [and source_ids] to index
        log.info('Formatting result...')
        merged_variant_table.set_index(('VEP', 'Feature'), append=True, inplace=True)
        if not source_ids:
            merged_variant_table.index.set_names(['SITE', 'ALLELE_NUM', 'Feature'], inplace=True)
            # merged_variant_table = merged_variant_table.reorder_levels(['SITE', 'ALLELE_NUM', 'Feature'])
            # merged_variant_table.sort_index(inplace=True)
        else:
            # Add SOURCE_ID to the variant table index
            merged_variant_table.set_index(('External', 'SOURCE_IDS'), append=True, inplace=True)
            merged_variant_table.index.set_names(['SITE', 'ALLELE_NUM', 'Feature', 'SOURCE_ID'], inplace=True)
            merged_variant_table = merged_variant_table.reorder_levels(['SOURCE_ID', 'SITE', 'ALLELE_NUM', 'Feature'])
            merged_variant_table.sort_index(inplace=True)

        return merged_variant_table

    def get_gnomad_variants(self, aln_info_table, include_other_info=False):
        """

        :param aln_info_table:
        :param include_other_info:
        :return:
        """
        # NB. gnomad fetcher is packed into a generator, which is extracted in the following list comp.
        sequence_variant_lists = [(row.seq_id, (x for r in row.genomic_ranges for x in self.fetch(*r)))
                                  for row in aln_info_table.dropna(subset=['genomic_ranges']).itertuples()]
        all_variants = [(variant, seq_id)
                        for seq_id, range_reader in tqdm.tqdm(sequence_variant_lists, desc='Loading variants...')
                        for variant in range_reader]

        # pass VCF records and source_ids
        n = 1000  # chunking seems to interact with redundant rows... Fix by adding chunk ID with `keys`
        try:
            variants_table = pd.concat([self.vcf_row_to_table(*list(zip(*all_variants[i:i + n])),
                                                              include_other_info=include_other_info)
                                        for i in tqdm.tqdm(range(0, len(all_variants), n), desc='Parsing variants...')],
                                       keys=list(range(0, len(all_variants), n)))
        except ValueError:
            # probably no variants leading to no objects to concatenate
            return None
        
        # Write alignment variants to a VCF
        # TODO: add alignment to file name? (needs refactoring...)
        with open(os.path.join('results', 'alignment_variants.vcf'), 'w') as vcf_out:
            vcf_writer = vcf.Writer(vcf_out, self)
            for v, _ in all_variants:
                vcf_writer.write_record(v)

        return variants_table
