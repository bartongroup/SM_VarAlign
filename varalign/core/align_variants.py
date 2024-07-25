import logging
import os

from Bio import AlignIO

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from numpy import vectorize

import pandas as pd

import seaborn as sns

import tqdm

from varalign.cli import cli
from varalign.config import defaults
from varalign.core import aacon
from varalign.core import alignments
from varalign.core import analysis_toolkit
from varalign.core import ensembl
from varalign.core import gnomad
from varalign.core import jalview
from varalign.core import occ_gmm
from varalign.core.utils import make_dir_if_needed

matplotlib.use("Agg")

# Constants
LOG_FILENAME = "align_variants.log"
RESULTS_PATH = "results"
DATA_PATH = os.path.join(".varalign", "aligned_variants_data")
DEFAULT_CANONICAL = eval(defaults.canonical)

# Set up logging
logging.basicConfig(
    filename=LOG_FILENAME,
    format="%(asctime)s %(name)s [%(levelname)-8s] - %(message)s",
    level=logging.INFO,
)
log = logging.getLogger(__name__)


def chunk_alignment(aln, n):
    """Return a generator that provides chunks of an alignment (sequence-wise)."""
    for i in range(0, len(aln), n):
        yield aln[i : i + n]


def chunk_table(table, n):
    """Return a generator that provides chunks of a DataFrame (row-wise)."""
    for i in range(0, len(table), n):
        yield table.iloc[i : i + n]


def save_table_and_log(method, path, description):
    """Save a table and log the action."""
    method(path)
    log.info(f"{description} saved to {path}")


def build_vep_filter(
    canonical=DEFAULT_CANONICAL,
    consequences=defaults.consequences,
    additional=defaults.additional,
):
    """Build a VEP filter string suitable for DataFrame.query()."""
    query = []
    if consequences != [""]:
        query.append(f"Consequence in {consequences}")
    if canonical:
        query.append('CANONICAL == "YES"')
    if additional != "":
        query.append(additional)
    return " & ".join(query)


def default_variant_filter(variants_table):
    """Apply standard filters to an alignment derived variant table."""
    filters = {
        "is_canonical": variants_table[("VEP", "CANONICAL")] == "YES",
        "is_ccds": variants_table[("VEP", "CCDS")] != "",
        "is_protein_coding": variants_table[("VEP", "BIOTYPE")]
        == "protein_coding",
        "at_protein_position": variants_table[("VEP", "Protein_position")]
        != "",
        "is_not_modifier": variants_table[("VEP", "IMPACT")] != "MODIFIER",
        "swissprot_matches_source": variants_table[
            "External", "SOURCE_ACCESSION"
        ]
        == variants_table["VEP", "SWISSPROT"],
        "trembl_matches_source": vectorize(lambda x, y: x in y)(
            variants_table[("External", "SOURCE_ACCESSION")],
            variants_table[("VEP", "TREMBL")],
        ),
    }
    filters["trembl_matches_source"][
        :
    ] = False  # OVERRIDE TREMBL TO KEEP ONLY SWISSPROT

    for name, condition in filters.items():
        log.info(f"{condition.sum()} variants pass {name} filter")

    combined_filter = (
        filters["is_canonical"]
        & filters["is_protein_coding"]
        & filters["is_not_modifier"]
        & filters["is_ccds"]
        & (
            filters["swissprot_matches_source"]
            | filters["trembl_matches_source"]
        )
        & filters["at_protein_position"]
    )

    return variants_table.loc[combined_filter].copy()


def map_uniprot_to_genome(uniprot, species="homo_sapiens", collapse=True):
    """Map a UniProt entry to the genome."""
    ensembl_genes = ensembl.get_xrefs(
        uniprot, species=species, features="gene"
    )
    if not ensembl_genes:
        return None

    ensembl_ranges = [ensembl.get_genomic_range(x) for x in ensembl_genes]
    for gene, gen_range in zip(ensembl_genes, ensembl_ranges):
        log.info(f"Mapped {uniprot} to {gene} on chr: {gen_range}")

    non_standard_ranges = [
        i
        for i, gen_range in enumerate(ensembl_ranges)
        if gen_range[0] not in ensembl.standard_regions
    ]
    for index in sorted(non_standard_ranges, reverse=True):
        del ensembl_genes[index]
        del ensembl_ranges[index]

    if not ensembl_ranges:
        log.warning(f"Could not map {uniprot} to the genome.")
        return None

    if collapse:
        ensembl_ranges = ensembl.merge_ranges(ensembl_ranges, min_gap=1000)

    return ensembl_ranges


def construct_mapping_table(alignment_info):
    """Construct a mapping table from alignment column to sequence residue."""
    mapping_table = pd.DataFrame(
        alignment_info["mapping"].tolist(),
        index=[alignment_info.index, alignment_info["seq_id"]],
    )
    mapping_table.reset_index(inplace=True)
    mapping_table = pd.melt(mapping_table, id_vars=["level_0", "seq_id"])
    mapping_table.dropna(subset=["value"], inplace=True)

    indexed_map_table = pd.DataFrame(
        mapping_table["value"].tolist(),
        columns=["Column", "Protein_position"],
        index=[mapping_table["seq_id"]],
    ).reset_index()
    indexed_map_table.set_index(["seq_id", "Protein_position"], inplace=True)
    indexed_map_table.index.rename(
        ["SOURCE_ID", "Protein_position"], inplace=True
    )
    indexed_map_table.columns = pd.MultiIndex.from_tuples(
        [("Alignment", col) for col in indexed_map_table.columns],
        names=["Type", "Field"],
    )
    return indexed_map_table


def calculate_column_occupancy(indexed_mapping_table):
    """Calculate column occupancy from a mapping table."""
    column_occupancy = (
        indexed_mapping_table[("Alignment", "Column")]
        .value_counts()
        .sort_index()
    )
    column_occupancy.name = "occupancy"
    return column_occupancy


def interpret_regression_results(
    regression_table, p_threshold=0.05, action=None
):
    """Provide human interpretation of regression results."""

    def create_result_string(control_name, pass_condition):
        return f'{control_name}... {"PASS" if pass_condition else "FAIL"}'

    negative_control_p = (
        regression_table.loc["filtered_synonymous", "pvalue"] > p_threshold
    )
    negative_control_m = (
        regression_table.loc["filtered_synonymous", "slope"] < 0
    )
    pass_negative = negative_control_p or negative_control_m

    positive_control_p = (
        regression_table.loc["filtered_missense", "pvalue"] < p_threshold
    )
    positive_control_m = regression_table.loc["filtered_missense", "slope"] > 0
    pass_positive = positive_control_p and positive_control_m

    results = [
        create_result_string(
            "Filtered synonymous vs. Shenkin (negative control)", pass_negative
        ),
        create_result_string(
            "Filtered missense vs Shenkin (positive control)", pass_positive
        ),
    ]

    if action:
        for result in results:
            action(result)
    else:
        return results


def write_variants_as_features(alignment_variant_table, feature_file_name):
    """Write a Jalview feature file marking up the variants in the alignment."""
    jalview.create_jalview_feature_file(
        {"missense_variant": "red", "synonymous_variant": "blue"},
        feature_file_name,
    )
    for (seq_id, consequence), variant_table in alignment_variant_table[
        "VEP"
    ].groupby(["SOURCE_ID", "Consequence"]):
        if consequence in ("missense_variant", "synonymous_variant"):
            residue_indexes = list(variant_table.index.get_level_values(1))
            variant_ids = list(variant_table["Existing_variation"])
            jalview.append_jalview_variant_features(
                seq_id.split("/")[0],
                residue_indexes,
                variant_ids,
                consequence,
                feature_file_name,
            )
    log.info(
        f"Wrote alignment variants to Jalview feature file {feature_file_name}"
    )


def get_genome_mappings(aln_info_table, species):
    """Map sequences to the genome and return genomic mappings."""
    genomic_ranges = [
        (row.seq_id, map_uniprot_to_genome(row.uniprot_id, species=species))
        for row in tqdm.tqdm(
            aln_info_table.itertuples(),
            total=len(aln_info_table),
            desc="Mapping sequences...",
        )
    ]
    if not genomic_ranges:
        log.error(
            "Failed to map any sequences to the genome... Are you sure there are human sequences?"
        )
        raise ValueError

    log.info(f"Mapped {len(genomic_ranges)} sequences to genome.")
    genomic_mapping_table = pd.DataFrame(
        genomic_ranges, columns=["seq_id", "genomic_ranges"]
    )
    return genomic_mapping_table


def map_variants_to_alignment(variants_df, residue_column_map):
    """Add alignment column numbers to a variant table."""
    variants_df.loc[:, ("VEP", "Protein_position")] = pd.to_numeric(
        variants_df.loc[:, ("VEP", "Protein_position")], errors="coerce"
    )
    variants_df.reset_index(["SITE", "ALLELE_NUM", "Feature"], inplace=True)
    variants_df.set_index(
        ("VEP", "Protein_position"), append=True, inplace=True
    )
    variants_df.index.set_names(
        ["SOURCE_ID", "Protein_position"], inplace=True
    )
    variants_df.sort_index(inplace=True)

    aligned_variants = residue_column_map.join(variants_df)
    aligned_variants.sort_index(inplace=True)
    return aligned_variants


def align_variants(
    aln_info_table, species="HUMAN", path_to_vcf=None, include_other_info=False
):
    """Align variants with the given alignment info table."""
    path_to_vcf = path_to_vcf or defaults.gnomad

    genomic_mapping_table = get_genome_mappings(aln_info_table, species)
    aln_info_table = aln_info_table.merge(
        genomic_mapping_table, on=["seq_id"], how="left"
    )

    parser = gnomad.Reader(
        filename=path_to_vcf, compressed=path_to_vcf.endswith("bgz") or None
    )
    variants_table = parser.get_gnomad_variants(
        aln_info_table, include_other_info=include_other_info
    )
    if variants_table.empty:
        log.warning("No variants found.")
        return variants_table

    source_uniprot_ids = aln_info_table.set_index("seq_id")["uniprot_id"]
    source_uniprot_ids.name = ("External", "SOURCE_ACCESSION")
    source_uniprot_ids.index.name = "SOURCE_ID"

    variants_table = variants_table.join(source_uniprot_ids)
    log.info(f"Variants before filtering:\t{len(variants_table)}")

    filtered_variants = default_variant_filter(variants_table)
    log.info(
        f'Redundant rows:\t{sum(filtered_variants.reset_index("Feature").index.duplicated())}'
    )
    filtered_variants.reset_index(level=0, drop=True, inplace=True)
    log.info(f"Total rows:\t{len(filtered_variants)}")

    indexed_map_table = construct_mapping_table(aln_info_table)
    aligned_variants = map_variants_to_alignment(
        filtered_variants, indexed_map_table
    )

    return aligned_variants


def run_aacon(alignment, results_prefix):
    """Run AACon and save the results."""
    conservation_methods = [x for x in aacon.aacon_methods if x != "LANDGRAF"]
    alignment_conservation = aacon.get_aacon(
        alignment, methods=conservation_methods
    )
    save_table_and_log(
        alignment_conservation.to_csv,
        f"{results_prefix}_aacon_scores.csv",
        "Formatted AACons results",
    )
    return alignment_conservation


def main(
    path_to_alignment,
    max_gaussians=5,
    n_groups=1,
    override=False,
    species="HUMAN",
):
    """Main function to align variants."""
    make_dir_if_needed(RESULTS_PATH)
    input_alignment_filename = os.path.basename(path_to_alignment)
    results_prefix = os.path.join(RESULTS_PATH, input_alignment_filename)
    make_dir_if_needed(DATA_PATH)
    data_prefix = os.path.join(DATA_PATH, input_alignment_filename)
    alignment = AlignIO.read(path_to_alignment, format="stockholm")

    is_data_available = all(
        [
            os.path.isfile(data_prefix + "_variants.p.gz"),
            os.path.isfile(data_prefix + "_info.p.gz"),
            os.path.isfile(data_prefix + "_mappings.p.gz"),
        ]
    )

    log.info("Generating alignment info table...")
    alignment_info = alignments.alignment_info_table(
        alignment, species
    )  # TODO: downstream this filters structural analysis too
    log.info(
        f"Alignment info table head:\n{alignment_info.head().to_string()}"
    )

    if override or not is_data_available:
        chunk_size = int(
            defaults.chunk_size
        )  # TODO: Optimise chunk size, consider N human sequences and other factors
        vartable_chunks = []
        for chunk in tqdm.tqdm(
            chunk_table(alignment_info, chunk_size),
            desc="Alignment chunks...",
            total=len(list(range(0, len(alignment_info), chunk_size))),
        ):
            try:
                _alignment_variant_table = align_variants(chunk)
            except AttributeError:
                continue
            vartable_chunks.append(_alignment_variant_table)
        alignment_variant_table = pd.concat(vartable_chunks)

        indexed_mapping_table = construct_mapping_table(alignment_info)
        save_table_and_log(
            alignment_info.to_pickle,
            data_prefix + "_info.p.gz",
            "Alignment info table pickle",
        )
        save_table_and_log(
            alignment_variant_table.to_pickle,
            data_prefix + "_variants.p.gz",
            "Alignment variant table pickle",
        )
        save_table_and_log(
            indexed_mapping_table.to_pickle,
            data_prefix + "_mappings.p.gz",
            "Alignment mapping table pickle",
        )
    else:
        log.info(f"Loading data for {path_to_alignment}...")
        alignment_info = pd.read_pickle(data_prefix + "_info.p.gz")
        alignment_variant_table = pd.read_pickle(
            data_prefix + "_variants.p.gz"
        )
        indexed_mapping_table = pd.read_pickle(data_prefix + "_mappings.p.gz")

    alignment_conservation = run_aacon(alignment, results_prefix)

    column_variant_counts = analysis_toolkit.count_column_variant_consequences(
        alignment_variant_table
    )
    save_table_and_log(
        column_variant_counts.to_csv,
        results_prefix + ".col_var_counts.csv",
        "Column variant counts",
    )

    rare_maf_threshold = 0.001
    is_rare = (
        alignment_variant_table[("Allele_INFO", "AF_POPMAX")]
        < rare_maf_threshold
    )
    column_rare_counts = analysis_toolkit.count_column_variant_consequences(
        alignment_variant_table[is_rare]
    )
    save_table_and_log(
        column_rare_counts.to_csv,
        results_prefix + ".col_rare_counts.csv",
        "Column rare variant counts",
    )

    is_missense = (
        alignment_variant_table[("VEP", "Consequence")] == "missense_variant"
    )
    column_missense_clinvar = analysis_toolkit.count_column_clinvar(
        alignment_variant_table[is_missense]
    )
    save_table_and_log(
        column_missense_clinvar.to_csv,
        results_prefix + ".col_mis_clinvar.csv",
        "Column missense variant ClinVar annotation frequencies",
    )

    is_synonymous = (
        alignment_variant_table[("VEP", "Consequence")] == "synonymous_variant"
    )
    column_synonymous_clinvar = analysis_toolkit.count_column_clinvar(
        alignment_variant_table[is_synonymous]
    )
    save_table_and_log(
        column_synonymous_clinvar.to_csv,
        results_prefix + ".col_syn_clinvar.csv",
        "Column synonymous variant ClinVar annotation frequencies",
    )

    column_occupancy = calculate_column_occupancy(
        indexed_mapping_table
    )  # TODO: consider adjusted count for unmapped sequences not covered in gnomAD
    column_summary = column_variant_counts.join(
        [column_missense_clinvar, column_occupancy, alignment_conservation]
    )

    gmms = occ_gmm._fit_mixture_models(
        column_summary["occupancy"], max_gaussians
    )
    M_best = occ_gmm._pick_best(gmms["models"], gmms["data"])
    subset_mask_gmm = occ_gmm._core_column_mask(M_best, gmms["data"], n_groups)
    column_summary = column_summary.assign(column_gmm_pass=subset_mask_gmm)

    variants_vs_occ = analysis_toolkit._comparative_regression(
        column_summary, "occupancy", filter_mask=subset_mask_gmm
    )
    save_table_and_log(
        variants_vs_occ.to_csv,
        results_prefix + ".variant_occ_regression.csv",
        "Variant vs. occupancy regression parameters",
    )

    shenkin_regressions = analysis_toolkit._comparative_regression(
        column_summary, "shenkin", filter_mask=subset_mask_gmm
    )
    save_table_and_log(
        shenkin_regressions.to_csv,
        results_prefix + ".variant_shenkin_regression.csv",
        "Variant vs. Shenkin regression parameters",
    )
    interpret_regression_results(shenkin_regressions, action=log.info)

    missense_scores = analysis_toolkit._column_variant_scores(
        column_summary[subset_mask_gmm],
        variant_class="missense_variant",
        occupancy="occupancy",
    )
    save_table_and_log(
        missense_scores.to_csv,
        results_prefix + ".col_missense_scores.csv",
        "Column missense scores",
    )
    column_summary = column_summary.join(missense_scores)
    column_summary = column_summary.join(
        column_summary.loc[subset_mask_gmm, "shenkin"].rank(pct=True),
        rsuffix="_percentile",
    )
    save_table_and_log(
        column_summary.to_csv,
        results_prefix + ".col_summary.csv",
        "Column summary data",
    )

    pdf = PdfPages(
        results_prefix + ".figures.pdf", metadata={"creationDate": None}
    )
    pdf.infodict().update(
        {
            "Title": f"Aligned Variant Diagnostics Plots for {path_to_alignment}",
            "Author": "align_variants.py",
        }
    )

    occ_gmm._gmm_plot(M_best, gmms["models"], gmms["data"])
    pdf.attach_note("Residue Occupancy GMM Diagnostics")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()

    fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
    column_summary.plot.scatter("occupancy", "missense_variant", ax=axs[0])
    column_summary.plot.scatter("occupancy", "synonymous_variant", ax=axs[1])
    axs[0].axvline(column_summary[subset_mask_gmm]["occupancy"].min())
    axs[1].axvline(column_summary[subset_mask_gmm]["occupancy"].min())
    fig.suptitle("N Variants vs. Occupancy")
    pdf.attach_note("N Variants vs. Occupancy")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()

    fig, axs = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True)
    sns.regplot(
        x="shenkin",
        y="missense_variant",
        data=column_summary[subset_mask_gmm],
        ax=axs[0],
    )
    pd.plotting.table(
        axs[0],
        shenkin_regressions.loc[["missense", "filtered_missense"]].round(2),
        loc="upper right",
        colWidths=[0.12] * 6,
        zorder=100,
    )
    sns.regplot(
        x="shenkin",
        y="synonymous_variant",
        data=column_summary[subset_mask_gmm],
        ax=axs[1],
    )
    pd.plotting.table(
        axs[1],
        shenkin_regressions.loc[["synonymous", "filtered_synonymous"]].round(
            2
        ),
        loc="upper right",
        colWidths=[0.12] * 6,
        zorder=100,
    )
    plt.title("N Variants vs. Shenkin")
    pdf.attach_note("N Variants vs. Shenkin")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()

    plot_data = column_summary[subset_mask_gmm].assign(
        pass_alpha=column_summary[subset_mask_gmm]["pvalue"] < 0.1
    )
    log.info(f"plot_data:\n{plot_data.head().to_string()}")

    protein_consequences = analysis_toolkit._aggregate_annotation(
        alignment_variant_table,
        ("VEP", "Consequence"),
        aggregate_by=["SOURCE_ID"],
    )
    protein_consequences.hist(facecolor="black", edgecolor="black")
    plt.title("Variants per Sequence")
    pdf.attach_note("Distribution of variants over alignment sequences")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()

    fig, axes = plt.subplots(1, 2)
    residue_counts = analysis_toolkit._aggregate_annotation(
        alignment_variant_table,
        ("VEP", "Consequence"),
        aggregate_by=["SOURCE_ID", "Protein_position"],
    )
    residue_counts = residue_counts.reindex(
        indexed_mapping_table.index
    ).fillna(0)
    residue_counts["missense_variant"].astype(int).value_counts().plot.bar(
        ax=axes[0], width=1, facecolor="black", edgecolor="black"
    )
    axes[0].set_title("Missense Variants per Residue")
    column_summary.loc[subset_mask_gmm, "missense_variant"].hist(
        ax=axes[1], facecolor="black", edgecolor="black"
    )
    axes[1].set_title("Missense Variants per Column")
    pdf.attach_note(
        "Distribution of variants over residues and alignment columns"
    )
    pdf.savefig(metadata={"creationDate": None})
    plt.close()
    pdf.close()

    umd_mask = column_summary.eval(
        "shenkin_percentile > 0.75 & oddsratio < 1 & pvalue < 0.1"
    )
    ume_mask = column_summary.eval(
        "shenkin_percentile > 0.75 & oddsratio > 1 & pvalue < 0.1"
    )
    cmd_mask = column_summary.eval(
        "shenkin_percentile < 0.25 & oddsratio < 1 & pvalue < 0.1"
    )
    cme_mask = column_summary.eval(
        "shenkin_percentile < 0.25 & oddsratio > 1 & pvalue < 0.1"
    )

    umd = umd_mask[umd_mask].index
    ume = ume_mask[ume_mask].index
    cmd = cmd_mask[cmd_mask].index
    cme = cme_mask[cme_mask].index

    indexed_mapping_table.reset_index().set_index(("Alignment", "Column")).loc[
        umd
    ].to_csv(results_prefix + ".umdres.csv")
    indexed_mapping_table.reset_index().set_index(("Alignment", "Column")).loc[
        ume
    ].to_csv(results_prefix + ".umeres.csv")
    indexed_mapping_table.reset_index().set_index(("Alignment", "Column")).loc[
        cmd
    ].to_csv(results_prefix + ".cmdres.csv")
    indexed_mapping_table.reset_index().set_index(("Alignment", "Column")).loc[
        cme
    ].to_csv(results_prefix + ".cmeres.csv")

    alignment_column_index = list(
        range(1, alignment.get_alignment_length() + 1)
    )
    jalview.marked_columns_track(
        umd_mask.reindex(alignment_column_index, fill_value=False),
        "UMD",
        "UMD columns at Shenkin PCR > 0.75 and missense OR < 1, p < 0.1",
        results_prefix + ".corners.ann",
    )
    jalview.marked_columns_track(
        ume_mask.reindex(alignment_column_index, fill_value=False),
        "UME",
        "UME columns at Shenkin PCR > 0.75 and missense OR > 1, p < 0.1",
        results_prefix + ".corners.ann",
        append=True,
    )
    jalview.marked_columns_track(
        cmd_mask.reindex(alignment_column_index, fill_value=False),
        "CMD",
        "CMD columns at Shenkin PCR < 0.25 and missense OR < 1, p < 0.1",
        results_prefix + ".corners.ann",
        append=True,
    )
    jalview.marked_columns_track(
        cme_mask.reindex(alignment_column_index, fill_value=False),
        "CME",
        "CME columns at Shenkin PCR < 0.25 and missense OR > 1, p < 0.1",
        results_prefix + ".corners.ann",
        append=True,
    )

    write_variants_as_features(
        alignment_variant_table, results_prefix + "_variant_features.feat"
    )

    log.info("DONE.")


if __name__ == "__main__":
    parameters = cli.align_variants_parser()
    main(**vars(parameters))
