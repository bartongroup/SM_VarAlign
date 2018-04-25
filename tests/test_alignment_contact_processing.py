import os

from unittest import TestCase

import pandas as pd

from varalign import prointvar_analysis

root = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(root, 'data')


class Test_alignment_contact_processing(TestCase):
    def setUp(self):
        """Setup"""
        # Test contacts
        self.table = pd.read_pickle(os.path.join(data_path, '4iqr_residue_B146_structure_table.p.gz'))
        # Test info table
        self.info = pd.read_pickle(os.path.join(data_path, 'p41235_pf00104.29_alignment_info.p.gz'))
        # Mapping table
        self.index = pd.read_pickle(os.path.join(data_path, 'p41235_pf00104.29_mappings.p.gz'))
        self.mappings = prointvar_analysis._format_mapping_table(self.info, self.index)  # or read cached...

    def test_retain_relevant_contacts(self):
        """Test `_filter_extra_domain_contacts` does not exclude relevant entries"""
        filtered = prointvar_analysis._filter_extra_domain_contacts(self.table, self.info)
        # ResName/ResNum format e.g. ASN144
        self.assertIn('PHE286', filtered['label_comp_id_B'].str.cat(filtered['PDB_dbResNum_B']).tolist())
        self.assertIn('ALA278', filtered['label_comp_id_A'].str.cat(filtered['PDB_dbResNum_A']).tolist())

    def test_alignment_column_merge_inclusion(self):
        """Test `_merge_alignment_columns_to_contacts` does not exclude relevant entries"""
        # merge columns
        mapped = prointvar_analysis._merge_alignment_columns_to_contacts(self.mappings, self.table)

        # Test
        # ResName/ResNum format e.g. ASN144
        self.assertIn('PHE286', mapped['label_comp_id_B'].str.cat(mapped['PDB_dbResNum_B']).tolist())
        self.assertIn('ALA278', mapped['label_comp_id_A'].str.cat(mapped['PDB_dbResNum_A']).tolist())

    def test_ab_sort_inclusion(self):
        """Test `_sort_ab_contacts` does not exclude relevant contacts"""
        # sorting is done on alignment columns so need to merge columns first
        mapped = prointvar_analysis._merge_alignment_columns_to_contacts(self.mappings, self.table)

        # sort (under test)
        sorted_ = prointvar_analysis._sort_ab_contacts(mapped)

        # Test
        # ResName/ResNum format e.g. ASN144
        # NB. PHE286 is converted from ATOM_A to ATOM_B
        self.assertIn('PHE286', sorted_['label_comp_id_A'].str.cat(sorted_['PDB_dbResNum_A']).tolist())
        self.assertIn('ALA278', sorted_['label_comp_id_A'].str.cat(sorted_['PDB_dbResNum_A']).tolist())

    def test_ab_sort_preserves_index(self):
        """Test `_sort_ab_contacts` does not filter anything"""
        mapped = prointvar_analysis._merge_alignment_columns_to_contacts(self.mappings, self.table)
        sorted_ = prointvar_analysis._sort_ab_contacts(mapped)  # sort (under test)
        self.assertListEqual(sorted_.index.tolist(), mapped.index.tolist())
