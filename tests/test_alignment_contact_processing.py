import os

from unittest import TestCase

import pandas as pd

from varalign import prointvar_analysis

root = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(root, 'data')


class Test_alignment_contact_processing(TestCase):
    def test_retain_relevant_contacts(self):
        """Test `_filter_extra_domain_contacts` does not exclude relevant entries"""
        table = pd.read_pickle(os.path.join(data_path, '4iqr_residue_B146_structure_table.p.gz'))
        info = pd.read_pickle(os.path.join(data_path, 'p41235_pf00104.29_alignment_info.p.gz'))
        filtered = prointvar_analysis._filter_extra_domain_contacts(table, info)
        # ResName/ResNum format e.g. ASN144
        self.assertIn('PHE286', filtered['label_comp_id_B'].str.cat(filtered['PDB_dbResNum_B']).tolist())
        self.assertIn('ALA278', filtered['label_comp_id_A'].str.cat(filtered['PDB_dbResNum_A']).tolist())
