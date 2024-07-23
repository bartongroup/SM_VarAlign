import os
from unittest import TestCase
from Bio import AlignIO
from varalign.core import aacon
from varalign.core.path_utils import get_issues_path

class TestGet_aacon(TestCase):
    def setUp(self):
        self.issues_path = get_issues_path()

    def test_aacon(self):
        # Reformat the test alignment
        alignment_path = os.path.join(self.issues_path, '006', 'peptide_fJ_with_ACCs.sto')
        alignment = AlignIO.read(alignment_path, 'stockholm')
        reformatted, _ = aacon._reformat_alignment_for_aacon(alignment)

        # Extract first sequence from each alignment
        s1 = str(alignment[0].seq).replace('-', '')
        s2 = str(reformatted[0].seq).replace('-', '')

        # Check they are equal
        message = 'Alignment reformatting has altered sequences.'
        self.assertEqual(s1, s2, message)
