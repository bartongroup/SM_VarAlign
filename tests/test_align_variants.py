import filecmp
import logging
import os
import shutil

from unittest import TestCase, expectedFailure
from varalign.six import add_move, MovedModule; add_move(MovedModule('mock', 'mock', 'unittest.mock'))
from varalign.six.moves import mock

from varalign.align_variants import main
from varalign.config import defaults as mock_defaults

root = os.path.abspath(os.path.dirname(__file__))
mock_defaults.gnomad = "{}/data/sample_swissprot_PF00001.18_full.vcf.gz".format(root)


@mock.patch("varalign.config.defaults", mock_defaults)
class TestAlign_Variants(TestCase):
    """Test the align_variants pipeline."""

    @classmethod
    def setUpClass(cls):
        """Execute pipeline"""
        # Set up test directory
        start_dir = os.getcwd()
        test_dir = os.path.join(os.path.dirname(__file__), 'tmp')
        os.makedirs(test_dir)
        os.chdir(test_dir)

        # Execute pipeline
        test_alignment = os.path.join(os.path.dirname(__file__), 'data', 'sample_swissprot_PF00001.18_full.sto')
        main(path_to_alignment=test_alignment, max_gaussians=5, n_groups=1, override=True, species='HUMAN')

        cls.start_dir = start_dir
        cls.test_dir = test_dir

    def setUp(self):
        """Set up test environment"""

        # Expected output
        output_files = ['_full.sto.umeres.csv',
                        '_full.sto.cmdres.csv',
                        '_full.sto.col_mis_clinvar.csv',
                        '_full.sto.umdres.csv',
                        '_full.sto.variant_occ_regression.csv',
                        '_full.sto.cmeres.csv',
                        '_full.sto.col_rare_counts.csv',
                        '_full.sto.col_var_counts.csv',
                        '_full.sto.variant_shenkin_regression.csv',
                        '_full.sto_aacon_scores.csv',
                        '_full.sto.col_syn_clinvar.csv',
                        '_full.sto.corners.ann',
                        '_full.sto.col_missense_scores.csv',
                        '_full.sto.col_summary.csv',
                        '_full.sto.figures.pdf',
                        '_full.sto_variant_features.feat']
        output_files = ['sample_swissprot_PF00001.18' + s for s in output_files]
        output_files.append('alignment_variants.vcf')
        output_files = [os.path.join('results', f) for f in output_files]

        # Compare output with expected
        standard_path = os.path.join(os.path.dirname(__file__), 'data', 'aligned_variants_test_expected')
        comparison = filecmp.cmpfiles(standard_path, TestAlign_Variants.test_dir, output_files)

        self.output_files = output_files
        self.comparison = comparison

    @classmethod
    def tearDownClass(cls):
        """Clean test environment"""

        os.chdir(TestAlign_Variants.start_dir)
        shutil.rmtree(TestAlign_Variants.test_dir)

    def test_expected_output_exists(self):
        cmpfiles_errors = self.comparison[2]
        message = 'The following file(s) appear to be missing: {}.'.format(cmpfiles_errors)
        self.assertFalse(cmpfiles_errors, message)

    def test_output_is_consistent(self):
        cmpfiles_mismatch = [f for f in self.comparison[1] if not f.endswith('.pdf')]  # Exclude pdf
        message = 'The following file(s) do not match their standards: {}'.format(cmpfiles_mismatch)
        self.assertFalse(cmpfiles_mismatch, message)

    @expectedFailure
    def test_pdf_output_is_consistent(self):
        cmpfiles_mismatch = [f for f in self.comparison[1] if f.endswith('.pdf')]  # PDF only
        message = 'The following file(s) do not match their standards: {}'.format(cmpfiles_mismatch)
        self.assertFalse(cmpfiles_mismatch, message)
