import filecmp
import os
import shutil

from unittest import TestCase

import pandas as pd

from varalign.prointvar_analysis import main


class TestProintvarAnalysis(TestCase):

    @classmethod
    def setUpClass(cls):
        """Execute pipeline"""
        # Set up test directory
        start_dir = os.getcwd()
        test_dir = os.path.join(os.path.dirname(__file__), 'tmp')
        os.makedirs(test_dir)
        os.chdir(test_dir)

        # Execute pipeline
        # TODO: This file doesn't actually exist but the path is used as a prefix to find the expected output... OK?
        test_alignment = os.path.join(os.path.dirname(__file__), 'data', 'aligned_variants_test_expected',
                                      'sample_swissprot_PF00001.18_full.sto')
        main(path_to_alignment=test_alignment, override=True, only_sifts_best=1, max_pdbs=2, n_proc=2)

        cls.start_dir = start_dir
        cls.test_dir = test_dir

    def setUp(self):
        """Set up test environment"""

        # Expected output
        output_files = ['.structural.figures.pdf', '_column_data.csv', '_seq_structure_mappings.csv']
        output_files = ['sample_swissprot_PF00001.18_full.sto' + s for s in output_files]
        output_files = [os.path.join('results', f) for f in output_files]

        # Compare output with expected
        standard_path = os.path.join(os.path.dirname(__file__), 'data', 'prointvar_analysis_test_expected')
        comparison = filecmp.cmpfiles(standard_path, TestProintvarAnalysis.test_dir, output_files)

        self.output_files = output_files
        self.standard_path = standard_path
        self.comparison = comparison

    @classmethod
    def tearDownClass(cls):
        """Clean test environment"""

        os.chdir(TestProintvarAnalysis.start_dir)
        shutil.rmtree(TestProintvarAnalysis.test_dir)

    def test_expected_output_exists(self):
        cmpfiles_errors = self.comparison[2]
        message = 'The following file(s) appear to be missing: {}.'.format(cmpfiles_errors)
        self.assertFalse(cmpfiles_errors, message)

    def test_output_is_consistent(self):
        # FIXME: This test passes if the output doesn't exist!
        cmpfiles_mismatch = self.comparison[1]
        message = 'The following file(s) do not match their standards: {}'.format(cmpfiles_mismatch)
        self.assertFalse(cmpfiles_mismatch, message)

    def test_tabular_output_is_stable(self):
        # Load numeric data from CSV files
        csv_output = [f for f in self.output_files if f.endswith('.csv')]
        csv_tables_standard = [pd.read_csv(os.path.join(self.standard_path, f)) for f in csv_output]
        csv_tables_test = [pd.read_csv(os.path.join(TestProintvarAnalysis.test_dir, f)) for f in csv_output]
        # Identify any tables with columns that don't match
        mismatched_tables = []
        for f, standard_df, test_df in zip(csv_output, csv_tables_standard, csv_tables_test):
            # NB. need to dropna since nan != nan
            mismatched_columns = [col for col in standard_df
                                  if not (standard_df[col].dropna() == test_df[col].dropna()).all()]
            if mismatched_columns:
                mismatched_tables.append({f: mismatched_columns})
        # Test that we didn't find any mismatches or report where they were
        message = 'The following table(s) / column(s) do not match their standards: {}'.format(mismatched_tables)
        self.assertFalse(mismatched_tables, message)

    # TODO: Test output consistency for numeric data at a lower precision? (Implement when add rounding)
