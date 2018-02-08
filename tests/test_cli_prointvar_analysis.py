from unittest import TestCase

from varalign.prointvar_analysis import cli


class TestCliAlignVariants(TestCase):
    """Test some aspects of the `prointvar_analysis` CLI."""

    def setUp(self):
        self.positionals = ['path_to_alignment']
        self.options = ['n_proc', 'override', 'only_sifts_best', 'max_pdbs']
        self.defaults = [1, False, False, None]

    def test_path_to_alignment(self):
        """Check first positional is called 'path_to_alignment'"""
        args = cli(['sample_swissprot_PF00001.18_full.sto'])
        self.assertEqual(args.path_to_alignment, 'sample_swissprot_PF00001.18_full.sto')

    def test_defaults(self):
        """Check defaults are as expected."""
        args = cli(['sample_swissprot_PF00001.18_full.sto'])  # dummy for positional
        mismatched = [k for k, v in zip(self.options, self.defaults) if not getattr(args, k) == v]
        message = 'Argument defaults are broken. Incorrect default(s) for {}'.format(mismatched)
        self.assertFalse(mismatched, message)

    def test_parsing(self):
        """Check we can modify parameters."""
        args = cli('--n_proc 2 --override sample_swissprot_PF00001.18_full.sto'.split())
        self.assertTrue(all([args.n_proc == 2, args.override]))
