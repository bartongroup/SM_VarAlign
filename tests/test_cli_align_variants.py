from unittest import TestCase

from varalign.cli import align_variants_parser


class TestCliAlignVariants(TestCase):
    """Test some aspects of the `align_variants` CLI."""

    def setUp(self):
        self.positionals = ['path_to_alignment']
        self.options = ['max_gaussians', 'n_groups', 'override', 'species']
        self.defaults = [5, 1, False, 'HUMAN']

    def test_path_to_alignment(self):
        """Check first positional is called 'path_to_alignment'"""
        args = align_variants_parser(['sample_swissprot_PF00001.18_full.sto'])
        self.assertEqual(args.path_to_alignment, 'sample_swissprot_PF00001.18_full.sto')

    def test_defaults(self):
        """Check defaults are as expected."""
        args = align_variants_parser(['sample_swissprot_PF00001.18_full.sto'])  # dummy for positional
        mismatched = [k for k, v in zip(self.options, self.defaults) if not getattr(args, k) == v]
        message = 'Argument defaults are broken. Incorrect default(s) for {}'.format(mismatched)
        self.assertFalse(mismatched, message)

    def test_parsing(self):
        """Check we can modify parameters."""
        args = align_variants_parser('--n_groups 2 --override sample_swissprot_PF00001.18_full.sto'.split())
        self.assertTrue(all([args.n_groups == 2, args.override]))
