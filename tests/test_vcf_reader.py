import os

from unittest import TestCase

from varalign import gnomad


class TestVCFReader(TestCase):
    """Test the PyVCF `vcf.Reader` subclass."""
    def test_clinvar_init(self):
        """Catch AttributeError during header parsing of ClinVar VCF (see issue #2)"""
        try:
            gnomad.Reader(filename=os.path.join(os.path.dirname(__file__), 'data', 'sample_clinvar_vep.vcf.gz'))
        except AttributeError:
            self.fail('gnomad.Reader.__init__ raised AttributeError for ClinVar VCF.')

    def test_clinvar_fetch(self):
        """Catch UnicodeDecodeError during fetching of ClinVar VCf"""
        parser = gnomad.Reader(filename=os.path.join(os.path.dirname(__file__), 'data', 'sample_clinvar_vep.vcf.gz'))
        try:
            [r for r in parser.fetch('1')]
        except UnicodeDecodeError:
            self.fail('`gnomad.Reader.fetch` raised UnicodeDecodeError whilst fetching variants from ClinVar VCF.')

    def test_default_encoding(self):
        """Ensure explicitly specified Reader.encoding is not overriden."""
        parser = gnomad.Reader(filename=os.path.join(os.path.dirname(__file__), 'data', 'sample_clinvar_vep.vcf.gz'))
        self.assertEqual(parser.encoding, 'utf_8')

    def test_encoding_can_be_specified(self):
        """Ensure explicitly specified Reader.encoding is not overriden."""
        parser = gnomad.Reader(filename=os.path.join(os.path.dirname(__file__), 'data', 'sample_clinvar_vep.vcf.gz'),
                               encoding='ascii')
        self.assertEqual(parser.encoding, 'ascii')
