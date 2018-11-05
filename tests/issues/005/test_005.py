import os
import shutil
import unittest

import varalign.align_variants

from varalign.six import add_move, MovedModule; add_move(MovedModule('mock', 'mock', 'unittest.mock'))
from varalign.six.moves import mock


@mock.patch("varalign.align_variants.defaults.gnomad", os.path.join(os.path.dirname(__file__), 'ttc6_sub.sto_variants_sorted_normed_dedup.vcf.gz'))
class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original_dir = os.getcwd()
        os.chdir(os.path.dirname(__file__))

    @classmethod
    def tearDownClass(cls):
        """Clean test environment"""
        shutil.rmtree('results', ignore_errors=True)
        shutil.rmtree(os.path.join('.varalign', 'aacon'), ignore_errors=True)
        shutil.rmtree(os.path.join('.varalign', 'aligned_variants_data'), ignore_errors=True)
        os.chdir(cls.original_dir)

    def test_something(self):
        test_alignment = os.path.join(os.path.dirname(__file__), 'ttc6_sub.sto')
        try:
            varalign.align_variants.main(path_to_alignment=test_alignment, max_gaussians=5, n_groups=1, override=True,
                                         species='HUMAN')
        except ValueError:
            self.fail('align_variants.main failed with value error.')


if __name__ == '__main__':
    unittest.main()
