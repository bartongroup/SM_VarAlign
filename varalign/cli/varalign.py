#!/usr/bin/env python
import varalign.cli.cli
from varalign.core import align_variants
from varalign.core import prointvar_analysis


def _filter_args(argsd, *whitelist):
    return {k: argsd[k] for k in whitelist}


def main():
    argsd = vars(varalign.cli.cli.varalign_parser())
    align_variants.main(**_filter_args(argsd, 'path_to_alignment', 'max_gaussians', 'n_groups', 'override', 'species'))
    prointvar_analysis.main(**_filter_args(argsd, 'path_to_alignment', 'override', 'only_sifts_best', 'max_pdbs',
                                           'n_proc'))


if __name__ == '__main__':
    main()
