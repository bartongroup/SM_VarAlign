import argparse
import logging

log = logging.getLogger(__name__)
log.setLevel('INFO')


_parent_parser = argparse.ArgumentParser(add_help=False)
_parent_parser.add_argument('path_to_alignment', type=str, help='Path to the alignment.')
_parent_parser.add_argument('--override', help='Override any previously generated files.', action='store_true')


def align_variants_parser(argv=None, logger=log):
    # TODO: Move this to appropriate script, its only still here whilst I update my usage on the cluster
    # CLI
    parser = argparse.ArgumentParser(description='Align variants to a Pfam alignment.', parents=[_parent_parser])
    parser.add_argument('--max_gaussians', type=int, default=5,
                        help='Maximum number of Gaussians for occupancy fitting.')
    parser.add_argument('--n_groups', type=int, default=1, help='Top Gaussians to select after occupancy fitting.')
    parser.add_argument('--species', type=str, help='Species (used for alignment filtering)', default='HUMAN')
    args = parser.parse_args(argv)

    if logger:
        # Log arguments
        for arg, value in sorted(vars(args).items()):
            logger.info("Command line argument %s: %r", arg, value)
        # TODO: or just `log.info(args)`

    return args


def prointvar_analysis_parser(argv=None, logger=log):
    # CLI
    parser = argparse.ArgumentParser(description='Structural properties of alignment columns.',
                                     parents=[_parent_parser])
    parser.add_argument('--n_proc', type=int, help='Number of processors.', default=1)
    parser_n_sifts_group = parser.add_mutually_exclusive_group()
    parser_n_sifts_group.add_argument('--only_sifts_best', help='Process only sifts best structure.',
                                      action='store_true')
    parser_n_sifts_group.add_argument('--max_pdbs', type=int,
                                      help='Maximum number of SIFTs PDB mappings to use for a sequence')
    args = parser.parse_args(argv)

    if logger:
        # Log arguments
        for arg, value in sorted(vars(args).items()):
            logger.info("Command line argument %s: %r", arg, value)
        # TODO: or just `log.info(args)`

    return args
