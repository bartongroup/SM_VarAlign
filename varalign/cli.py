import argparse
import logging

log = logging.getLogger(__name__)
log.setLevel('INFO')


# Common parser
_parent_parser = argparse.ArgumentParser(add_help=False)
_parent_parser.add_argument('path_to_alignment', type=str, help='Path to the alignment.')
_parent_parser.add_argument('--override', help='Override any previously generated files.', action='store_true')

# align_variants parser
_align_variants_parser = argparse.ArgumentParser(add_help=False)
_align_variants_parser.add_argument('--max_gaussians', type=int, default=5,
                                    help='Maximum number of Gaussians for occupancy fitting.')
_align_variants_parser.add_argument('--n_groups', type=int, default=1, help='Top Gaussians to select after occupancy fitting.')
_align_variants_parser.add_argument('--species', type=str, help='Species (used for alignment filtering)', default='HUMAN')

# prointvar_analysis parser
_prointvar_analysis_parser = argparse.ArgumentParser(add_help=False)
_prointvar_analysis_parser.add_argument('--n_proc', type=int, help='Number of processors.', default=1)
_parser_n_sifts_group = _prointvar_analysis_parser.add_mutually_exclusive_group()
_parser_n_sifts_group.add_argument('--only_sifts_best', help='Process only sifts best structure.',
                                   action='store_true')
_parser_n_sifts_group.add_argument('--max_pdbs', type=int,
                                   help='Maximum number of SIFTs PDB mappings to use for a sequence')


def _parse_args_and_log(parser, argv, logger):
    """
    Parse arguments and log their values.

    :param parser: ArgumentParser.
    :param argv: None or list of arguments.
    :param logger: Logger.
    :return:
    """
    args = parser.parse_args(argv)
    if logger:
        # Log arguments
        for arg, value in sorted(vars(args).items()):
            logger.info("Command line argument %s: %r", arg, value)
        # TODO: or just `log.info(args)`
    return args


def align_variants_parser(argv=None, logger=log):
    parser = argparse.ArgumentParser(description='Align variants to a Pfam alignment.',
                                     parents=[_parent_parser, _align_variants_parser])
    return _parse_args_and_log(parser, argv, logger)


def prointvar_analysis_parser(argv=None, logger=log):
    parser = argparse.ArgumentParser(description='Structural properties of alignment columns.',
                                     parents=[_parent_parser, _prointvar_analysis_parser])
    return _parse_args_and_log(parser, argv, logger)


def varalign_parser(argv=None, logger=log):
    parser = argparse.ArgumentParser(description=('Analyse variants, conservation and structural features of a Pfam '
                                                  'alignment.'),
                                     parents=[_parent_parser, _align_variants_parser, _prointvar_analysis_parser])
    return _parse_args_and_log(parser, argv, logger)
