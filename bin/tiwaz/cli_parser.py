import argparse

def default_cli_parser(parser_help):
    """Create an argparse object with the common flags."""

    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument('--run',
                        action='store_true',
                        help="Run the required simulations.")

    parser.add_argument('-f', '--force',
                        action='store_true',
                        help="Overwrite directories if need be.")

    parser.add_argument('--post-process',
                        action='store_true',
                        help="Perform post-processing.")

    parser.add_argument('--copy-to-paper',
                        action='store_true',
                        help="Copy all results from this run to the paper.")

    parser.add_argument('--reference',
                        action='store_true',
                        help="Also compute the reference solution.")

    return parser
