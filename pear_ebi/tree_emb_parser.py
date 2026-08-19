import argparse

import pear_ebi


def parser(argv=None):
    """Parses PEAR's command line arguments.

    Args:
        argv (list, optional): argument list to parse. Defaults to None, which means
            sys.argv[1:], i.e. the real command line. Passing a list makes the parser
            testable and lets callers drive PEAR programmatically; without it there
            was no way to exercise this function at all.

    Returns:
        argparse.Namespace: the parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="PEAR",
        description=f"PEAR-EBI v{pear_ebi.__version__} | \
        Phylogeny Embedding and Approximate Representation",
        epilog="Author: Andrea Rubbi & others - Goldman Group | EMBL-European Bioinformatics Institute",
    )

    # argparse's own version action, replacing a commented-out -v block. A citable
    # scientific tool needs a machine-readable way to record which version produced a
    # result; there was none.
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"pear_ebi {pear_ebi.__version__}",
    )

    """parser.add_argument(
        "-v",
        dest="version",
        action="store_true",
        help="version of PEAR",
        required=False,
    )"""

    parser.add_argument(
        type=str,
        dest="input",
        metavar="input",
        help="input file : one or more tree sets in Newick format",
        nargs="*",  # "?",
    )
    parser.add_argument(
        "-o",
        type=str,
        dest="output",
        metavar="output",
        help="output file : storage of distance matrix",
        required=False,
    )
    parser.add_argument(
        "-i",
        "--interactive",
        dest="interactive_mode",
        help="run the program in interactive mode - only the input file, distance matrix, output file, and metadata arguments will be considered",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--dM",
        type=str,
        dest="distance_matrix",
        metavar="distance_matrix",
        help="distance matrix : specify file containing a precomputed distance matrix",
        required=False,
    )
    parser.add_argument(
        "--meta",
        type=str,
        dest="metadata",
        metavar="metadata",
        help="metadata : csv file with metadata for each tree",
        required=False,
    )
    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        type=str,
        help="calculates tree distances using specified method (hashrf_RF, hashrf_wRF, smart_RF, tqdist_quartet, tqdist_triplet)",
        required=False,
    )
    parser.add_argument(
        "--pcoa",
        dest="pcoa",
        metavar="N",
        type=int,
        help="embedding using PCoA: specify number of coordinate dimensions N (int) to be calculated",
        required=False,
    )
    parser.add_argument(
        "--tsne",
        dest="tsne",
        metavar="N",
        type=int,
        help="embedding using t-SNE: specify number of coordinate dimensions N (int) to be calculated",
        required=False,
    )
    parser.add_argument(
        "-p",
        "--plot",
        dest="plot",
        action="store_true",
        help="plot embedding in 2 or 3 dimensions",
        required=False,
    )

    """parser.add_argument(
        "-subset",
        "--s",
        dest="subset",
        type=int,
        help="extract subset of collection",
        required=False,
    )"""
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        type=str,
        help="toml config file",
        required=False,
    )
    """parser.add_argument(
        "-report",
        "--r",
        action="store_true",
        dest="report",
        help="print long quality report of embedding",
        required=False,
    )"""
    parser.add_argument(
        "-q",
        "--quality",
        action="store_true",
        dest="quality",
        help="assess quality of embedding",
        required=False,
    )
    parser.add_argument(
        "--dir",
        dest="dir",
        type=str,
        help="directory with files",
        required=False,
    )
    parser.add_argument(
        "--pattern",
        dest="pattern",
        type=str,
        help="pattern of files in directory",
        required=False,
    )

    return parser.parse_args(argv)
