import argparse

import pear_ebi


def parser():
    """Generates parser for PEAR program

    Returns:
        arg parser: PEAR parser
    """
    parser = argparse.ArgumentParser(
        prog="PEAR",
        description=f"PEAR-EBI v{pear_ebi.__version__} | \
        Phylogeny Embedding and Approximate Representation \n \
        Calculates Robison-Foulds distances between large set of trees",
        epilog="Author: Andrea Rubbi - Goldman Group | European Bioinformatics Institute",
    )

    parser.add_argument(
        "-v",
        dest="version",
        action="store_true",
        help="version of PEAR",
        required=False,
    )

    parser.add_argument(
        type=str,
        dest="input",
        metavar="input",
        help="input file : tree set in newic format",
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
        "--i",
        dest="interactive_mode",
        help="run the program in interactive mode",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--d",
        "-dM",
        type=str,
        dest="distance_matrix",
        metavar="distance_matrix",
        help="distance matrix : file of the distance matrix",
        required=False,
    )
    parser.add_argument(
        "-meta",
        type=str,
        dest="metadata",
        metavar="metadata",
        help="metadata : csv file with metadata for each tree",
        required=False,
    )
    parser.add_argument(
        "-m",
        "--m",
        dest="method",
        type=str,
        help="calculates tree distances using specified method (hashrf, weighted_hashrf, days_RF, quartet, triplet)",
        required=False,
    )
    parser.add_argument(
        "-pca",
        dest="pca",
        type=int,
        help="embedding using PCoA: select #principal components",
        required=False,
    )
    parser.add_argument(
        "-tsne",
        dest="tsne",
        type=int,
        help="embedding using t-SNE: select #final dimensions",
        required=False,
    )
    parser.add_argument(
        "-plot",
        "--p",
        dest="plot",
        action="store_true",
        help="plot embedding in 2 or 3 dimensions",
        required=False,
    )

    parser.add_argument(
        "-subset",
        "--s",
        dest="subset",
        type=int,
        help="extract subset of collection",
        required=False,
    )
    parser.add_argument(
        "-config",
        "--c",
        dest="config",
        type=str,
        help="toml config file",
        required=False,
    )
    parser.add_argument(
        "-report",
        "--r",
        action="store_true",
        dest="report",
        help="print long quality report of embedding",
        required=False,
    )
    parser.add_argument(
        "-quality",
        "--q",
        action="store_true",
        dest="quality",
        help="asess quality of embedding",
        required=False,
    )
    parser.add_argument(
        "-dir",
        dest="dir",
        type=str,
        help="directory with files",
        required=False,
    )
    parser.add_argument(
        "-pattern",
        dest="pattern",
        type=str,
        help="pattern of files in directory",
        required=False,
    )

    return parser.parse_args()
