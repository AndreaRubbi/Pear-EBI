import argparse


def parser():
    """Generates parser for PEAR program

    Returns:
        arg parser: PEAR parser
    """
    parser = argparse.ArgumentParser(
        prog="PEAR",
        description="PEAR-EBI | \
        Phylogeny Embedding and Approximate Representation \n \
        Calculates RF distances between large set of trees",
        epilog="Author: Andrea Rubbi - Goldman Group",
    )
    parser.add_argument(
        type=str,
        dest="input",
        metavar="input",
        help="input file : tree set in newic format",
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
        "-v",
        dest="verbose",
        action="store_true",
        help="print info while running",
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
        type=str,
        dest="distance_matrix",
        metavar="distance_matrix",
        help="distance matrix : file of the distance matrix",
        required=False,
    )
    parser.add_argument(
        "-m",
        type=str,
        dest="metadata",
        metavar="metadata",
        help="metadata : csv file with metadata for each tree",
        required=False,
    )
    parser.add_argument(
        "-hashrf",
        "--h",
        dest="hashrf",
        action="store_true",
        help="calculates tree distances using hashrf",
        required=False,
    )
    parser.add_argument(
        "-hashrf_w",
        "--hw",
        dest="hashrf_weighted",
        action="store_true",
        help="calculates tree distances using weighted hashrf",
        required=False,
    )
    parser.add_argument(
        "-daysRF",
        "--drf",
        dest="days_RF",
        action="store_true",
        help="calculates tree distances using an implementation of Day's algorithm",
        required=False,
    )
    parser.add_argument(
        "-pca",
        dest="pca",
        type=int,
        help="embedding using PCoA: select 2 or 3 principal components",
        required=False,
    )
    parser.add_argument(
        "-tsne",
        dest="tsne",
        type=int,
        help="embedding using t-SNE: select 2 or 3 final dimensions",
        required=False,
    )
    parser.add_argument(
        "-plot",
        dest="plot",
        type=int,
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

    return parser.parse_args()
