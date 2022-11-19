import argparse

def parser():
    parser = argparse.ArgumentParser(
        prog = 'TreeEmbedding',
        description='TreeEmbedding-EBI | \
        Calculates RF distances between large set of trees',
        epilog = 'Author: Andrea Rubbi - Goldman Group')
    parser.add_argument(type=str, dest='input', metavar = 'input',
                    help='input file : tree set in newic format')
    parser.add_argument('-o', type=str, dest='output', metavar = 'output',
                    help='output file : storage of distance matrix', required= False)
    parser.add_argument('--i', dest='interactive_mode',
                    help='run the program in interactive mode', 
                    required= False, action= 'store_true')
    parser.add_argument('-d', type=str, dest='distance_matrix', metavar = 'distance_matrix',
                    help='distance matrix : file of the distance matrix', required= False)
    parser.add_argument('-m', type=str, dest='metadata', metavar = 'metadata',
                    help='metadata : csv file with metadata for each tree', required= False)
    parser.add_argument('--hashrf', '--h', dest='hashrf', action= 'store_true',
                    help='calculates tree distances using hashrf', required= False)
    parser.add_argument('--pca', dest='pca', type= int,
                    help='embedding using pca: select 2 or 3 principal components;\n \
                        in case of 2pc, a 2D graph is produced; \n \
                            in case of #pc > 2, a 2D and a 3D graphs are produced', required= False)
    parser.add_argument('--tsne', dest='tsne', type= int,
                    help='embedding using tsne: select 2 or 3 final dimensions;\n \
                        in case of 2D, a 2D graph is produced; \n \
                            in case of #D > 2, a 2D and a 3D graphs are produced', required= False)
    
    return parser.parse_args()