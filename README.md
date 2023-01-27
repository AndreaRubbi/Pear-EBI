# TreeEmbedding ![Logo Goldman](./goldman_logo.jpg | width=200)
<img src="./goldman_logo.jpg" width="200" height="200">
### by Goldman group - EMBL-EBI

TreeEmbedding can:
1. Compute the distance matrix given a set of phylogenetic trees;
2. Embed and represent the distance matrix in 2D or 3D.

## TreeEmbedding usage
Please, run `python3 ./main.py --help` to see the complete list of required arguments and flags.

### Simple usage
`python3 ./main.py examples_trees_sets/beast_trees_500_tips.txt --hashrf --pca 3`
this script calculates the distances between the trees in the file "beast_trees_500_tips.txt", which contains 2001 phylogenetic trees. 
The flag "--hashrf" defines the use of HasRF algorithm to perform this task [See HashRF documentation](https://code.google.com/archive/p/hashrf/).
The flag "--pca 3" indicates to the program to perform a Principal Coordinate Analysis to embed the distance matrix in 3 dimensions, returning the final coordinates in the lower-dimensional space and plotting them in 2D and 3D.

### Interactive mode
`python3 ./main.py examples_trees_sets/beast_trees_500_tips.txt --i`
this script launches the program in interactive mode. Once the program starts, it is going to guide you through its usage thanks to an intuitive interface.
