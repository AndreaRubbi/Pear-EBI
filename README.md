# <font color='green'>P</font>hylogeny <font color='green'>E</font>mbedding & <font color='green'>A</font>pproximate <font color='green'>R</font>epresentation <img src="LOGO_PEAR.png" width="100" height="100">
### Goldman group - European Bioinformatics Institute <img src="goldman_logo.jpg" width="12" height="12">

PEAR can:
1. Compute the distance matrix given a set of phylogenetic trees;
2. Embed and represent the distance matrix in 2D or 3D.

## PEAR usage
Pear is both a python software and library. It can be installed with `python -m pip install pear_ebi` or downloaded from <a href="https://github.com/AndreaRubbi/Pear-EBI">Github</a>. Pear is currently compatible with Linux and Mac OSs. 

### PEAR as a python library
Once installed, Pear can be used to upload newick trees in python and represent them in embedded spaces. We recommend to use it on either jupyter notebook or lab, as these tools allow for more interaction with the graphs. On these platforms, the user is allowed to interact with widgets that allows to modify several parameteres of the plots. For specific uses and applications, see the <a href='https://github.com/AndreaRubbi/Pear-EBI/tree/pear_ebi/examples_tree_sets'>examples</a>.

### PEAR as a program
Run `python3 PEAR.py --help` to see the complete list of arguments and flags.
#### Simple usage
`python3 PEAR.py examples_trees_sets/beast_trees/beast_run1.trees -hashrf -pca 3 -plot 2`
this script calculates the unweighted <a href='https://doi.org/10.1016/0025-5564(81)90043-2'>Robison Foulds</a> distances between the trees in the file "beast_run1.trees", which contains 1001 phylogenetic trees. 
The flag "-hashrf" defines the use of [HasRF](https://code.google.com/archive/p/hashrf/) algorithm to perform this task.
The flag "-pca 3" indicates to the program to perform a Principal Coordinate Analysis to embed the distance matrix in 3 dimensions, returning the final coordinates in the lower-dimensional space and plotting them 3D thanks to "-plot 2".
#### Config file
A standard config toml file can be used for specific emebddings of multiple sets of trees. Instances of toml files are reported in ##EXAMPLE FOLDER##....
#### Interactive mode
`python3 PEAR.py examples_trees_sets/beast_trees/beast_run1.trees --i`
this script launches the program in interactive mode. Once the program starts, it is going to guide you through its usage thanks to an intuitive interface.


