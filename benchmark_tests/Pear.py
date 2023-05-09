import os
import sys

from pear_ebi import tree_set

args = sys.argv
print(args[1], args[2])
set = tree_set.tree_set(args[1])
set.calculate_distances(args[2])
set.plot_2D("pca")
print("Done!")
