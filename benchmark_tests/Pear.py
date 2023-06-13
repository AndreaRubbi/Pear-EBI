import os
import sys

from pear_ebi import tree_set

args = sys.argv
print(args[1])
set = tree_set.tree_set(args[1])
set.calculate_distances("hashrf_RF")
set.plot_2D("pca")
print("Done!")
