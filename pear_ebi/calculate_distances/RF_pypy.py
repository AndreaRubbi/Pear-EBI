import sys

try:
    from .maple_RF import *
except ImportError:
    from maple_RF import *

if __name__ == "__main__":
    pckl = open("Trees.pckl", "rb")
    trees = pickle.load(pckl)
    pckl.close()

    i = int(sys.argv[1])
    inputTree = trees[i : i + 1]
    inputRFtrees = trees[i + 1 :]

    tree1 = readNewick(inputTree)[0]
    tree1_prep = prepareTreeComparison(tree1, rooted=False)
    otherTrees = readNewick(inputRFtrees)

    RF_distances = list()
    for tree in otherTrees:
        res = RobinsonFouldsWithDay1985(tree, tree1_prep, rooted=False)
        RF_distances.append(res[0])
    print("\n", RF_distances, "\n")
