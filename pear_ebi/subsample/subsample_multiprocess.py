import inspect
import os
import random
import sys
from multiprocessing import Manager, Process

import numpy as np

try:
    from ..calculate_distances import maple_RF
except:
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)
    from calculate_distances import maple_RF

# TODO - I am trying to make this function recursive and
# TODO - running on multiple processes
#! PROBLEM: cannot share information between sub-processes
#! therefore it is hard to control the sampling


def subsample(file, n_trees, n_required, subp=True):
    """subsample a set of trees considering their
       distribution in the n_trees dimensional space.
       It tries to maximize the distance between the
       points in the sample considering the pairwise distance
       with respect to the furthest points found at a certain step.
       If the distance of a sample point P is not greater than
       the one between MD1 & MD2, then a random value is
       retrieved from a uniform distribution {0,1}. If the
       value is greater than 0.5, then the point is kept,
       else discarded. This allows to sample also considering
       the density of the points.

    Args:
        file (str): name of file containing the set of trees in newick format.
        n_trees (int): number of trees in set.
        n_required (int): number of trees in subsample.

    Returns:
        interesting points (list): list of trees subsampled.
        idxs (list): list of indexes of the trees subsampled.
    """
    global trees
    with open(file, "r") as f:
        trees = {tree: idx for idx, tree in enumerate(list(f.readlines()))}
        f.close()

    MD1_tree = random.sample(trees.keys(), 1)[0]
    MD1_idx = trees[MD1_tree]
    MD1 = maple_RF.readNewick(MD1_tree)[0]
    MD1_prep = maple_RF.prepareTreeComparison(MD1, rooted=False)
    trees.pop(MD1_tree)

    MD2_tree = random.sample(trees.keys(), 1)[0]
    MD2_idx = trees[MD2_tree]
    MD2 = maple_RF.readNewick(MD2_tree)[0]
    MD2_prep = maple_RF.prepareTreeComparison(MD2, rooted=False)
    trees.pop(MD2_tree)
    print(MD1)
    print(MD2)
    # with Manager() as manager:
    #    interesting_points = manager.list([MD1_tree, MD2_tree])
    #    idxs = manager.list([MD1_idx, MD2_idx])
    #    trees = manager.dict(trees)
    global interesting_points, idxs
    interesting_points = [MD1_tree, MD2_tree]
    idxs = [MD1_idx, MD2_idx]
    d_MD1_MD2 = maple_RF.RobinsonFouldsWithDay1985(MD2, MD1_prep, rooted=False)[0]

    def sample_points(MD1_prep, MD2_prep, d_MD1_MD2):
        while len(interesting_points) < n_required:
            P_tree = random.sample(trees.keys(), 1)[0]
            P_idx = trees[P_tree]
            P = maple_RF.readNewick(P_tree)[0]
            d_MD1_P = maple_RF.RobinsonFouldsWithDay1985(P, MD1_prep, rooted=False)[0]
            d_MD2_P = maple_RF.RobinsonFouldsWithDay1985(P, MD2_prep, rooted=False)[0]
            if d_MD1_P > d_MD1_MD2:
                d_MD1_MD2 = d_MD1_P
                trees.pop(P_tree)
                MD2, MD2_prep = P, maple_RF.prepareTreeComparison(P, rooted=False)
                interesting_points.append(P_tree)
                idxs.append(P_idx)
            elif d_MD2_P > d_MD1_MD2:
                d_MD1_MD2 = d_MD2_P
                trees.pop(P_tree)
                MD1, MD1_prep = P, maple_RF.prepareTreeComparison(P, rooted=False)
                interesting_points.append(P_tree)
                idxs.append(P_idx)
            elif random.randint(0, 100) > 50:  # alpha:
                trees.pop(P_tree)
                if random.randint(0, 10) > 5:
                    MD2_sub, MD2_prep_sub = P, maple_RF.prepareTreeComparison(
                        P, rooted=False
                    )
                    d_MD1_MD2_sub = d_MD1_P
                    args = (MD1_prep, MD2_prep_sub, d_MD1_MD2_sub)
                    sample_points(*args)
                    # p = Process(target = sample_points, args = args, daemon=False)
                    # p.start()
                    # p.join()
                else:
                    d_MD1_MD2_sub = d_MD2_P
                    MD1_sub, MD1_prep_sub = P, maple_RF.prepareTreeComparison(
                        P, rooted=False
                    )
                    args = (MD1_prep_sub, MD2_prep, d_MD1_MD2_sub)
                    sample_points(*args)
                    # p = Process(target = sample_points, args = args, daemon=False)
                    # p.start()
                    # p.join()

                interesting_points.append(P_tree)
                idxs.append(P_idx)
        return interesting_points, idxs

    return sample_points(MD1_prep, MD2_prep, d_MD1_MD2)


if __name__ == "__main__":
    args = sys.argv[1:]
    res = subsample(args[0], int(args[1]), int(args[2]))
    print("\n", res[0], "\n", res[1], "\n")
