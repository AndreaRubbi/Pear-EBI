__author__ = "Andrea Rubbi"
__credits__ = {
    "Nicola de Maio": "EMBL-EBI",
    "Myrthe Willemsen": "Freie Universitat Berlin",
}
""" maple_RF is the adapted version of the implementation of the Day's algorithm
	for the Robison Foulds computation originally designed by Nicola de Maio for MAPLE.
	Further improvements thanks to Myrthe's work made the implementation more flexible.
	MAPLE code is available on GitHub: https://github.com/NicolaDM/MAPLE"""

import multiprocessing as mp
import os
import pickle
import shutil
import subprocess
import sys
import time
import warnings

import numpy as np
import pandas as pd
import rich
from rich.console import Console


class Tree(object):
    def __init__(self, name="", children=None, dist=0.000033):
        if name != "":
            self.name = name
        self.dist = dist
        self.replacements = 0
        self.children = list()
        self.up = None
        self.dirty = True

        if children is not None:
            for child in children:
                self.add_child(child)

    def __repr__(self):
        try:
            return str(self.name)
        except AttributeError:
            try:
                return self.name
            except AttributeError:
                return ""

    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)


# function to read input newick string
def readNewick(tree_list, defaultBLen=0.000033, normalizeInputBLen=1.0):
    """From a list of strings defining phylogenetic trees
    in newick format, returns a list of Tree instances
    containing the biartitions of the tree in a
    modularized form, suitable for comparisons with other trees.

    Args:
            tree_list (list): list of trees in newick format.
            defaultBLen (float, optional): default branch length. Defaults to 0.000033.
            normalizeInputBLen (float, optional): value used to normalize branch lenghts. Defaults to 1.0.

    Returns:
            trees: list of Tree instances
    """
    trees = list()

    tree_list = [tree_list] if type(tree_list) != type(list()) else tree_list
    for tree in tree_list:
        nwString = tree.replace("\n", "")
        branchLengthsIn = ":" in nwString
        index = 0
        node = Tree()
        name = distStr = ""
        finished = False
        while index < len(nwString):
            # start:bipartition
            if nwString[index] == "(":
                newNode = Tree()
                newNode.minorSequences = list()
                node.add_child(newNode)
                newNode.up = node
                node = newNode
                index += 1
            # end:tree
            elif nwString[index] == ";":
                trees.append(node)
                finished = True
                break
            # start
            elif nwString[index] == "[":
                # end
                while nwString[index] != "]":
                    index += 1
                index += 1
            # delimiter:branch_length
            elif nwString[index] == ":":
                index += 1
                while (
                    nwString[index] != ","
                    and nwString[index] != ")"
                    and nwString[index] != ";"
                ):
                    distStr += nwString[index]
                    index += 1
            # delimiter:nodes
            elif nwString[index] == ",":
                if name != "":
                    node.name, name = name, ""
                if distStr != "":
                    node.dist = float(distStr) * normalizeInputBLen
                    if node.dist < 0.0:
                        warnings.warn(
                            f"Warning: negative branch length in the input tree: {str(distStr)} ; converting it to positive."
                        )
                        node.dist = abs(node.dist)
                    distStr = ""
                else:
                    node.dist = defaultBLen

                newNode = Tree()
                newNode.minorSequences = list()
                node = node.up
                node.add_child(newNode)
                newNode.up = node
                node = newNode
                index += 1
            # end:bipartition
            elif nwString[index] == ")":
                if name != "":
                    node.name, name = name, ""
                if distStr != "":
                    node.dist = float(distStr) * normalizeInputBLen
                    distStr = ""
                else:
                    node.dist = defaultBLen
                node = node.up
                index += 1
            else:
                name += nwString[index]
                index += 1

        assert finished, "Error, final character ';' not found in newick tree."

    return trees


# Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
# prepareTreeComparison prepares the data to compare trees to a reference t1.
# RobinsonFouldsWithDay1985(), below, is the implementation of Day's algorithm.


def prepareTreeComparison(t1, rooted=False, minimumBLen=0.000006):
    """Prepares newick tree for comparison

    Args:
        t1 (str): input tree
        rooted (bool, optional): set to True if t1 is rooted, default ro False
        minimumBLen (float, optional): minimum value for branch length, default to 6E-6

    Returns:
        tree metrics (tuple of values) (leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths,): tree metrics among which the RFL.
    """
    """ Additional info from original program by NdM:
    For every branch that is shared among two trees,
    add the absolute difference of the two branches to the RFL.
    (This is first stored seperately as the KL value first).
    For every branch that is only in either of the two trees, add the length to the RFL distance.
    For the tree that is 'prepared for comparison', usually the true tree, two hash tables are stored
    (python dictionaries), containing all branch lengths:
            - One contains all leaf-branches. These are indexed by the leaf's name,
              usually a number (node.name).
            - The other contains all inner branches, indexed by the <Left, Right> values
              of the subtree that is attached by the branch of interest.

    Furthermore, it calculates the sum of all the inner branch lengths.
    Whenever in the second function of the Day's algorithm a branch is found in both trees,
    the absolute difference between the two branches is computed,
    and the result is the difference between that value and the difference with
    the branch length of the true tree.
    This calculation is not performed for leaf branches,
    as we can be sure that a leaf branch exists in both trees (if higher than the minimum length).

    Similar to the normal RF distance function, branches with length smaller
    than the minimum branch length are not taken into account.
    However, cases may occur that a branch exists in both trees,
    but in one of them it is lower than the minimum branch length,
    making the RFL less accurate."""

    # dictionary of values given to sequence names
    leafNameDict = dict()
    # list of sequence names sorted according to value
    leafNameDictReverse = list()
    # table containing clusters in the tree
    nodeTable = list()
    # dictionaries with branch lengths
    branchLengthDict, leafDistDict = dict(), dict()
    sumBranchLengths = 0  # sum of all branch lengths, except from the leaves
    # if comparing as unrooted trees, calculate tot num of leaves (using a postorder traversal), which will become useful later
    if not rooted:
        nLeaves = 0
        node = t1
        movingFrom = 0
        while node != t1.up:
            if (
                movingFrom == 0
            ):  # 0 means reaching node from parent, 1 means coming back from a child
                if len(node.children) == 0:
                    nLeaves += 1
                    nextNode = node.up
                    movingFrom = 1
                    nodeTable.append([0, 0])
                else:
                    nextNode = node.children[0]
                    movingFrom = 0
                    node.exploredChildren = 0
            else:
                nChildren = len(node.children)
                node.exploredChildren += 1
                if node.exploredChildren == nChildren:
                    nextNode = node.up
                    movingFrom = 1
                else:
                    nextNode = node.children[node.exploredChildren]
                    movingFrom = 0
            node = nextNode

    # implementing a non-recursive postorder traversal to assign values to internal nodes to fill nodeTable
    leafCount = 0
    node = t1
    movingFrom = 0
    lastL = float("inf")
    lastR = float("-inf")
    lastDesc = 0
    numBranches = 0
    while node != t1.up:
        if (
            movingFrom == 0
        ):  # 0 means reaching node from parent, 1 means coming back from a child
            if len(node.children) == 0:
                node.name = (node.name).replace("?", "_").replace("&", "_")
                leafNameDict[node.name] = leafCount
                leafNameDictReverse.append(node.name)
                if rooted:
                    nodeTable.append([0, 0])
                lastL = leafCount
                lastR = leafCount
                lastDesc = 1
                leafCount += 1
                nextNode = node.up
                movingFrom = 1
                leafDistDict[node.name] = node.dist
            else:
                node.exploredChildren = 0
                node.maxSoFar = float("-inf")
                node.minSoFar = float("inf")
                node.nDescendants = 0
                nextNode = node.children[0]
                movingFrom = 0
        else:
            nChildren = len(node.children)
            node.exploredChildren += 1
            if lastL < node.minSoFar:
                node.minSoFar = lastL
            if lastR > node.maxSoFar:
                node.maxSoFar = lastR
            node.nDescendants += lastDesc
            if node.exploredChildren == nChildren:
                nextNode = node.up
                movingFrom = 1
                lastL = node.minSoFar
                lastR = node.maxSoFar
                lastDesc = node.nDescendants
                if node != t1:
                    sumBranchLengths += node.dist
                if node == t1:
                    nodeTable[lastR][0] = lastL
                    nodeTable[lastR][1] = lastR
                else:
                    if (not rooted) and (node.up == t1) and (len(t1.children) == 2):
                        if node == t1.children[1]:
                            currentBL = node.dist + t1.children[0].dist
                            addBranch = True
                        else:
                            addBranch = False
                    else:
                        currentBL = node.dist
                        addBranch = True
                    if addBranch and currentBL > minimumBLen:
                        numBranches += 1
                        if rooted or lastL > 0:
                            if node == node.up.children[-1]:
                                if nodeTable[lastL][0] == 0 and nodeTable[lastL][1] == 0:
                                    nodeTable[lastL][0] = lastL
                                    nodeTable[lastL][1] = lastR
                                else:
                                    nodeTable[lastR][0] = lastL
                                    nodeTable[lastR][1] = lastR
                            else:
                                nodeTable[lastR][0] = lastL
                                nodeTable[lastR][1] = lastR
                            branchLengthDict[(lastL, lastR)] = currentBL
                        # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
                        else:
                            flippedL = lastR + 1
                            flippedR = nLeaves - 1
                            nodeTable[flippedL][0] = flippedL
                            nodeTable[flippedL][1] = flippedR
                            branchLengthDict[(flippedL, flippedR)] = currentBL
            else:
                nextNode = node.children[node.exploredChildren]
                movingFrom = 0
        node = nextNode
    return (
        leafNameDict,
        nodeTable,
        leafCount,
        numBranches,
        leafDistDict,
        branchLengthDict,
        sumBranchLengths,
    )


# Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
# this function compares the current tree t2 to a previous one for which prepareTreeComparison() was run.
# Example usage: leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(phyloTrue,rooted=False)
# numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(phyloEstimated,leafNameDict,nodeTable,leafCount,numBranches,rooted=False)
def RobinsonFouldsWithDay1985(t2, t1, rooted=False, minimumBLen=0.000006):
    """Computes Robison Foulds distances using Day's algorithm

    Args:
        t2 (string): newick tree to be compared to t1
        t1 (tuple): t1 after preprocessing using prepareTreeComparison()
        rooted (bool, optional): True if t2 is rooted. Defaults to False.
        minimumBLen (float, optional): minimum branch length. Defaults to 0.000006.

    Returns:
        comparison parameters (tuple) (numDiffs, float(numDiffs) / (normalization), leafCount, foundBranches, missedBranches, (numBranches - foundBranches), RFL): results of comparison, among which the RF distance (numDiffs).
    """
    (
        leafNameDict,
        nodeTable,
        leafCount,
        numBranches,
        leafDistDict,
        branchLengthDict,
        sumBranchLengths,
    ) = t1
    # implementing a non-recursive postorder traversal to check branch existance in the reference tree
    node = t2
    # branches in reference tree that are also in t2
    foundBranches = 0
    # branches in t2 that are not found in the reference
    missedBranches = 0
    movingFrom = 0
    lastL = float("inf")
    lastR = float("-inf")
    lastDesc = 0
    visitedLeaves = 0
    RFL = sumBranchLengths
    KF = 0
    while node != t2.up:
        if (
            movingFrom == 0
        ):  # 0 means reaching node from parent, 1 means coming back from a child
            if len(node.children) == 0:
                node.name = (node.name).replace("?", "_").replace("&", "_")
                if node.name in leafNameDict:
                    leafNum = leafNameDict[node.name]
                else:
                    print(node.name + " not in reference tree - aborting RF distance")
                    return None, None, None, None, None, None
                lastL = leafNum
                lastR = leafNum
                lastDesc = 1
                nextNode = node.up
                movingFrom = 1
                visitedLeaves += 1
                trueDist = leafDistDict[node.name]
                KF += abs(
                    trueDist - node.dist
                )  # As described, I have not added branch lengths of leaf nodes to sumBranchLengths, so no need for: RFL=- trueDist
            else:
                node.exploredChildren = 0
                node.maxSoFar = float("-inf")
                node.minSoFar = float("inf")
                node.nDescendants = 0
                nextNode = node.children[0]
                movingFrom = 0
        else:
            nChildren = len(node.children)
            node.exploredChildren += 1
            if lastL < node.minSoFar:
                node.minSoFar = lastL
            if lastR > node.maxSoFar:
                node.maxSoFar = lastR
            node.nDescendants += lastDesc
            if node.exploredChildren == nChildren:
                nextNode = node.up
                movingFrom = 1
                lastL = node.minSoFar
                lastR = node.maxSoFar
                lastDesc = node.nDescendants
                if node != t2:
                    if (not rooted) and (node.up == t2) and (len(t2.children) == 2):
                        if node == t2.children[1]:
                            currentBL = node.dist + t2.children[0].dist
                            searchBranch = True
                        else:
                            searchBranch = False
                    else:
                        currentBL = node.dist
                        searchBranch = True
                    if searchBranch and currentBL > minimumBLen:
                        if (lastR + 1 - lastL) == lastDesc:
                            if rooted or lastL > 0:
                                if (
                                    nodeTable[lastL][0] == lastL
                                    and nodeTable[lastL][1] == lastR
                                ):
                                    foundBranches += 1
                                    trueDist = branchLengthDict[(lastL, lastR)]
                                    KF += abs(trueDist - currentBL)
                                    RFL -= trueDist
                                elif (
                                    nodeTable[lastR][0] == lastL
                                    and nodeTable[lastR][1] == lastR
                                ):
                                    foundBranches += 1
                                    trueDist = branchLengthDict[(lastL, lastR)]
                                    KF += abs(trueDist - currentBL)
                                    RFL -= trueDist
                                else:
                                    missedBranches += 1
                                    RFL += currentBL
                            # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
                            else:
                                flippedL = lastR + 1
                                flippedR = leafCount - 1
                                if (
                                    nodeTable[flippedL][0] == flippedL
                                    and nodeTable[flippedL][1] == flippedR
                                ):
                                    foundBranches += 1
                                    trueDist = branchLengthDict[(flippedL, flippedR)]
                                    KF += abs(trueDist - currentBL)
                                    RFL -= trueDist
                                elif (
                                    nodeTable[flippedR][0] == flippedL
                                    and nodeTable[flippedR][1] == flippedR
                                ):
                                    foundBranches += 1
                                    trueDist = branchLengthDict[(flippedL, flippedR)]
                                    KF += abs(trueDist - currentBL)
                                    RFL -= trueDist
                                else:
                                    missedBranches += 1
                                    RFL += currentBL
                        else:
                            missedBranches += 1
                            RFL += currentBL
            else:
                nextNode = node.children[node.exploredChildren]
                movingFrom = 0
        node = nextNode
    if visitedLeaves < leafCount:
        print(
            f"There are leaves in the reference that have not been found in this new tree - leafCount {str(leafCount)} visitedLeaves {str(visitedLeaves)}"
        )
        return None, None, None, None, None, None
    # first value is number of differences, second value is max number of differences just in case one wants the normalized values;
    # the other values are there just in case one wants more detail.
    numDiffs = (numBranches - foundBranches) + missedBranches
    RFL += KF
    if rooted:
        normalization = numBranches + leafCount - 2
    else:
        normalization = numBranches + leafCount - 3
    return (
        numDiffs,
        float(numDiffs) / (normalization),
        leafCount,
        foundBranches,
        missedBranches,
        (numBranches - foundBranches),
        RFL,
    )


def calculate_distance_matrix(file, n_trees, output_file):
    """Computes the whole pipeline that calculates the pairwise distances in a collection of trees

    Args:
        file (str): file containing the newick trees
        n_trees (int): number of trees (or lines) in file
        output_file (str): output file for distance matrix

    Returns:
        distance_matrix (np.array): distance matrix
    """
    console = Console()
    global func_pool

    def func_pool(i):
        return subprocess.check_output(
            ["pypy3", f"{current}/RF_pypy.py", str(i)], universal_newlines=True
        )

    current = os.path.dirname(os.path.realpath(__file__))

    with open(file, "r") as f:
        trees = list(f.read().splitlines())
        f.close()
    distance_matrix = list()
    if shutil.which("pypy3") is not None:
        pckl = open("Trees.pckl", "wb")
        pickle.dump(trees, pckl)
        pckl.close()
        workers = os.cpu_count()
        if "sched_getaffinity" in dir(os):
            workers = len(os.sched_getaffinity(0))
        try:
            with mp.Pool(workers) as pool:
                arg_pool = list(range(len(trees)))
                results = list(pool.imap(func_pool, arg_pool))
            distance_matrix_upper = list(
                map(lambda res: eval(res.split("\n")[1].strip()), results)
            )

            distance_matrix = np.zeros((n_trees, n_trees))
            for i, line in enumerate(distance_matrix_upper):
                distance_matrix[i, i + 1 :] = line

            distance_matrix_lower = distance_matrix.transpose()
        except:
            print("Multiprocessing failed - trying on single core")
            print(
                "Suggestion: compute the distance matrix on a computing unit with more cores or analyze a smaller dataset"
            )

            distance_matrix = np.zeros((n_trees, n_trees))
            for i, tree in enumerate(trees):
                # status.update(f"[bold red]{i+1}/{len(trees)} [bold blue]Computing Robison Foulds distances...")
                inputTree = trees[i : i + 1]
                inputRFtrees = trees[i + 1 :]
                tree1 = readNewick(inputTree)[0]
                tree1_prep = prepareTreeComparison(tree1, rooted=False)
                otherTrees = readNewick(inputRFtrees)
                RF_distances = list()
                for tree in otherTrees:
                    res = RobinsonFouldsWithDay1985(tree, tree1_prep, rooted=False)
                    RF_distances.append(res[0])
                distance_matrix[i, i + 1 :] = RF_distances

            distance_matrix_lower = distance_matrix.transpose()

    else:
        distance_matrix = np.zeros((n_trees, n_trees))
        for i, tree in enumerate(trees):
            # status.update(f"[bold red]{i+1}/{len(trees)} [bold blue]Computing Robison Foulds distances...")
            inputTree = trees[i : i + 1]
            inputRFtrees = trees[i + 1 :]
            tree1 = readNewick(inputTree)[0]
            tree1_prep = prepareTreeComparison(tree1, rooted=False)
            otherTrees = readNewick(inputRFtrees)
            RF_distances = list()
            for tree in otherTrees:
                res = RobinsonFouldsWithDay1985(tree, tree1_prep, rooted=False)
                RF_distances.append(res[0])
            distance_matrix[i, i + 1 :] = RF_distances

        distance_matrix_lower = distance_matrix.transpose()

    distance_matrix = pd.DataFrame(distance_matrix + distance_matrix_lower)
    distance_matrix.to_csv(output_file, header=False, index=False)
    return distance_matrix.values
