#


## Tree
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/maple_RF.py/#L47)
```python 
Tree(
   name = '', children = None, dist = 3.3e-05
)
```




**Methods:**


### .add_child
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/maple_RF.py/#L70)
```python
.add_child(
   node
)
```


----


### readNewick
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/maple_RF.py/#L76)
```python
.readNewick(
   tree_list, defaultBLen = 3.3e-05, normalizeInputBLen = 1.0
)
```

---
From a list of strings defining phylogenetic trees
in newick format, returns a list of Tree instances
containing the bipartitions of the tree in a
modularized form, suitable for comparisons with other trees.


**Args**

* **tree_list** (list) : list of trees in newick format.
* **defaultBLen** (float, optional) : default branch length. Defaults to 0.000033.
* **normalizeInputBLen** (float, optional) : value used to normalize branch lenghts. Defaults to 1.0.


**Returns**

* **trees**  : list of Tree instances


----


### prepareTreeComparison
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/maple_RF.py/#L177)
```python
.prepareTreeComparison(
   t1, rooted = False, minimumBLen = 6e-06
)
```

---
Prepares newick tree for comparison


**Args**

* **t1** (str) : input tree
* **rooted** (bool, optional) : set to True if t1 is rooted, default ro False
* **minimumBLen** (float, optional) : minimum value for branch length, default to 6E-6


**Returns**

* **metrics** (tuple of values) (leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths,) : tree metrics among which the RFL.


----


### RobinsonFouldsWithDay1985
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/maple_RF.py/#L353)
```python
.RobinsonFouldsWithDay1985(
   t2, t1, rooted = False, minimumBLen = 6e-06
)
```

---
Computes Robinson Foulds distances using Day's algorithm


**Args**

* **t2** (string) : newick tree to be compared to t1
* **t1** (tuple) : t1 after preprocessing using prepareTreeComparison()
* **rooted** (bool, optional) : True if t2 is rooted. Defaults to False.
* **minimumBLen** (float, optional) : minimum branch length. Defaults to 0.000006.


**Returns**

* **parameters** (tuple) (numDiffs, float(numDiffs) / (normalization), leafCount, foundBranches, missedBranches, (numBranches - foundBranches), RFL) : results of comparison, among which the RF distance (numDiffs).


----


### calculate_distance_matrix
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/maple_RF.py/#L515)
```python
.calculate_distance_matrix(
   file, n_trees, output_file
)
```

---
Computes the whole pipeline that calculates the pairwise distances in a collection of trees


**Args**

* **file** (str) : file containing the newick trees
* **n_trees** (int) : number of trees (or lines) in file
* **output_file** (str) : output file for distance matrix


**Returns**

* **distance_matrix** (np.array) : distance matrix

