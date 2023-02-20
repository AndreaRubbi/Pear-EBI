#


## Tree
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/maple_RF.py/#L25)
```python 
Tree(
   name = '', children = None, dist = 3.3e-05
)
```




**Methods:**


### .add_child
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/maple_RF.py/#L48)
```python
.add_child(
   node
)
```


----


### readNewick
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/maple_RF.py/#L54)
```python
.readNewick(
   tree_list, defaultBLen = 3.3e-05, normalizeInputBLen = 1.0
)
```

---
From a list of strings defining phylogenetic trees
in newick format, returns a list of Tree instances
containing the biartitions of the tree in a
modularized form, suitable for comparisons with other trees.


**Args**

* **tree_list** (list) : list of trees in newick format.
* **defaultBLen** (float, optional) : default branch length. Defaults to 0.000033.
* **normalizeInputBLen** (float, optional) : value used to normalize branch lenghts. Defaults to 1.0.


**Returns**

* **trees**  : list of Tree instances


----


### prepareTreeComparison
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/maple_RF.py/#L155)
```python
.prepareTreeComparison(
   t1, rooted = False, minimumBLen = 6e-06
)
```

---
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

---
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
making the RFL less accurate.

:param t1: input tree;
:param rooted: set to True if t1 is rooted, default ro False;
:param minimumBLen: minimum value for branch length, default to 6E-6;

:return: tree metrics, among which the RFL.

----


### RobinsonFouldsWithDay1985
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/maple_RF.py/#L327)
```python
.RobinsonFouldsWithDay1985(
   t2, t1, rooted = False, minimumBLen = 6e-06
)
```


----


### calculate_distance_matrix
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/maple_RF.py/#L478)
```python
.calculate_distance_matrix(
   file, n_trees, output_file
)
```

