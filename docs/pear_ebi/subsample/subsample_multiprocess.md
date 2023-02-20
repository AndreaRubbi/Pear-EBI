#


### subsample
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/subsample/subsample_multiprocess.py/#L23)
```python
.subsample(
   file, n_trees, n_required, subp = True
)
```

---
subsample a set of trees considering their
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


**Args**

* **file** (str) : name of file containing the set of trees in newick format.
* **n_trees** (int) : number of trees in set.
* **n_required** (int) : number of trees in subsample.


**Returns**

* **points** (list) : list of trees subsampled.
* **idxs** (list) : list of indexes of the trees subsampled.
