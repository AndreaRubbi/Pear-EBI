#


### quartet
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/tqdist.py/#L118)
```python
.quartet(
   file, n_trees, output_file
)
```

---
Computes quartet distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (numpy.ndarray) : computed distance matrix


**Raises**

* **PearExecutableError**  : if tqDist cannot be run, fails, or produces a matrix
    that does not match n_trees.


----


### triplet
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/tqdist.py/#L136)
```python
.triplet(
   file, n_trees, output_file
)
```

---
Computes triplet distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (numpy.ndarray) : computed distance matrix


**Raises**

* **PearExecutableError**  : if tqDist cannot be run, fails, or produces a matrix
    that does not match n_trees.

