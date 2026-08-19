#


### hashrf
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/hashrf.py/#L183)
```python
.hashrf(
   file, n_trees, output_file
)
```

---
Computes unweighted Robinson Foulds distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (numpy.ndarray) : computed distance matrix


**Raises**

* **PearExecutableError**  : if hashrf cannot be run, or fails, or produces
    output that does not match n_trees.


----


### hashrf_weighted
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/hashrf.py/#L205)
```python
.hashrf_weighted(
   file, n_trees, output_file
)
```

---
Computes weighted Robinson Foulds distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (numpy.ndarray) : computed distance matrix


**Raises**

* **PearExecutableError**  : if hashrf cannot be run, or fails, or produces
    output that does not match n_trees.

