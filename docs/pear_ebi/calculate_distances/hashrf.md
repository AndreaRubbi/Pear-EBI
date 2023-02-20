#


### bash_command
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/hashrf.py/#L14)
```python
.bash_command(
   cmd
)
```


----


### hashrf
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/hashrf.py/#L37)
```python
.hashrf(
   file, n_trees, output_file
)
```

---
compute unweighted Robison Foulds distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (pandas.DataFrame) : computed distance matrix


----


### hashrf_weighted
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/hashrf.py/#L80)
```python
.hashrf_weighted(
   file, n_trees, output_file
)
```

---
compute weighted Robison Foulds distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (pandas.DataFrame) : computed distance matrix

