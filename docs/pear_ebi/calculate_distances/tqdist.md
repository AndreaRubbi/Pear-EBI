#


### bash_command
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/tqdist.py/#L14)
```python
.bash_command(
   cmd
)
```


----


### quartet
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/tqdist.py/#L48)
```python
.quartet(
   file, n_trees, output_file
)
```

---
compute quartet distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (pandas.DataFrame) : computed distance matrix


----


### triplet
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/calculate_distances/tqdist.py/#L91)
```python
.triplet(
   file, n_trees, output_file
)
```

---
compute triplet distances


**Args**

* **file** (str) : name of input file with phylogenetic trees in newick format
* **n_trees** (int) : number of trees in file
* **output_file** (str) : name of output file that will contain the distance matrix


**Returns**

* **distance_matrix** (pandas.DataFrame) : computed distance matrix
