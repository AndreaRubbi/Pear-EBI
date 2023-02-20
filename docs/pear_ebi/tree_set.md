#


## tree_set
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L68)
```python
tree_set(
   file, output_file = None, distance_matrix = None, metadata = None
)
```


---
Class for the analysis of a set of phylogenetic trees


**Methods:**


### .calculate_distances
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L144)
```python
.calculate_distances(
   method
)
```


### .embed
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L161)
```python
.embed(
   method, dimensions, quality = False
)
```


### .plot_2D
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L211)
```python
.plot_2D(
   method, save = False, name_plot = None, static = False, plot_meta = 'SET-ID',
   plot_set = None, select = False, same_scale = False
)
```


### .plot_3D
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L265)
```python
.plot_3D(
   method, save = False, name_plot = None, static = False, plot_meta = 'SET-ID',
   plot_set = None, select = False, same_scale = False
)
```


### .get_subset
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L321)
```python
.get_subset(
   n_required, method = 'sequence'
)
```


----


## set_collection
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L392)
```python
set_collection(
   collection = list(), file = 'Set_collection',
   output_file = './Set_collection_distance_matrix', distance_matrix = None
)
```




**Methods:**


### .concatenate
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/tree_set.py/#L510)
```python
.concatenate(
   other
)
```
