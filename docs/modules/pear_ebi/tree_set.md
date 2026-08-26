#


## tree_set
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L433)
```python 
tree_set(
   file, output_file = None, distance_matrix = None, metadata = None
)
```


---
Class for the analysis of a set of phylogenetic trees


**Methods:**


### .tool_input
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L586)
```python
.tool_input()
```

---
Path to hand to the native tools.

The newline-normalised copy when there is one, otherwise self.file. Always
use this rather than self.file when invoking hashrf, tqDist or maple_RF:
tqDist and maple_RF read one tree per line and drop a final tree that has
no trailing newline.

There is no __del__ any more. It existed only to unlink the old
delete=False temp file, and __del__ is not guaranteed to run -- notably at
interpreter exit, and in notebooks where the object stays referenced.
TemporaryDirectory handles cleanup through weakref.finalize instead.

### .calculate_distances
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L602)
```python
.calculate_distances(
   method
)
```

---
Computes tree_set distance matrix with method of choice


**Args**

* **method** (str) : method/algorithm used to compute distance matrix


### .embed
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L625)
```python
.embed(
   method, dimensions, quality = False, report = False, output = None
)
```

---
Compute embedding with n-dimensions and method of choice


**Args**

* **method** (str) : method of choice to embed data
* **dimensions** (_type_) : number of dimensions/components
* **quality** (bool, optional) : returns quality report and self.emb_quality. Defaults to False.


### .plot_2D
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L768)
```python
.plot_2D(
   method, save = False, name_plot = None, static = False, plot_meta = 'SET-ID',
   plot_set = None, select = False, same_scale = False
)
```

---
Plot 2D embedding performed with method of choice


**Args**

* **method** (str) : embedding method
* **save** (bool, optional) : save plot HTML. Defaults to False.
* **name_plot** (str, optional) : name of plot's file. Defaults to None.
* **static** (bool, optional) : return less interactive plot. Defaults to False.
* **plot_meta** (str, optional) : meta-variable used to color the points. Defaults to "SET-ID".
* **plot_set** (list, optional) : list of sets to plot from set_collection. Defaults to None.
* **select** (bool, optional) : return set of buttons to show or hide specific traces. Defaults to False.
* **same_scale** (bool, optional) : use same color_scale for all traces when scale is continuous. Defaults to False.


**Raises**

* **ValueError**  : method can only be either pcoa or tsne for now


**Returns**

* **plot**  : either interactive or not


### .plot_3D
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L890)
```python
.plot_3D(
   method, save = False, name_plot = None, static = False, plot_meta = 'SET-ID',
   plot_set = None, select = False, same_scale = False, z_axis = None
)
```

---
Plot 3D embedding performed with method of choice


**Args**

* **method** (str) : embedding method
* **save** (bool, optional) : save plot HTML. Defaults to False.
* **name_plot** (str, optional) : name of plot's file. Defaults to None.
* **static** (bool, optional) : return less interactive plot. Defaults to False.
* **plot_meta** (str, optional) : meta-variable used to color the points. Defaults to "SET-ID".
* **plot_set** (list, optional) : list of sets to plot from set_collection. Defaults to None.
* **select** (bool, optional) : return set of buttons to show or hide specific traces. Defaults to False.
* **same_scale** (bool, optional) : use same color_scale for all traces when scale is continuous. Defaults to False.


**Raises**

* **ValueError**  : method can only be either pcoa or tsne for now


**Returns**

* **plot**  : either interactive or not


### .get_subset
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L1017)
```python
.get_subset(
   n_required, method = 'sequence'
)
```

---
Gets subset of phylogenetic trees


**Args**

* **n_required** (int) : number of points to extract
* **method** (str, optional) : method used to extract points ('sequence', 'random', 'syst'). Defaults to "sequence".


**Returns**

* **plots**  : 2D and 3D embedding plots of subset


----


## set_collection
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L1143)
```python 
set_collection(
   collection = list(), file = 'Set_collection_', output_file = None,
   distance_matrix = None, metadata = None
)
```




**Methods:**


### .calculate_distances
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L1288)
```python
.calculate_distances(
   method
)
```

---
Computes tree_set distance matrix with method of choice


**Args**

* **method** (str) : method/algorithm used to compute distance matrix


### .concatenate
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_set.py/#L1412)
```python
.concatenate(
   other
)
```

---
Concatenates two collections or collection and tree_set


**Args**

* **other** (tree_set ot set_collection) : tree_set ot set_collection


**Returns**

* **set_collection**  : concatenated set_collection

