#


### plot_embedding
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/embeddings/graph/graph.py/#L17)
```python
.plot_embedding(
   data, metadata, dimensions, save = False, name_plot = 'Tree_embedding',
   static = False, plot_meta = 'SET-ID', plot_set = None, select = False,
   same_scale = False
)
```

---
Plot embedding of distance matrix - in 2D or 3D


**Args**

* **data** (pandas.DataFrame) : embedding of distance matrix
* **metadata** (pandas.DataFrame) : metadata of tree_set or set_collection
* **dimensions** (int) : number of dimensions - either 2 or 3
* **save** (bool, optional) : save plot in pdf format. Defaults to False.
* **name_plot** (str, optional) : name of plot. Defaults to 'Tree_embedding'.
* **static** (bool, optional) : if True, returns a less interactive format of plot. Defaults to False.
* **plot_meta** (str, optional) : defines the meta-feature used to color the points. Defaults to 'SET-ID'.
* **plot_set** (list or str, optional) : sets in set_collection to be plotted. Defaults to None.
* **select** (bool, optional) : if True, generates widgets that allow to show or hide uo to 16 set traces. Defaults to False.


**Returns**

* **image**  : plot with related widgets - interactive or static format

