#


### tsne
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/master/pear_ebi/embeddings/tSNE_e.py/#L12)
```python
.tsne(
   distance_matrix, n_dimensions, metadata = None, quality = False
)
```

---
embed distance_matrix in n_components with t-Stochastic Neighbor Embedding


**Args**

* **distance_matrix** (pandas.DataFrame) : distance_matrix
* **n_dimensions** (int) : number of desired dimensions
* **metadata** (pandas.DataFrame, optional) : metadata of elements. Defaults to None.


**Returns**

* **Distances_embedded_ND** (numpy.array) : distances embedded in n_dimensions

