#


### tsne
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/embeddings/tSNE_e.py/#L26)
```python
.tsne(
   distance_matrix, n_dimensions, metadata = None, quality = False, report = False,
   output = './t-SNE_Embedding.csv'
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

