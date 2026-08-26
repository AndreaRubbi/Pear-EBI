#


### lle
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/embeddings/LLE_e.py/#L12)
```python
.lle(
   distance_matrix, n_components, metadata = None, quality = False, report = False,
   output = './LLE_Embedding.csv'
)
```

---
embed distance_matrix in n_components with Locally Linear Embedding


**Args**

* **distance_matrix** (pandas.DataFrame) : distance_matrix
* **n_components** (int) : number of desired components
* **metadata** (pandas.DataFrame, optional) : metadata of elements. Defaults to None.


**Returns**

* **components** (numpy.array) : Embedding of distance matrix
