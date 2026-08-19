#


### pcoa
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/embeddings/PCoA_e.py/#L12)
```python
.pcoa(
   distance_matrix, n_components, metadata = None, quality = False, report = False,
   output = './PCoA_Embedding.csv'
)
```

---
embed distance_matrix in n_components with Principal Coordinate Analysis


**Args**

* **distance_matrix** (pandas.DataFrame) : distance_matrix
* **n_components** (int) : number of desired components
* **metadata** (pandas.DataFrame, optional) : metadata of elements. Defaults to None.


**Returns**

* **components** (numpy.array) : principal coordinates(components) of distance matrix

