#


### euclidean_distance
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/embeddings/emb_quality.py/#L5)
```python
.euclidean_distance(
   distances
)
```

---
Computes euclidean distances - norm of vectors


**Args**

* **distances** (np.array) : n-dimensional coordinates of points


**Returns**

* **euclidean_distances** (np.array) : array with euclidean distances


----


### pear_correlation
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/embeddings/emb_quality.py/#L18)
```python
.pear_correlation(
   distance_matrix, embedding
)
```

---
Computes pearson correlation between euclidean distances
in different dimensionalities


**Args**

* **distance_matrix** (np.array) : n-dimensional distance matrix
* **embedding** (np.array) : (m < n)-dimensional embedding


**Returns**

* **correlation** (float) : pearson correlation


----


### DRM
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/embeddings/emb_quality.py/#L34)
```python
.DRM(
   distance_matrix, embedding, inverse_emb
)
```

---
Return DRM object from https://github.com/zhangys11/pyDRMetrics


**Args**

* **distance_matrix** (np.array) : n-dimensional distance matrix
* **embedding** (np.array) : (m < n)-dimensional embedding
* **inverse_emb** (np.array) : reverse fit of model on embeddings


**Returns**

* **DRM**  : DRM object with quality metrics
