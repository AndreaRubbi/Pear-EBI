__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from sklearn.manifold import Isomap

from .emb_quality import DRM, pear_correlation


# ──────────────────────────────────────────────────────────────────────────────
# ─── ISOMAP N COMPONENTS ──────────────────────────────────────────────────────
def isomap(distance_matrix, n_components, metadata=None, quality=False, report=False):
    """embed distance_matrix in n_components with Isomap

    Args:
        distance_matrix (pandas.DataFrame): distance_matrix
        n_components (int): number of desired components
        metadata (pandas.DataFrame, optional): metadata of elements. Defaults to None.

    Returns:
        components (numpy.array): embedding of distance matrix
    """
    embedding = Isomap(n_components=n_components)
    components = embedding.fit_transform(distance_matrix)
    pd.DataFrame(components).to_csv("./ISOMAP_Embedding.csv")

    if report:
        Xr = None
        qu_re = DRM(distance_matrix, components, Xr)
    else:
        qu_re = None

    if quality:
        return (
            components,
            pear_correlation(distance_matrix, components),
            qu_re,
        )

    return components
