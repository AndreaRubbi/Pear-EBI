__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from sklearn.manifold import locally_linear_embedding

from .emb_quality import DRM, pear_correlation


# ──────────────────────────────────────────────────────────────────────────────
# ─── LLE N COMPONENTS ─────────────────────────────────────────────────────────
def lle(distance_matrix, n_components, metadata=None, quality=False, report=False):
    """embed distance_matrix in n_components with Locally Linear Embedding

    Args:
        distance_matrix (pandas.DataFrame): distance_matrix
        n_components (int): number of desired components
        metadata (pandas.DataFrame, optional): metadata of elements. Defaults to None.

    Returns:
        components (numpy.array): Embedding of distance matrix
    """
    embedding, _err_ = locally_linear_embedding(
        distance_matrix, n_neighbors=5, n_components=n_components
    )
    pd.DataFrame(embedding).to_csv("./ISOMAP_Embedding.csv")

    if report:
        Xr = None
        qu_re = DRM(distance_matrix, embedding, Xr)
    else:
        qu_re = None

    if quality:
        return (
            embedding,
            pear_correlation(distance_matrix, embedding),
            qu_re,
        )

    return embedding
