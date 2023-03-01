__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from .emb_quality import DRM, pear_correlation


# ──────────────────────────────────────────────────────────────────────────────
# ─── t-SNE ND ─────────────────────────────────────────────────────────────────
def tsne(distance_matrix, n_dimensions, metadata=None, quality=False, report=False):
    """embed distance_matrix in n_components with t-Stochastic Neighbor Embedding

    Args:
        distance_matrix (pandas.DataFrame): distance_matrix
        n_dimensions (int): number of desired dimensions
        metadata (pandas.DataFrame, optional): metadata of elements. Defaults to None.

    Returns:
        Distances_embedded_ND (numpy.array): distances embedded in n_dimensions
    """
    method = "barnes_hut"
    if n_dimensions > 3:
        method = "exact"

    tsne = TSNE(
        n_components=n_dimensions,
        method=method,
        init="random",
        learning_rate=200.0,
        perplexity=3,
    )

    Distances_embedded_ND = tsne.fit_transform(distance_matrix.values.astype(np.float32))
    pd.DataFrame(Distances_embedded_ND).to_csv("./t-SNE_Embedding.csv")

    if report:
        # Xr = tsne.inverse_transform(Distances_embedded_ND)
        Xr = None
        qu_re = DRM(distance_matrix, Distances_embedded_ND, Xr)
    else:
        qu_re = None

    if quality:
        return (
            Distances_embedded_ND,
            pear_correlation(distance_matrix, Distances_embedded_ND),
            qu_re,
        )

    return Distances_embedded_ND
