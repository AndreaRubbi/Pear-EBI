__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from .emb_quality import DRM, pear_correlation


# ──────────────────────────────────────────────────────────────────────────────
# ─── t-SNE ND ─────────────────────────────────────────────────────────────────
def _perplexity_for(distance_matrix, preferred=3.0):
    """Return a legal perplexity for this many samples.

    sklearn requires perplexity < n_samples. The value was hardcoded at 3, so
    `--tsne` raised a raw ValueError for three trees or fewer. This clamps rather
    than rescales: whenever the original 3 is legal it is used unchanged, so
    embeddings of existing tree sets are numerically identical to before.
    """
    n_samples = len(distance_matrix)
    if n_samples <= 2:
        return 1.0
    return float(min(preferred, n_samples - 1))


def tsne(
    distance_matrix,
    n_dimensions,
    metadata=None,
    quality=False,
    report=False,
    output="./t-SNE_Embedding.csv",
):
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
        perplexity=_perplexity_for(distance_matrix),
    )

    Distances_embedded_ND = tsne.fit_transform(distance_matrix.astype(np.float32))
    pd.DataFrame(Distances_embedded_ND).to_csv(output, header=False, index=False)

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
