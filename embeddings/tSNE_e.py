__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
# ──────────────────────────────────────────────────────────────────────────────
# ─── t-SNE ND ─────────────────────────────────────────────────────────────────
def tsne(distance_matrix, n_dimensions, metadata=None):
    """embed distance_matrix in n_components with t-Stochastic Neighbor Embedding

    Args:
        distance_matrix (pandas.DataFrame): distance_matrix
        n_dimensions (int): number of desired dimensions
        metadata (pandas.DataFrame, optional): metadata of elements. Defaults to None.

    Returns:
        Distances_embedded_ND (numpy.array): distances embedded in n_dimensions 
    """    
    method = 'barnes_hut'
    if n_dimensions > 3: method = 'exact'
    
    Distances_embedded_ND = TSNE(n_components=n_dimensions, method=method,
            init='random', learning_rate=200.0, perplexity=3).fit_transform(distance_matrix.values.astype(np.float32))
    pd.DataFrame(Distances_embedded_ND).to_csv('./t-SNE_Embedding.csv')

    return Distances_embedded_ND
