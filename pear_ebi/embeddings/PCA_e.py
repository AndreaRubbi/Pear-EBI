__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
# ──────────────────────────────────────────────────────────────────────────────
# ─── PCA N COMPONENTS ─────────────────────────────────────────────────────────
def pca(distance_matrix, n_components, metadata=None):
    """embed distance_matrix in n_components with Principal Coordinate Analysis

    Args:
        distance_matrix (pandas.DataFrame): distance_matrix
        n_components (int): number of desired components
        metadata (pandas.DataFrame, optional): metadata of elements. Defaults to None.

    Returns:
        components (numpy.array): principal coordinates(components) of distance matrix 
    """
    pca = PCA(n_components=n_components)
    components = pca.fit_transform(distance_matrix)
    total_var = pca.explained_variance_ratio_.sum() * 100
    pd.DataFrame(components).to_csv('./PCA_Embedding.csv')

    return components