import numpy as np
import pandas as pd
from pyDRMetrics.pyDRMetrics import DRMetrics


def euclidean_distance(distances):
    """Computes euclidean distances - norm of vectors

    Args:
        distances (np.array): n-dimensional coordinates of points

    Returns:
        euxlidean_distances (np.array): array with euclidean distances
    """
    EC = np.array([np.linalg.norm(a - b) for a in distances for b in distances])
    return EC


def pear_correlation(distance_matrix, embedding):
    """Computes pearson correlation between euclidean distances
    in different dimensionalities

    Args:
        distance_matrix (np.array): n-dimensional distance matrix
        embedding (np.array): (m < n)-dimensional embedding

    Returns:
        correlation (float): pearson correlation
    """
    euclidean_dM = euclidean_distance(distance_matrix)
    euclidean_emb = euclidean_distance(embedding)
    return np.corrcoef(euclidean_dM, euclidean_emb)


def DRM(distance_matrix, embedding, inverse_emb):
    """Return DRM object from https://github.com/zhangys11/pyDRMetrics

    Args:
        distance_matrix (np.array): n-dimensional distance matrix
        embedding (np.array): (m < n)-dimensional embedding
        inverse_emb (np.array): reverse fit of model on embeddings

    Returns:
        DRM: DRM object with quality metrics
    """
    drm = DRMetrics(distance_matrix, embedding, inverse_emb)
    return drm
