import numpy as np
import pandas as pd


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
    # Imported here rather than at module scope on purpose. pyDRMetrics 0.0.7 is
    # unmaintained and declares no dependencies, so an import failure in it used to
    # break `import pear_ebi.tree_set` outright -- every embedding module imports
    # this one. DRM() is the only thing that needs it, and it is reached only when
    # quality reporting is requested.
    from pyDRMetrics.pyDRMetrics import DRMetrics

    drm = DRMetrics(distance_matrix, embedding, inverse_emb)
    return drm
