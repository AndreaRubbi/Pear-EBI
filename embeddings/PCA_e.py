__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
# ──────────────────────────────────────────────────────────────────────────────
# ─── PCA N COMPONENTS ─────────────────────────────────────────────────────────
def pca(distance_matrix, n_components, metadata=None):
    pca = PCA(n_components=n_components)
    components = pca.fit_transform(distance_matrix)
    total_var = pca.explained_variance_ratio_.sum() * 100
    pd.DataFrame(components).to_csv('./PCA_Embedding.csv')

    return components