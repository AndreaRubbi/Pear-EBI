__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import plotly.express as px
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from .graph import graph
# ──────────────────────────────────────────────────────────────────────────────
# ─── PCA N COMPONENTS ─────────────────────────────────────────────────────────
'''def pca(distance_matrix, n_components):
    pca = PCA(n_components=n_components)
    components = pca.fit_transform(distance_matrix)
    total_var = pca.explained_variance_ratio_.sum() * 100
    pd.DataFrame(components).to_csv('./PCA_Embedding.csv')

    if n_components >= 3:    
        fig = px.scatter_3d(
            components, x=0, y=1, z=2,
            title=f'Total Explained Variance: {total_var:.2f}%',
            labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
        )
        fig.write_html("./graph_PCA3.html")

    fig = px.scatter(components, x=0, y=1, 
            title=f'Total Explained Variance: {total_var:.2f}%') #color=Distances[WHAT]
    fig.write_html("./graph_PCA2.html")
    
    return components, fig'''
        
# ──────────────────────────────────────────────────────────────────────────────


def pca(distance_matrix, n_components, metadata=None):
    pca = PCA(n_components=n_components)
    components = pca.fit_transform(distance_matrix)
    total_var = pca.explained_variance_ratio_.sum() * 100
    pd.DataFrame(components).to_csv('./PCA_Embedding.csv')
    
    fig = graph.plot_embedding(components, metadata, n_components, name_plot='Tree_embedding')

    return components, fig