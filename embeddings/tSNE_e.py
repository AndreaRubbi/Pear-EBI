__author__ = "Andrea Rubbi"
# ──────────────────────────────────────────────────────────────────────────────
import plotly.express as px
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
# ──────────────────────────────────────────────────────────────────────────────
# ─── t-SNE ND ─────────────────────────────────────────────────────────────────
def tsne(distance_matrix, n_dimensions):
    method = 'barnes_hut'
    if n_dimensions > 3: method = 'exact'
    
    Distances_embedded_ND = TSNE(n_components=n_dimensions, method=method,
            init='random', perplexity=3).fit_transform(distance_matrix.values.astype(np.float32))
    pd.DataFrame(Distances_embedded_ND).to_csv('./t-SNE_Embedding.csv')
    
    if n_dimensions >= 3:    
        fig = px.scatter_3d(Distances_embedded_ND, x=0, y=1, z=2,
                 labels={'0': 'D1', '1': 'D2', '2': 'D3'},
                 title=f't-SNE 3D Embedding') #color=Distances[WHAT]
        fig.write_html("./graph_t-SNE3.html")
          
    fig = px.scatter(Distances_embedded_ND, x=0, y=1, 
                    title=f't-SNE 2D Embedding') #color=Distances[WHAT]
    fig.write_html("./graph_t-SNE2.html")
        
    return Distances_embedded_ND

# ──────────────────────────────────────────────────────────────────────────────

