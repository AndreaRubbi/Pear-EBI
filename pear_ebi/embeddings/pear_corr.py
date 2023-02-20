import numpy as np
import pandas as pd


def euclidean_distance(distances):
    EC = np.array([np.linalg.norm(a - b) for a in distances for b in distances])
    return EC


def pear_correlation(distance_matrix, embedding):
    euclidean_dM = euclidean_distance(distance_matrix)
    euclidean_emb = euclidean_distance(embedding)
    return np.corrcoef(euclidean_dM, euclidean_emb)
