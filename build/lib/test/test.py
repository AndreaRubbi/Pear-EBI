import os
import subprocess
import sys
import unittest
import warnings

import numpy as np
import pandas as pd

# getting the name of the directory
current = os.path.dirname(os.path.realpath(__file__))

# Getting the parent directory name
parent = os.path.dirname(current)

# adding the parent directory to
# the sys.path.
sys.path.append(parent)

from pear_ebi.tree_set import set_collection, tree_set

DIR = "../examples_tree_sets/beast_trees/"

FILE1 = "beast_run1.trees"
FILE2 = "beast_run2.trees"

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)


class TestPEAR(unittest.TestCase):
    def setUp(self):
        self.Set = tree_set(os.path.join(DIR, FILE1))
        self.Collection = set_collection([os.path.join(DIR, F) for F in (FILE1, FILE2)])

    def test_init(self):
        self.assertIsInstance(self.Set, tree_set)

        self.assertIsInstance(self.Collection, set_collection)

    def test_calculate_distances(self):
        methods = [
            "hashrf_RF",
            "hashrf_wRF",
            "smart_RF",
            "tqdist_quartet",
            "tqdist_triplet",
        ]

        self.Set.calculate_distances("hashrf_RF")
        dM = pd.read_csv("./beast_run1_distance_matrix.csv", header=None).values
        np.testing.assert_almost_equal(self.Set.distance_matrix, dM, decimal=1)
        cmds = ["rm -f ./beast_run1_distance_matrix.csv"]
        for cmd in cmds:
            subprocess.run(
                ["/bin/bash", "-c", cmd],
                stderr=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
            )

    def test_compute_embeddings(self):
        methods = ["pca", "tsne", "isomap", "lle"]
        self.Set.embed("pca", 3)
        embedding = pd.read_csv("./PCA_Embedding.csv", header=None).values
        np.testing.assert_almost_equal(self.Set.embedding_pca3D, embedding, decimal=1)
        cmds = ["rm -f ./beast_run1_distance_matrix.csv"]
        for cmd in cmds:
            subprocess.run(
                ["/bin/bash", "-c", cmd],
                stderr=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
            )

    def test_graph(self):
        self.Set.plot_2D("pca", static=True)


if __name__ == "__main__":
    import os
    import subprocess
    import sys
    import unittest
    import warnings

    import numpy as np
    import pandas as pd

    # getting the name of the directory
    current = os.path.dirname(os.path.realpath(__file__))

    # Getting the parent directory name
    parent = os.path.dirname(current)

    # adding the parent directory to
    # the sys.path.
    sys.path.append(parent)

    from pear_ebi.tree_set import set_collection, tree_set

    DIR = "../examples_tree_sets/beast_trees/"

    FILE1 = "beast_run1.trees"
    FILE2 = "beast_run2.trees"

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    unittest.main()
