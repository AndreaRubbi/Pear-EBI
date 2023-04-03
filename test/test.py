import os
import sys
import unittest

import pear_ebi

DIR = "../examples_tree_sets/beast_trees/"

FILE1 = "beast_run1.trees"
FILE2 = "beast_run2.trees"


class TestPEAR(unittest.TestCase):
    def setUp(self):
        self.Set = pear_ebi.tree_set(os.path.join(DIR, FILE1))
        self.Collection = pear_ebi.set_collection(
            [os.path.join(DIR, F) for F in (FILE1, FILE2)]
        )

    def test_init(self):
        self.assertIsInstance(self.Set, pear_ebi.tree_set)

        self.assertIsInstance(self.Collection, pear_ebi.set_collection)

    def test_calculate_distances(self):
        methods = {
            "hashrf": pear_ebi.calculate_distances.hashrf.hashrf,
            "hashrf_weighted": pear_ebi.calculate_distances.hashrf.hashrf_weighted,
            "days_RF": pear_ebi.calculate_distances.maple_RF.calculate_distance_matrix,
            "quartet": pear_ebi.calculate_distances.tqdist.quartet,
            "triplet": pear_ebi.calculate_distances.tqdist.triplet,
            "None": None,
        }
        pass

    def test_compute_embeddings(self):
        pass

    pass


if __name__ == "__main__":
    unittest.main()
