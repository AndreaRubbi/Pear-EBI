from pyTQDist import *

import unittest
import logging
import sys

class TestTQDist(unittest.TestCase):

    def testTripletDistance(self):
        log = logging.getLogger("python_test")
        log.debug("Testing tripletDistance: ")

        self.assertEqual(tripletDistance("../trees/tree_ab-c.new", "../trees/tree_ac-b.new"), 1)
        self.assertEqual(tripletDistance("../trees/tree_ac-b.new", "../trees/tree_ab-c.new"), 1)

        self.assertEqual(tripletDistance("../trees/test_tree1.new", "../trees/test_tree2.new"), 26)
        self.assertEqual(tripletDistance("../trees/test_tree2.new", "../trees/test_tree1.new"), 26)

        self.assertEqual(tripletDistance("../trees/test_tree3.new", "../trees/test_tree4.new"), 187793)
        self.assertEqual(tripletDistance("../trees/test_tree4.new", "../trees/test_tree3.new"), 187793)

        log.debug("ok.\n")
        

    def testAllPairsTripletDistance(self):
        log = logging.getLogger("python_test")
        log.debug("Testing allPairsTripletDistance: ")
        
        result = allPairsTripletDistance("../trees/two_trees.new")
        self.assertEqual(result[1][0], 26)
        self.assertEqual(result[0][0], 0 )
        self.assertEqual(result[2][1], 26)
        
        log.debug("ok.\n")

    def testPairsTripletDistance(self):
        log = logging.getLogger("python_test")
        log.debug("Testing pairsTripletDistance: ")

        result = pairsTripletDistance("../trees/pairs1.new", "../trees/pairs2.new")
        self.assertEqual(result[0], 26)
        self.assertEqual(result[1], 26)
        self.assertEqual(result[2], 187793)
        self.assertEqual(result[3], 187793)

        log.debug("ok.\n")

    def testQuartetDistance(self):
        log = logging.getLogger("python_test")
        log.debug("Testing quartetDistance: ")
        
        self.assertEqual(quartetDistance("../trees/quartet1.new", "../trees/quartet2.new"), 1)
        self.assertEqual(quartetDistance("../trees/quartet2.new", "../trees/quartet1.new"), 1)

        self.assertEqual(quartetDistance("../trees/test_tree1.new", "../trees/test_tree2.new"), 26)
        self.assertEqual(quartetDistance("../trees/test_tree2.new", "../trees/test_tree1.new"), 26)

        self.assertEqual(quartetDistance("../trees/test_tree3.new", "../trees/test_tree4.new"), 5485860)
        self.assertEqual(quartetDistance("../trees/test_tree4.new", "../trees/test_tree3.new"), 5485860)

        log.debug("ok.\n")

    def testPairsQuartetDistance(self):
        log = logging.getLogger("python_test")
        log.debug("Testing pairsQuartetDistance: ")

        result = pairsQuartetDistance("../trees/pairs1.new", "../trees/pairs2.new")
        self.assertEqual(result[0], 26)
        self.assertEqual(result[1], 26)
        self.assertEqual(result[2], 5485860)
        self.assertEqual(result[3], 5485860)

        log.debug("ok.\n")

    def testAllPairsQuartetDistance(self):
        log = logging.getLogger("python_test")
        log.debug("Testing allPairsQuartetDistance: ")
    
        result = allPairsQuartetDistance("../trees/five_trees.new")
        self.assertEqual(result[1][0], 1)
        self.assertEqual(result[0][0], 0)
        self.assertEqual(result[2][1], 1)
        self.assertEqual(result[3][2], 1)
        
        log.debug("ok.\n")

if __name__ == "__main__":
    logging.basicConfig( stream=sys.stderr )
    logging.getLogger( "python_test").setLevel( logging.DEBUG )
    unittest.main()
