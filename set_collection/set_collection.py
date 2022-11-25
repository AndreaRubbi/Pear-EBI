import os
import pandas as pd
import numpy as np

from ..tree_set import tree_set


class set_collection(tree_set.tree_set):
    def __init__(self, collection = list(), file = 'Set_collection', 
                 output_file = "./Set_collection_distance_matrix.csv"):
        self.file = file
        self.distance_matrix = None
        
        if file != 'Set_collection' and output_file == "./Set_collection_distance_matrix.csv": 
            self.output_file = "./{file}_distance_matrix.csv".format(file = os.path.splitext(os.path.basename(self.file))[0])
        else: self.output_file = output_file
        
        if isinstance(collection, tree_set.tree_set): 
            self.collection = [collection]
            with open(self.file, 'w') as trees:
                with open(collection.file, 'r') as file:
                        trees.write(file.read())
                        file.close()
                trees.close()
            
        elif len(collection) > 0:
            for element in collection: assert isinstance(element, tree_set.tree_set), 'Every element in a set_collection must be a tree_set'
            self.collection = collection
            with open(self.file, 'w') as trees:
                for set in collection:
                    with open(set.file, 'r') as file:
                        trees.write(file.read())
                        file.close()
                trees.close()
        
        else: self.collection = collection
        
        with open(self.file, 'r') as f:
            self.n_trees = len(f.readlines())
            f.close()
            
        self.data = dict()
        
        self.metadata = pd.DataFrame()
        
        for set in self.collection:
            key = os.path.splitext(os.path.basename(set.file))[0]
            
            metadata = set.metadata
            metadata['SET-ID'] = np.array([key] * set.n_trees)
            
            self.metadata = pd.concat([self.metadata, metadata])
            
            self.data[key] = {'metadata':metadata, 
                                 'n_trees':set.n_trees}
            
        
    
    def __add__(self, other):
        try: 
            assert isinstance(other, tree_set.tree_set)
            return set_collection(self.collection + [other])
            
        except: 
            try: 
                assert isinstance(other, set_collection)
                return set_collection(self.collection + other.collection)
            except: 
                for element in other: assert isinstance(element, tree_set.tree_set), 'You can concatenate a set_collection only with another set_collection, a tree_set, or a list of tree_set'
                return set_collection(self.collection + other)
                
    def __str__(self):
        computed = 'not computed'
        if type(self.distance_matrix) != type(None): computed = 'computed'
            
        summary = f"─────────────────────────────\n Tree set collection containing {self.n_trees} trees;\n File: {self.file};\n Distance matrix: {computed}.\n───────────────────────────── \n"
        for key,value in self.data.items():
            summary += f"{key}; Containing {value['n_trees']} trees. \n" 
        
        return summary
        
    
    def concatenate(self, other):
        try: 
            assert isinstance(other, tree_set.tree_set)
            return set_collection(self.collection + [other])
            
        except: 
            try: 
                assert isinstance(other, set_collection)
                return set_collection(self.collection + other.collection)
            except: 
                for element in other: assert isinstance(element, tree_set.tree_set), 'You can concatenate a set_collection only with another set_collection, a tree_set, or a list of tree_set'
                return set_collection(self.collection + other)
            
            
            
        
        

         