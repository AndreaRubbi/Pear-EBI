import os, sys
import pandas as pd
import numpy as np
# getting the name of the directory
current = os.path.dirname(os.path.realpath(__file__))
 
# Getting the parent directory name
parent = os.path.dirname(current)
 
# adding the parent directory to
# the sys.path.
sys.path.append(parent)
from .embeddings import PCA_e, tSNE_e
from .calculate_distances import hashrf
from .interactive_mode import interactive
# ────────────────────────────────────────────────────────── TREE_SET CLASS ─────
class tree_set():
    # ─── INIT ──────────────────────────────────────────────────────────────────
    def __init__(self,
                file, 
                output_file=None,
                distance_matrix=None,
                metadata=None):
        
        self.file = file
        self.output_file = output_file
        self.distance_matrix = distance_matrix
        self.metadata = metadata 
        
        if self.output_file == None: self.output_file = "./{file}_distance_matrix.csv".format(file = os.path.splitext(os.path.basename(self.file))[0])
        
        if self.distance_matrix != None:
            try: pd.read_csv(self.distance_matrix) 
            except: print("There's an error with the Distance Matrix file - please check the correct location and name of the .csv file"), exit() 
              
        if self.metadata != None:
            try: pd.read_csv(self.metadata) 
            except: print("There's an error with the Metadata file - please check the correct location and name of the .csv file"), exit() 
        
        with open(file, 'r') as f:
            self.n_trees = len(f.readlines())
            f.close()
            
    # ─── STR ───────────────────────────────────────────────────────────────────
    def __str__(self):
        computed = 'not computed'
        if type(self.distance_matrix) != type(None): computed = 'computed'
            
        return(f"─────────────────────────────\n Tree set containing {self.n_trees} trees;\n File: {self.file};\n Distance matrix: {computed}.\n───────────────────────────── \n")
    
    # ─── CALCULATE DISTANCES ───────────────────────────────────────────────────
    def calculate_distances(self, method):
        methods = {'hashrf' : hashrf.hashrf,
                   'None' : None, }
        
        self.distance_matrix = methods[method](self.file, self.n_trees, self.output_file)
        print(f'{method} | Done!')
    
    # ─── EMBED ─────────────────────────────────────────────────────────────────
    def embed(self, method, dimensions):
        methods = {'pca' : PCA_e.pca,
                   'tsne' : tSNE_e.tsne,
                   'None' : None, }
        
        if dimensions == 1: raise Exception('Cannot embed in 1 dimension')
        if dimensions > self.n_trees: raise Exception('Cannot embed in #dimensions greater than #trees')
        
        if type(self.distance_matrix) == None: 
            raise Exception('Distance matrix has to be computed or inserted by the user prior the embedding')
        if method == 'pca':
            self.embedding_pca, fig = methods[method](self.distance_matrix, dimensions, self.metadata)
        if method == 'tsne':
            self.embedding_tsne, fig = methods[method](self.distance_matrix, dimensions, self.metadata)
        
        return fig
            
# ──────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────── SET_COLLECTION CLASS ─────
class set_collection(tree_set):
    def __init__(self, collection = list(), file = 'Set_collection', 
                 output_file = "./Set_collection_distance_matrix.csv"):
        self.file = file
        self.distance_matrix = None
        
        if file != 'Set_collection' and output_file == "./Set_collection_distance_matrix.csv": 
            self.output_file = "./{file}_distance_matrix.csv".format(file = os.path.splitext(os.path.basename(self.file))[0])
        else: self.output_file = output_file
        
        if isinstance(collection, tree_set): 
            self.collection = [collection]
            with open(self.file, 'w') as trees:
                with open(collection.file, 'r') as file:
                        trees.write(file.read())
                        file.close()
                trees.close()
            
        elif len(collection) > 0:
            for element in collection: assert isinstance(element, tree_set), 'Every element in a set_collection must be a tree_set'
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
            if metadata == None: metadata = pd.DataFrame()
            metadata['SET-ID'] = np.array([key] * set.n_trees)
            
            self.metadata = pd.concat([self.metadata, metadata])
            
            self.data[key] = {'metadata':metadata, 
                                 'n_trees':set.n_trees}
            
        
    
    def __add__(self, other):
        try: 
            assert isinstance(other, tree_set)
            return set_collection(self.collection + [other])
            
        except: 
            try: 
                assert isinstance(other, set_collection)
                return set_collection(self.collection + other.collection)
            except: 
                for element in other: assert isinstance(element, tree_set), 'You can concatenate a set_collection only with another set_collection, a tree_set, or a list of tree_set'
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
            assert isinstance(other, tree_set)
            return set_collection(self.collection + [other])
            
        except: 
            try: 
                assert isinstance(other, set_collection)
                return set_collection(self.collection + other.collection)
            except: 
                for element in other: assert isinstance(element, tree_set), 'You can concatenate a set_collection only with another set_collection, a tree_set, or a list of tree_set'
                return set_collection(self.collection + other)
            
            
            
        
        

         