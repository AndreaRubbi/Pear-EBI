__author__ = 'Andrea Rubbi'

# ─── About ────────────────────────────────────────────────────────────────────
'''This file defines the tree_set and set_collection classes.
    tree_set contains the information relative to a single set of phyogenetic trees
    in newick format. It allows to compute the distance matrix using different methods and metrics.
    The distance matrix can be then embedded using different methods and subsequently plotted in 2D or 3D.
    A distance matrix and metadata can be given as .csv files. Moreover, metadata is modified
    in order to give information regarding the name of the tree set and the index (or step) of each tree.
    Please note that, once an instance of a class is generated, its metadata dataframe should not be substituted
    as this would invalidate it for the plotting functions. Addition of columns and features is possible by
    accessing the dataframe and modifying it as a pandas.DataFrame instance. 
    set_collection behaves similarly to set_collection. Matter of fact, it is a subclass of the latter and therefore
    shares most of its methods. Its purpose is to analyze concurrently multiple instances of tree_sets and plot their 
    relative distance in a common embedding. Examples of possible applications are present at: ###LINK###'''
# ──────────────────────────────────────────────────────────────────────────────

__copyright__ = "..."
__credits__ = ["Andrea Rubbi", "Lukas Weilguny", "Nick Goldman", "Nicola de Maio"]

__license__ = "..."
__version__ = "1.0.1"
__maintainer__ = "Andrea Rubbi"
__institute__ = 'EMBL-EBI'
__email__ = 'andrear@ebi.ac.uk'
__status__ = "Production"

# ──────────────────────────────────────────────────────────────────────────────
import os, sys
import pandas as pd
import numpy as np
import uuid

#? rich is a very nice library that allows to 
#? easily format the output of console
#? https://github.com/Textualize/rich
from rich.console import Console
from rich import print

# getting the name of the directory
current = os.path.dirname(os.path.realpath(__file__))
 
# Getting the parent directory name
parent = os.path.dirname(current)
 
# adding the parent directory to
# the sys.path.
sys.path.append(parent)

# importing other modules
from .embeddings import PCA_e, tSNE_e
from .calculate_distances import hashrf
from .interactive_mode import interactive
from .embeddings.graph import graph

# ────────────────────────────────────────────────────────── TREE_SET CLASS ─────
class tree_set():
    """Class for the analysis of a set of phylogenetic trees"""
    # Console from rich -> takes control of console output
    console = Console()
    # ─── INIT ──────────────────────────────────────────────────────────────────
    def __init__(self,
                file, 
                output_file=None,
                distance_matrix=None,
                metadata=None):
        """ Initialize tree_set
        
        file: mandatory - file with set of phylogenetic trees in newick format
        output_file: facultative - specifies output_file of distance matrix
        distance_matrix: facultative - specifies file with already-computed distance matrix
        metadata: dìfacultative - specifies file containing additional information for each tree in set.
        It should contain a column for each feature, a row for each tree (blank row if no info) """
        
        self.file = file
        self.output_file = output_file
        self.distance_matrix = distance_matrix
        self.metadata = metadata
        self.embedding_pca2D = None
        self.embedding_tsne2D = None
        self.embedding_pca3D = None
        self.embedding_tsne3D = None
                 
        if self.output_file == None: self.output_file = "./{file}_distance_matrix.csv".format(file = os.path.splitext(os.path.basename(self.file))[0])
        
        with open(file, 'r') as f:
            self.n_trees = len(f.readlines())
            f.close()
        
        if type(self.distance_matrix) != type(None):
            try: pd.read_csv(self.distance_matrix) 
            except: print("There's an error with the Distance Matrix file - please check the correct location and name of the .csv file"), exit() 
              
        if type(self.metadata) != type(None):
            try: self.metadata = pd.read_csv(self.metadata) 
            except: print("There's an error with the Metadata file - please check the correct location and name of the .csv file"), exit() 
        
        else: self.metadata = pd.DataFrame()
        self.metadata['SET-ID'] = [os.path.splitext(os.path.basename(self.file))[0] for i in range(self.n_trees)]
        self.metadata['STEP'] = [i for i in range(self.n_trees)]
        self.sets = np.unique(self.metadata['SET-ID'])
            
    # ─── STR ───────────────────────────────────────────────────────────────────
    def __str__(self):
        computed = 'not computed'
        if type(self.distance_matrix) != type(None): computed = 'computed'
            
        return(f"─────────────────────────────\n Tree set containing {self.n_trees} trees;\n File: {self.file};\n Distance matrix: {computed}.\n───────────────────────────── \n")
    
    # ─── CALCULATE DISTANCES ───────────────────────────────────────────────────
    def calculate_distances(self, method):
        methods = {'hashrf' : hashrf.hashrf,
                   'hashrf_weighted' : hashrf.hashrf_weighted,
                   'None' : None, }
        
        with self.console.status("[bold green]Calculating distances...") as status:
            self.distance_matrix = methods[method](self.file, self.n_trees, self.output_file)
        print(f'[bold red]{method} | Done!')
    
    # ─── EMBED ─────────────────────────────────────────────────────────────────
    def embed_2D(self, method):
        methods = {'pca' : PCA_e.pca,
                   'tsne' : tSNE_e.tsne,
                   'None' : None, }
        
        dimensions = 2
        
        if type(self.distance_matrix) == type(None): self.calculate_distances('hashrf')
            #raise Exception('Distance matrix has to be computed or inserted by the user prior the embedding')
        if method == 'pca':
            self.embedding_pca2D = methods[method](self.distance_matrix, dimensions, self.metadata)
        if method == 'tsne':
            self.embedding_tsne2D = methods[method](self.distance_matrix, dimensions, self.metadata)
    
    def embed_3D(self, method):
        methods = {'pca' : PCA_e.pca,
                   'tsne' : tSNE_e.tsne,
                   'None' : None, }
        
        dimensions = 3
        
        if type(self.distance_matrix) == type(None): self.calculate_distances('hashrf') 
            #raise Exception('Distance matrix has to be computed or inserted by the user prior the embedding')
        if method == 'pca':
            self.embedding_pca3D = methods[method](self.distance_matrix, dimensions, self.metadata)
        if method == 'tsne':
            self.embedding_tsne3D = methods[method](self.distance_matrix, dimensions, self.metadata)
    
    # ─── PLOT EMBEDDING ─────────────────────────────────────────────────────────
    
    def plot_2D(self, method, save=False, name_plot=None, static=False, plot_meta = 'SET-ID', plot_set = None, select = False):
        if type(plot_set) == type(None): plot_set = self.sets
        if method == 'pca':
            if name_plot == None: name_plot='PCA_2D'
            if type(self.embedding_pca2D) == type(None): self.embed_2D('pca')
            fig = graph.plot_embedding(self.embedding_pca2D, self.metadata, 2, save, name_plot, static, plot_meta, plot_set, select)
        
        elif method == 'tsne':
            if name_plot == None: name_plot='TSNE_2D'
            if type(self.embedding_tsne2D) == type(None): self.embed_2D('tsne')
            fig = graph.plot_embedding(self.embedding_tsne2D, self.metadata, 2, save, name_plot, static, plot_meta, plot_set, select)
        
        else: raise ValueError("'method' can only be either 'pca' or 'tsne' ")
        
        return fig
    
    def plot_3D(self, method, save=False, name_plot=None, static=False, plot_meta = 'SET-ID', plot_set = None, select = False):
        if type(plot_set) == type(None): plot_set = self.sets
        if method == 'pca':
            if name_plot == None: name_plot='PCA_3D'
            if type(self.embedding_pca3D) == type(None): self.embed_3D('pca')
            fig = graph.plot_embedding(self.embedding_pca3D, self.metadata, 3, save, name_plot, static, plot_meta, plot_set, select)
        
        elif method == 'tsne':
            if name_plot == None: name_plot='TSNE_3D'
            if type(self.embedding_tsne3D) == type(None): self.embed_3D('tsne')
            fig = graph.plot_embedding(self.embedding_tsne3D, self.metadata, 3, save, name_plot, static, plot_meta, plot_set, select)
        
        else: raise ValueError("'method' can only be either 'pca' or 'tsne' ")
        
        return fig
            
            
# ──────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────── SET_COLLECTION CLASS ─────
class set_collection(tree_set):
    # NB: set_collection is a sub_class of tree_set
    # therefore, most methods are shared between these two classes
    def __init__(self, collection = list(), file = 'Set_collection', 
                 output_file = "./Set_collection_distance_matrix"):
        """ Initialize set_collection
        
        collection: facultative - tree_set or list of tree_sets
        NB: if no collection is given an empty set_collection is generated 
        file: facultative - file with set of phylogenetic trees in newick format
        output_file: facultative - specifies output_file of distance matrix"""
        
        self.id = uuid.uuid4()
        self.file = file + str(self.id)
        self.distance_matrix = None
        self.embedding_pca2D = None
        self.embedding_tsne2D = None
        self.embedding_pca3D = None
        self.embedding_tsne3D = None
        
        if self.file != 'Set_collection' + str(self.id) and output_file == "./Set_collection_distance_matrix.csv": 
            self.output_file = "./{file}_distance_matrix.csv".format(file = os.path.splitext(os.path.basename(self.file))[0])
        else: self.output_file = output_file + str(self.id) + ".csv"
        
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
            if type(metadata) == type(None): metadata = pd.DataFrame()
            metadata['SET-ID'] = np.array([key] * set.n_trees)
            
            self.metadata = pd.concat([self.metadata, metadata])
            
            self.data[key] = {'metadata':metadata, 
                                 'n_trees':set.n_trees}
        
        self.sets = np.unique(self.metadata['SET-ID'])
            
        
    # the result of addition between two collections 
    # is the concatenation of the two collections
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
        
    # concatenate is a more formal method to concatenate collections
    # using this allows for more clarity in the codebase
    #! less error prone !#
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
            
            
            
        
        

         