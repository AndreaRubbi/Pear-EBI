__author__ = "Andrea Rubbi"

# ─── About ────────────────────────────────────────────────────────────────────
"""This file defines the tree_set and set_collection classes.
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
    relative distance in a common embedding. Examples of possible applications are present at: ###LINK###"""
# ──────────────────────────────────────────────────────────────────────────────

__copyright__ = "2023-present Andrea Rubbi and other contributors"
__credits__ = ["Andrea Rubbi", "Lukas Weilguny", "Nick Goldman", "Nicola de Maio"]

__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Andrea Rubbi"
__institute__ = "EMBL-EBI"
__email__ = "andrear@ebi.ac.uk"
__status__ = "Production"

# ──────────────────────────────────────────────────────────────────────────────
import os
import random
import shutil
import subprocess
import sys
import time
import uuid
import warnings

import numpy as np
import pandas as pd
from rich import print
# ? rich is a very nice library that allows to
# ? easily format the output of console
# ? https://github.com/Textualize/rich
from rich.console import Console

# getting the name of the directory
current = os.path.dirname(os.path.realpath(__file__))

# Getting the parent directory name
parent = os.path.dirname(current)

# adding the parent directory to
# the sys.path.
sys.path.append(parent)

# importing other modules
try:
    from .calculate_distances import hashrf, maple_RF, tqdist
    from .embeddings import Isomap_e, LLE_e, PCA_e, tSNE_e
    from .embeddings.graph import graph
    from .interactive_mode import interactive
    from .subsample import subsample
except:
    sys.exit("Error")


# ────────────────────────────────────────────────────────── TREE_SET CLASS ─────
class tree_set:
    """Class for the analysis of a set of phylogenetic trees"""

    # Console from rich -> takes control of console output
    console = Console()

    # ─── INIT ──────────────────────────────────────────────────────────────────
    def __init__(self, file, output_file=None, distance_matrix=None, metadata=None):
        """Initialize tree_set

        file: mandatory - file with set of phylogenetic trees in newick format
        output_file: facultative - specifies output_file of distance matrix
        distance_matrix: facultative - specifies file with already-computed distance matrix
        metadata: facultative - specifies file containing additional information for each tree in set.
        It should contain a column for each feature, a row for each tree (blank row if no info)
        """

        self.file = file
        self.output_file = output_file
        self.distance_matrix = distance_matrix
        self.metadata = metadata
        self.embedding_pca2D = None
        self.embedding_tsne2D = None
        self.embedding_pca3D = None
        self.embedding_tsne3D = None

        if self.output_file == None:
            self.output_file = "./{file}_distance_matrix.csv".format(
                file=os.path.splitext(os.path.basename(self.file))[0]
            )

        self.size = int(f"{os.path.getsize(file)/(1<<30):,.0f}")
        # if self.size > 3: sys.exit(f'File is too large: {self.size} GB')

        try:
            self.n_trees = int(
                subprocess.check_output(["wc", "-l", self.file]).decode().split(" ")[0]
            )
        except:
            with open(file, "r") as f:
                self.n_trees = len(f.readlines())
                f.close()

        if type(self.distance_matrix) != type(None):
            try:
                pd.read_csv(self.distance_matrix)
            except:
                sys.exit(
                    "There's an error with the Distance Matrix file - please check the correct location and name of the .csv file"
                )

        if type(self.metadata) != type(None):
            try:
                self.metadata = pd.read_csv(self.metadata)
            except:
                sys.exit(
                    "There's an error with the Metadata file - please check the correct location and name of the .csv file"
                )

        else:
            self.metadata = pd.DataFrame()
        self.metadata["SET-ID"] = [
            os.path.splitext(os.path.basename(self.file))[0] for i in range(self.n_trees)
        ]
        self.metadata["STEP"] = [i for i in range(self.n_trees)]
        self.sets = np.unique(self.metadata["SET-ID"])

    # ─── STR ───────────────────────────────────────────────────────────────────
    def __str__(self):
        """Returns string representation of tree_set

        Returns:
            __str__: summary of tree_set
        """
        computed = "not computed"
        if type(self.distance_matrix) != type(None):
            computed = "computed"

        return f"─────────────────────────────\n Tree set containing {self.n_trees} trees;\n File: {self.file};\n Distance matrix: {computed}.\n───────────────────────────── \n"

    # ─── CALCULATE DISTANCES ───────────────────────────────────────────────────
    def calculate_distances(self, method):
        """Computes tree_set distance matrix with method of choice

        Args:
            method (str): method/algorithm used to compute distance matrix
        """
        methods = {
            "hashrf": hashrf.hashrf,
            "hashrf_weighted": hashrf.hashrf_weighted,
            "days_RF": maple_RF.calculate_distance_matrix,
            "quartet": tqdist.quartet,
            "triplet": tqdist.triplet,
            "None": None,
        }

        with self.console.status("[bold green]Calculating distances...") as status:
            self.distance_matrix = methods[method](
                self.file, self.n_trees, self.output_file
            )
        print(f"[bold blue]{method} | Done!")

    # ─── EMBED ─────────────────────────────────────────────────────────────────
    def embed(self, method, dimensions, quality=False, report=False):
        """Compute embedding with n-dimensions and method of choice

        Args:
            method (str): method of choice to embed data
            dimensions (_type_): number of dimensions/components
            quality (bool, optional): returns quality report and self.emb_quality. Defaults to False.
        """
        methods = {
            "pca": PCA_e.pca,
            "tsne": tSNE_e.tsne,
            "isomap": Isomap_e.isomap,
            "lle": LLE_e.lle,
            "None": None,
        }

        if type(self.distance_matrix) == type(None):
            self.calculate_distances("hashrf")

        dim = dimensions if dimensions > 2 else 3

        with self.console.status("[bold green]Embedding distances...") as status:
            embedding = methods[method](
                self.distance_matrix,
                dimensions,
                self.metadata,
                quality=quality if not report else True,
                report=report,
            )
        print(f"[bold blue]{method} | Done!")

        if quality:
            if method == "pca":
                embedding, var, corr, self.emb_quality = embedding
                print(
                    f"With {dimensions} components/dimensions, the explained variance is {var:.2f},\n with an estimated correlation {corr[0, 1]:.2f} with the {self.n_trees}-dimensional coordinates"
                )
            else:
                embedding, corr, self.emb_quality = embedding
                print(
                    f"With {dimensions} components/dimensions, the estimated correlation with the {self.n_trees}-dimensional coordinates is {corr[0, 1]:.2f}"
                )

        if method == "pca":
            self.embedding_pca = embedding
            self.embedding_pca3D = embedding[:, :4]
            self.embedding_pca2D = embedding[:, :3]
        elif method == "tsne":
            if dimensions > 3:
                warnings.warn(
                    "t-SNE with more than 3 dimensions can be considerably slow"
                )
            self.embedding_tsne = embedding
            self.embedding_tsne3D = embedding[:, :4]
            self.embedding_tsne2D = embedding[:, :3]
        elif method == "isomap":
            if dimensions > 3:
                warnings.warn(
                    "Isomap with more than 3 dimensions can be considerably slow"
                )
            self.embedding_isomap = embedding
            self.embedding_isomap3D = embedding[:, :4]
            self.embedding_isomap2D = embedding[:, :3]

        elif method == "lle":
            if dimensions > 3:
                warnings.warn("LLE with more than 3 dimensions can be considerably slow")
            self.embedding_lle = embedding
            self.embedding_lle3D = embedding[:, :4]
            self.embedding_lle2D = embedding[:, :3]

    # ─── PLOT EMBEDDING ─────────────────────────────────────────────────────────

    def plot_2D(
        self,
        method,
        save=False,
        name_plot=None,
        static=False,
        plot_meta="SET-ID",
        plot_set=None,
        select=False,
        same_scale=False,
    ):
        """Plot 2D embedding performed with method of choice

        Args:
            method (str): embedding method
            save (bool, optional): save plot HTML. Defaults to False.
            name_plot (str, optional): name of plot's file. Defaults to None.
            static (bool, optional): return less interactive plot. Defaults to False.
            plot_meta (str, optional): meta-variale used to color the points. Defaults to "SET-ID".
            plot_set (list, optional): list of sets to plot from set_collection. Defaults to None.
            select (bool, optional): return set of buttons to show or hide specific traces. Defaults to False.
            same_scale (bool, optional): use same color_scale for all traces when scale is continuous. Defaults to False.

        Raises:
            ValueError: method can only be either pca or tsne for now

        Returns:
            plot: either interactive or not
        """
        if type(plot_set) == type(None):
            plot_set = self.sets
        if method == "pca":
            if name_plot == None:
                name_plot = "PCA_2D"
            if type(self.embedding_pca2D) == type(None):
                self.embed("pca", 2)
            fig = graph.plot_embedding(
                self.embedding_pca2D,
                self.metadata,
                2,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        elif method == "tsne":
            if name_plot == None:
                name_plot = "TSNE_2D"
            if type(self.embedding_tsne2D) == type(None):
                self.embed("tsne", 2)
            fig = graph.plot_embedding(
                self.embedding_tsne2D,
                self.metadata,
                2,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        elif method == "isomap":
            if name_plot == None:
                name_plot = "ISOMAP_2D"
            if type(self.embedding_isomap2D) == type(None):
                self.embed("isomap", 2)
            fig = graph.plot_embedding(
                self.embedding_isomap2D,
                self.metadata,
                2,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        elif method == "lle":
            if name_plot == None:
                name_plot = "LLE_2D"
            if type(self.embedding_lle2D) == type(None):
                self.embed("lle", 2)
            fig = graph.plot_embedding(
                self.embedding_lle2D,
                self.metadata,
                2,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        else:
            raise ValueError("'method' can only be either 'pca' or 'tsne' ")

        return fig

    def plot_3D(
        self,
        method,
        save=False,
        name_plot=None,
        static=False,
        plot_meta="SET-ID",
        plot_set=None,
        select=False,
        same_scale=False,
    ):
        """Plot 3D embedding performed with method of choice

        Args:
            method (str): embedding method
            save (bool, optional): save plot HTML. Defaults to False.
            name_plot (str, optional): name of plot's file. Defaults to None.
            static (bool, optional): return less interactive plot. Defaults to False.
            plot_meta (str, optional): meta-variale used to color the points. Defaults to "SET-ID".
            plot_set (list, optional): list of sets to plot from set_collection. Defaults to None.
            select (bool, optional): return set of buttons to show or hide specific traces. Defaults to False.
            same_scale (bool, optional): use same color_scale for all traces when scale is continuous. Defaults to False.

        Raises:
            ValueError: method can only be either pca or tsne for now

        Returns:
            plot: either interactive or not
        """
        if type(plot_set) == type(None):
            plot_set = self.sets
        if method == "pca":
            if name_plot == None:
                name_plot = "PCA_3D"
            if type(self.embedding_pca3D) == type(None):
                self.embed("pca", 3)
            fig = graph.plot_embedding(
                self.embedding_pca3D,
                self.metadata,
                3,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        elif method == "tsne":
            if name_plot == None:
                name_plot = "TSNE_3D"
            if type(self.embedding_tsne3D) == type(None):
                self.embed("tsne", 3)
            fig = graph.plot_embedding(
                self.embedding_tsne3D,
                self.metadata,
                3,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        elif method == "isomap":
            if name_plot == None:
                name_plot = "ISOMAP_3D"
            if type(self.embedding_isomap3D) == type(None):
                self.embed("isomap", 3)
            fig = graph.plot_embedding(
                self.embedding_isomap3D,
                self.metadata,
                3,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        elif method == "lle":
            if name_plot == None:
                name_plot = "LLE_3D"
            if type(self.embedding_lle3D) == type(None):
                self.embed("lle", 3)
            fig = graph.plot_embedding(
                self.embedding_lle3D,
                self.metadata,
                3,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
            )

        else:
            raise ValueError("'method' can only be either 'pca' or 'tsne' ")

        return fig

    # ─── GET SUBSET ───────────────────────────────────────────────────────

    def get_subset(self, n_required, method="sequence"):
        """Gets subset of phylogenetic trees

        Args:
            n_required (int): number of points to extract
            method (str, optional): method used to extact points ('sequence', 'random', 'syst'). Defaults to "sequence".

        Returns:
            subset plots: 2D and 3D embedding plots of subset
        """
        console = Console()
        with console.status("[bold blue]Extracting subsample...") as status:
            if method == "syst":
                if shutil.which("pypy3") is not None:
                    command = [
                        "pypy3",
                        f"{current}/subsample/subsample.py",
                        self.file,
                        str(self.n_trees),
                        str(n_required),
                    ]
                    res = subprocess.check_output(command, universal_newlines=True).split(
                        "\n"
                    )
                    subsample_trees, idxs = eval(res[3]), eval(res[4])
                else:
                    console.log(
                        "[bold red]Could not find pypy3 on your sytem PATH - using python3..."
                    )
                    subsample_trees, idxs = subsample.subsample(
                        self.file, self.n_trees, n_required, subp=False
                    )

            else:
                with open(self.file, "r") as f:
                    trees = list(enumerate(f.readlines()))
                    f.close()

                if method == "random":
                    selection = random.sample(trees, n_required)
                    subsample_trees, idxs = list(
                        map(lambda elem: elem[1], selection)
                    ), list(map(lambda elem: elem[0], selection))
                elif method == "sequence":
                    step = self.metadata.shape[0] // n_required
                    idxs = [step * (i + 1) - 1 for i in range(n_required)]
                    subsample_trees = [trees[i][1] for i in idxs]

                else:
                    sys.exit(f"Method {method} not available for subsampling")

            file_sub = f"SUBSAMPLE"
            with open(file_sub, "w") as f:
                for i in subsample_trees:
                    f.write(i)
                f.close()
            # print(len(subsample_trees), len(idxs))
            status.update("[bold green]Calculating distances...")
            dM = hashrf.hashrf(file_sub, n_required, file_sub + "_distances.csv")
            components = PCA_e.pca(dM, 3)
            status.update(f"[bold blue] Done!")
            time.sleep(0.2)

        sorted_elements = sorted(enumerate(idxs), key=lambda x: x[1])
        idxs_sorted, order = list(map(lambda x: x[1], sorted_elements)), list(
            map(lambda x: x[0], sorted_elements)
        )
        comp_sorted = np.array([components[line, :] for line in order])
        # dM_sorted = np.array([dM[line,:] for line in order])
        SetID_sub = self.metadata["SET-ID"][idxs_sorted]
        meta_sub = pd.DataFrame({"SET-ID": SetID_sub, "STEP": idxs_sorted})

        fig1, fig2 = graph.plot_embedding(comp_sorted, meta_sub, 2), graph.plot_embedding(
            comp_sorted, meta_sub, 3
        )
        return fig1, fig2


# ──────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────── SET_COLLECTION CLASS ─────
class set_collection(tree_set):
    # NB: set_collection is a sub_class of tree_set
    # therefore, most methods are shared between these two classes
    def __init__(
        self,
        collection=list(),
        file="Set_collection_",
        output_file=None,
        distance_matrix=None,
        metadata=None,
    ):
        """Initialize set_collection

        collection: facultative - tree_set or list of tree_sets
        NB: if no collection is given an empty set_collection is generated
        file: facultative - file with set of phylogenetic trees in newick format
        output_file: facultative - specifies output_file of distance matrix"""

        self.id = uuid.uuid4()
        self.file = file + str(self.id)
        self.distance_matrix = (
            pd.read_csv(distance_matrix, header=None).values
            if distance_matrix
            else distance_matrix
        )
        self.embedding_pca2D = None
        self.embedding_tsne2D = None
        self.embedding_pca3D = None
        self.embedding_tsne3D = None

        if self.file != "Set_collection_" + str(self.id) and output_file is None:
            self.output_file = "{file}_distance_matrix.csv".format(
                file=os.path.splitext(os.path.basename(self.file))[0]
            )
        elif output_file is None:
            self.output_file = "Set_collection_distance_matrix_" + str(self.id) + ".csv"
        else:
            if output_file[-4:] == ".csv":
                self.output_file = output_file[:-4] + "_" + str(self.id) + ".csv"
            else:
                self.output_file = output_file + "_" + str(self.id) + ".csv"

        if isinstance(collection, tree_set):
            self.collection = [collection]
            with open(self.file, "w") as trees:
                with open(collection.file, "r") as file:
                    trees.write(file.read())
                    file.close()
                trees.close()

        elif len(collection) > 0:
            for element in collection:
                assert isinstance(
                    element, tree_set
                ), "Every element in a set_collection must be a tree_set"
            self.collection = collection
            with open(self.file, "w") as trees:
                for set in collection:
                    with open(set.file, "r") as file:
                        trees.write(file.read())
                        file.close()
                trees.close()

        else:
            self.collection = collection

        with open(self.file, "r") as f:
            self.n_trees = len(f.readlines())
            f.close()

        self.data = dict()

        self.metadata = pd.DataFrame()

        for set in self.collection:
            key = os.path.splitext(os.path.basename(set.file))[0]

            metadata = set.metadata
            if type(metadata) == type(None):
                metadata = pd.DataFrame()
            metadata["SET-ID"] = np.array([key] * set.n_trees)

            self.metadata = pd.concat([self.metadata, metadata])

            self.data[key] = {"metadata": metadata, "n_trees": set.n_trees}

        self.metadata.reset_index(drop=True, inplace=True)

        self.sets = np.unique(self.metadata["SET-ID"])

    # the result of addition between two collections
    # is the concatenation of the two collections
    def __add__(self, other):
        """Concatenates two collectionsor collection and tree_set

        Args:
            other (tree_set ot set_colletion): tree_set ot set_colletion

        Returns:
            set_collection: concatenated set_collection
        """
        try:
            assert isinstance(other, set_collection)
            return set_collection(self.collection + other.collection)
        except:
            try:
                assert isinstance(other, tree_set)
                return set_collection(self.collection + [other])
            except:
                for element in other:
                    assert isinstance(
                        element, tree_set
                    ), "You can concatenate a set_collection \
                        only with another set_collection, a tree_set,\
                            or a list of tree_set"
                return set_collection(self.collection + other)

    def __str__(self):
        computed = "not computed"
        if type(self.distance_matrix) != type(None):
            computed = "computed"

        summary = f"─────────────────────────────\
            \n Tree set collection containing {self.n_trees} trees;\
            \n File: {self.file};\n Distance matrix: {computed}.\
                \n───────────────────────────── \n"
        for key, value in self.data.items():
            summary += f"{key}; Containing {value['n_trees']} trees. \n"

        return summary

    # concatenate is a more formal method to concatenate collections
    # using this allows for more clarity in the codebase
    def concatenate(self, other):
        """Concatenates two collectionsor collection and tree_set

        Args:
            other (tree_set ot set_colletion): tree_set ot set_colletion

        Returns:
            set_collection: concatenated set_collection
        """
        try:
            assert isinstance(other, tree_set)
            return set_collection(self.collection + [other])

        except:
            try:
                assert isinstance(other, set_collection)
                return set_collection(self.collection + other.collection)
            except:
                for element in other:
                    assert isinstance(
                        element, tree_set
                    ), "You can concatenate a set_collection only with another \
                        set_collection, a tree_set, or a list of tree_set"
                return set_collection(self.collection + other)
