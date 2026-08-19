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
import tempfile
import time
import uuid
import warnings
import re

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

# silencing some warnings
from scipy.sparse import SparseEfficiencyWarning

from .calculate_distances import _exec, hashrf, maple_RF, tqdist
from .embeddings import Isomap_e, LLE_e, PCoA_e, tSNE_e
from .embeddings.graph import graph
# from .interactive_mode import interactive
from .subsample import subsample

# importing other modules



warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)


# ─── NEWICK READING / NORMALISATION ───────────────────────────────────────────
# Every path that needs to know how many trees a file holds, or to hand a file to
# one of the native tools, goes through _read_trees()/_write_trees(). Before this,
# tree counting was spread over six sites using three different idioms (semicolon
# counts, `wc -l`, and len(readlines())), which is how the reported bug arose: a
# file whose last tree has no trailing newline counts as 0 lines under `wc -l`,
# and 0 was then passed to hashrf.
#
# Correct newline termination is not cosmetic. hashrf's parser is whitespace
# agnostic and tolerates several trees on one line *given the right count*, but
# tqDist's is strictly line-based (NewickParser.cpp:86-100): it getline()s, and a
# final line with no "\n" sets eofbit so the loop breaks before using it, silently
# dropping the last tree. maple_RF splits on lines too. So the tools need one tree
# per line, newline-terminated.


def _split_newick(text):
    """Split a Newick stream into tree strings on top-level ';' terminators.

    A plain text.count(";") over-counts, because a ';' can also appear inside a
    [...] comment or a '...' quoted label, and NEXUS files use ';' to end their own
    statements. This walks the text tracking comment depth and quoting, and drops
    newlines and tabs outside quoted labels so each tree comes back on one line.
    """
    trees = []
    buf = []
    comment_depth = 0
    in_quote = False
    i = 0
    n = len(text)

    while i < n:
        ch = text[i]

        if in_quote:
            buf.append(ch)
            if ch == "'":
                # '' is an escaped quote within a quoted label
                if i + 1 < n and text[i + 1] == "'":
                    buf.append(text[i + 1])
                    i += 2
                    continue
                in_quote = False
            i += 1
            continue

        if ch in "\r\n\t":
            # Structural whitespace only: real whitespace inside a quoted label is
            # preserved by the branch above.
            i += 1
            continue

        if comment_depth:
            buf.append(ch)
            if ch == "[":
                comment_depth += 1
            elif ch == "]":
                comment_depth -= 1
            i += 1
            continue

        if ch == "'":
            in_quote = True
            buf.append(ch)
        elif ch == "[":
            comment_depth = 1
            buf.append(ch)
        elif ch == ";":
            buf.append(";")
            trees.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
        i += 1

    tail = "".join(buf).strip()
    if tail:
        # Trailing content with no ';' at all: a tree missing its terminator.
        trees.append(tail + ";")

    # Drop statements that are only a terminator or whitespace.
    return [t for t in trees if t.strip(" \t\r\n;")]


def _extract_nexus_trees(text):
    """Pull the Newick strings out of the TREES block of a NEXUS file.

    BEAST and friends emit NEXUS with a header, a `translate` table and one
    `tree STATE_n = [&R] (...);` statement per tree. Reusing _split_newick means
    the ';' ending `begin trees;`, each translate entry and `end;` are separated
    out as their own statements and simply discarded here.

    The translate table is deliberately not applied: RF, quartet and triplet
    distances depend only on how labels partition, so leaving the numeric labels
    in place gives identical distances and avoids a whole class of mis-mapping.
    """
    trees = []
    for statement in _split_newick(text):
        head = statement.lstrip().lower()
        if not (head.startswith("tree ") or head.startswith("utree ")):
            continue
        _, sep, rhs = statement.partition("=")
        if not sep:
            continue
        rhs = rhs.strip()
        # Skip rooting/annotation comments such as [&R] that precede the topology.
        while rhs.startswith("["):
            close = rhs.find("]")
            if close == -1:
                break
            rhs = rhs[close + 1 :].lstrip()
        if rhs:
            trees.append(rhs if rhs.endswith(";") else rhs + ";")
    return trees


def _read_trees(path):
    """Read a Newick or NEXUS tree file and return one normalised string per tree.

    Each returned string is a single line ending in ';'. This is the single source
    of truth for how many trees a file contains.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        content = fh.read()

    if content.lstrip()[:6].upper() == "#NEXUS":
        trees = _extract_nexus_trees(content)
    else:
        trees = _split_newick(content)

    normalised = []
    for tree in trees:
        tree = tree.strip()
        if not tree.endswith(";"):
            tree += ";"
        normalised.append(tree)
    return normalised


def _write_trees(handle, trees):
    """Write tree strings one per line, each newline-terminated.

    Returns the number of trees written.
    """
    count = 0
    for tree in trees:
        handle.write(tree.rstrip("\n") + "\n")
        count += 1
    return count


def _clean_metadata_df(df):
    """Normalize metadata DataFrame:
    - Strip column names
    - Drop rows that are entirely empty/NaN
    - If the first data row equals the column names (duplicate header), drop it
    - Reset index
    """
    try:
        # Ensure columns are stripped strings
        df.columns = [str(c).strip() for c in df.columns]

        # Replace empty-string-only cells with NaN for detection
        df = df.replace("", np.nan)

        # Drop rows that are completely NaN
        df = df.dropna(axis=0, how="all")

        # Reset index to standard 0..n-1
        df = df.reset_index(drop=True)

        # Detect duplicate header row: compare first row values (as strings) to column names
        if df.shape[0] > 0:
            first_row = [str(x).strip() for x in df.iloc[0].astype(object).tolist()]
            cols = [str(c).strip() for c in df.columns.tolist()]
            # If identical, drop the first row
            if first_row == cols:
                df = df.iloc[1:].reset_index(drop=True)

        return df
    except Exception:
        return df


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
        self.embedding_pcoa2D = None
        self.embedding_tsne2D = None
        self.embedding_pcoa3D = None
        self.embedding_tsne3D = None

        if self.output_file == None:
            self.output_file = "./{file}_distance_matrix.csv".format(
                file=os.path.splitext(os.path.basename(self.file))[0]
            )

        self.size = os.path.getsize(file) // (1 << 30)
        # if self.size > 3: sys.exit(f'File is too large: {self.size} GB')

        # _read_trees() is the single source of truth for the tree count, replacing
        # the previous four-way mixture of semicolon counting, `wc -l` and
        # len(readlines()). `wc -l` is what produced the reported bug: it counts
        # "\n" bytes, so a lone tree with no trailing newline counted as 0 trees and
        # 0 was handed to hashrf.
        trees = _read_trees(file)
        self.n_trees = len(trees)
        if self.n_trees == 0:
            sys.exit(
                f"No phylogenetic trees found in {file}.\n"
                "Expected a Newick file (one or more ';'-terminated trees) or a "
                "NEXUS file with a trees block."
            )

        # The native tools need one tree per line, newline-terminated, so a
        # normalised copy is written to a temp directory. self.file deliberately
        # keeps pointing at the user's path: it used to be reassigned to the temp
        # file, which then leaked into every derived name -- __str__, the embedding
        # CSV, all eight plot filenames and the _SUBSAMPLE file (which ended up
        # written into /tmp). Use tool_input() to get the path for the tools.
        #
        # TemporaryDirectory cleans up through weakref.finalize, which also runs at
        # interpreter exit; the previous delete=False temp file relied on __del__
        # and leaked whenever that did not run.
        self._tmpdir = tempfile.TemporaryDirectory(prefix="pear_ebi_")
        self._normalized_file = os.path.join(
            self._tmpdir.name, os.path.basename(file) + ".normalised"
        )
        with open(self._normalized_file, "w", encoding="utf-8") as fh:
            _write_trees(fh, trees)

        if type(self.distance_matrix) != type(None):
            try:
                self.distance_matrix = pd.read_csv(
                    self.distance_matrix, header=None, index_col=None
                ).values
                #    header=0,
                #    index_col=0,
                #    dtype=np.float32,
                # self.distance_matrix.columns = list(range(self.distance_matrix.shape[1]))
            except:
                sys.exit(
                    "There's an error with the Distance Matrix file - please check the correct location and name of the .csv file"
                )

        if type(self.metadata) != type(None):
            try:
                self.metadata = pd.read_csv(self.metadata)
                # Clean metadata: remove empty trailing rows and duplicate header rows
                try:
                    self.metadata = _clean_metadata_df(self.metadata)
                except Exception:
                    pass
            except Exception:
                sys.exit(
                    "There's an error with the Metadata file - please check the correct location and name of the .csv file"
                )

        else:
            self.metadata = pd.DataFrame()
            # self.file is the user's path (the normalised copy lives separately,
            # see tool_input), so SET-ID matches the name used in config files.
            set_id_base = os.path.splitext(os.path.basename(self.file))[0]
            self.metadata["SET-ID"] = [set_id_base for i in range(self.n_trees)]
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

    def tool_input(self):
        """Path to hand to the native tools.

        The newline-normalised copy when there is one, otherwise self.file. Always
        use this rather than self.file when invoking hashrf, tqDist or maple_RF:
        tqDist and maple_RF read one tree per line and drop a final tree that has
        no trailing newline.

        There is no __del__ any more. It existed only to unlink the old
        delete=False temp file, and __del__ is not guaranteed to run -- notably at
        interpreter exit, and in notebooks where the object stays referenced.
        TemporaryDirectory handles cleanup through weakref.finalize instead.
        """
        return getattr(self, "_normalized_file", None) or self.file

    # ─── CALCULATE DISTANCES ───────────────────────────────────────────────────
    def calculate_distances(self, method):
        """Computes tree_set distance matrix with method of choice

        Args:
            method (str): method/algorithm used to compute distance matrix
        """
        methods = {
            "hashrf_RF": hashrf.hashrf,
            "hashrf_wRF": hashrf.hashrf_weighted,
            "smart_RF": maple_RF.calculate_distance_matrix,
            "tqdist_quartet": tqdist.quartet,
            "tqdist_triplet": tqdist.triplet,
            "None": None,
        }

        with self.console.status("[bold green]Calculating distances...") as status:
            self.distance_matrix = methods[method](
                self.tool_input(), self.n_trees, self.output_file
            )
        print(f"[bold blue]{method} | Done!")

    # ─── EMBED ─────────────────────────────────────────────────────────────────
    def embed(self, method, dimensions, quality=False, report=False, output=None):
        """Compute embedding with n-dimensions and method of choice

        Args:
            method (str): method of choice to embed data
            dimensions (_type_): number of dimensions/components
            quality (bool, optional): returns quality report and self.emb_quality. Defaults to False.
        """
        methods = {
            "pcoa": PCoA_e.pcoa,
            "tsne": tSNE_e.tsne,
            "isomap": Isomap_e.isomap,
            "lle": LLE_e.lle,
            "None": None,
        }

        if type(self.distance_matrix) == type(None):
            self.calculate_distances("hashrf_RF")

        # Sanity check: embedding routines expect self.distance_matrix rows to
        # match the number of metadata rows. If they don't, try to reconcile by
        # trimming metadata to the distance matrix size and warn the user.
        n_data = (
            self.distance_matrix.shape[0]
            if hasattr(self.distance_matrix, "shape")
            else len(self.distance_matrix)
        )

        if self.metadata.shape[0] != n_data:
            # This used to trim the metadata down to the matrix and carry on with a
            # yellow warning, which meant a wrong tree count reached the plots
            # looking plausible -- trees silently dropped off the end, and every
            # point mislabelled if the mismatch was not at the tail. There is no
            # safe way to guess which rows to discard, so stop.
            per_set = ""
            if hasattr(self, "data"):
                per_set = "\n  Per-set tree counts:\n" + "\n".join(
                    f"    - {name}: {info['n_trees']}"
                    for name, info in self.data.items()
                )
            sys.exit(
                f"Shape mismatch: the distance matrix has {n_data} rows but the "
                f"metadata has {self.metadata.shape[0]}.\n"
                f"  These must agree -- one row per tree.\n"
                f"  Most often this means a precomputed distance matrix (-d/--dM) was "
                f"reused with a different set of trees, or a metadata CSV has the "
                f"wrong number of rows.{per_set}"
            )

        if dimensions < 2:
            sys.exit("Dimensions of embedding must be greater or equal to 2")

        if output == None:
            output = f"./{os.path.splitext(os.path.basename(self.file))[0]}_{str.upper(method)}_embedding.csv"
        with self.console.status("[bold green]Embedding distances...") as status:
            embedding = methods[method](
                self.distance_matrix,
                dimensions,
                self.metadata,
                quality=quality if not report else True,
                report=report,
                output=output,
            )
        print(f"[bold blue]{method} | Done!")

        if quality:
            if method == "pcoa":
                embedding, var, corr, self.emb_quality = embedding
                print(
                    f"With {dimensions} components/dimensions, the explained variance is {var:.2f},\n with an estimated correlation {corr[0, 1]:.2f} with the {self.n_trees}-dimensional coordinates"
                )
            else:
                embedding, corr, self.emb_quality = embedding
                print(
                    f"With {dimensions} components/dimensions, the estimated correlation with the {self.n_trees}-dimensional coordinates is {corr[0, 1]:.2f}"
                )

        if method == "pcoa":
            self.embedding_pcoa = embedding
            if dimensions == 2:
                self.embedding_pcoa2D = embedding
            if dimensions == 3:
                self.embedding_pcoa2D = embedding[:, :3]
                self.embedding_pcoa3D = embedding
            if dimensions > 3:
                self.embedding_pcoa3D = embedding[:, :4]
                self.embedding_pcoa2D = embedding[:, :3]

        elif method == "tsne":
            if dimensions > 3:
                warnings.warn(
                    "t-SNE with more than 3 dimensions can be considerably slow"
                )
            self.embedding_tsne = embedding
            if dimensions == 2:
                self.embedding_tsne2D = embedding
            if dimensions == 3:
                self.embedding_tsne2D = embedding[:, :3]
                self.embedding_tsne3D = embedding
            if dimensions > 3:
                self.embedding_tsne3D = embedding[:, :4]
                self.embedding_tsne2D = embedding[:, :3]

        elif method == "isomap":
            if dimensions > 3:
                warnings.warn(
                    "Isomap with more than 3 dimensions can be considerably slow"
                )
            self.embedding_isomap = embedding
            if dimensions == 2:
                self.embedding_isomap2D = embedding
            if dimensions == 3:
                self.embedding_isomap2D = embedding[:, :3]
                self.embedding_isomap3D = embedding
            if dimensions > 3:
                self.embedding_isomap3D = embedding[:, :4]
                self.embedding_isomap2D = embedding[:, :3]

        elif method == "lle":
            if dimensions > 3:
                warnings.warn("LLE with more than 3 dimensions can be considerably slow")
            self.embedding_lle = embedding
            if dimensions == 2:
                self.embedding_lle2D = embedding
            if dimensions == 3:
                self.embedding_lle2D = embedding[:, :3]
                self.embedding_lle3D = embedding
            if dimensions > 3:
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
            ValueError: method can only be either pcoa or tsne for now

        Returns:
            plot: either interactive or not
        """

        # you can surely write something better here @andrear
        if type(plot_set) == type(None):
            plot_set = self.sets
        if method == "pcoa":
            if name_plot == None:
                name_plot = (
                    f"./{os.path.splitext(os.path.basename(self.file))[0]}_PCOA_2D"
                )
            if type(self.embedding_pcoa2D) == type(None):
                self.embed("pcoa", 2)
            fig = graph.plot_embedding(
                self.embedding_pcoa2D,
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
                name_plot = (
                    f"./{os.path.splitext(os.path.basename(self.file))[0]}_TSNE_2D"
                )
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
                name_plot = (
                    f"./{os.path.splitext(os.path.basename(self.file))[0]}_ISOMAP_2D"
                )
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
                name_plot = f"./{os.path.splitext(os.path.basename(self.file))[0]}_LLE_2D"
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
            raise ValueError("'method' can only be either 'pcoa' or 'tsne' ")

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
        z_axis=None,
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
            ValueError: method can only be either pcoa or tsne for now

        Returns:
            plot: either interactive or not
        """
        if type(plot_set) == type(None):
            plot_set = self.sets
        if method == "pcoa":
            if name_plot == None:
                name_plot = (
                    f"./{os.path.splitext(os.path.basename(self.file))[0]}_PCOA_3D"
                )
            if type(self.embedding_pcoa3D) == type(None):
                self.embed("pcoa", 3)
            fig = graph.plot_embedding(
                self.embedding_pcoa3D,
                self.metadata,
                3,
                save,
                name_plot,
                static,
                plot_meta,
                plot_set,
                select,
                same_scale,
                z_axis,
            )

        elif method == "tsne":
            if name_plot == None:
                name_plot = (
                    f"./{os.path.splitext(os.path.basename(self.file))[0]}_TSNE_3D"
                )
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
                z_axis,
            )

        elif method == "isomap":
            if name_plot == None:
                name_plot = (
                    f"./{os.path.splitext(os.path.basename(self.file))[0]}_ISOMAP_3D"
                )
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
                z_axis,
            )

        elif method == "lle":
            if name_plot == None:
                name_plot = f"./{os.path.splitext(os.path.basename(self.file))[0]}_LLE_3D"
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
                z_axis,
            )

        else:
            raise ValueError("'method' can only be either 'pcoa' or 'tsne' ")

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
        # tool_input() rather than file: the normalised copies are one tree per
        # line, so the line-based reading further down cannot lose a tree whose
        # source file had no trailing newline.
        if isinstance(self, set_collection):
            files = [TS.tool_input() for TS in self.collection]
        else:
            files = [self.tool_input()]
        console = Console()
        with console.status("[bold blue]Extracting subsample...") as status:
            if method == "syst":
                if shutil.which("pypy3") is not None:
                    command = [
                        "pypy3",
                        f"{current}/subsample/subsample.py",
                        str(files),
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
                        str(files), self.n_trees, n_required, subp=False
                    )

            else:
                trees = list()
                last_max = 0
                for file in files:
                    trees_file = _read_trees(file)
                    n_t = len(trees_file)
                    trees.extend(list(enumerate(trees_file, start=last_max)))
                    last_max += n_t

                if method == "random":
                    selection = random.sample(trees, n_required)
                    subsample_trees, idxs = list(
                        map(lambda elem: elem[1], selection)
                    ), list(map(lambda elem: elem[0], selection))
                elif method == "sequence":
                    if n_required < 1 or n_required > len(trees):
                        sys.exit(
                            f"Cannot take a subsample of {n_required} from "
                            f"{len(trees)} trees."
                        )
                    # Guarded because step == 0 (n_required > n_trees) previously
                    # produced idxs = [-1] * n_required, i.e. the last tree repeated,
                    # with no warning; n_required == 0 raised ZeroDivisionError.
                    step = len(trees) // n_required
                    idxs = [step * (i + 1) - 1 for i in range(n_required)]
                    subsample_trees = [trees[i][1] for i in idxs]

                else:
                    sys.exit(f"Method {method} not available for subsampling")

            # Written into a temp dir: this used to be f"{self.file}_SUBSAMPLE", and
            # since self.file had been reassigned to the normalisation temp file, the
            # subsample and its distance matrix were written into /tmp under a random
            # name and never cleaned up.
            self._subsample_tmpdir = tempfile.TemporaryDirectory(prefix="pear_ebi_sub_")
            file_sub = os.path.join(
                self._subsample_tmpdir.name,
                os.path.basename(self.file) + "_SUBSAMPLE",
            )
            with open(file_sub, "w", encoding="utf-8") as f:
                _write_trees(f, subsample_trees)
            # print(len(subsample_trees), len(idxs))
            status.update("[bold green]Calculating distances...")
            dM = hashrf.hashrf(file_sub, n_required, file_sub + "_distances.csv")
            components = PCoA_e.pcoa(dM, 3)
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
        self.file = file + str(self.id) if file == "Set_collection_" else file
        self.distance_matrix = (
            pd.read_csv(distance_matrix, header=None, index_col=None).values  #
            if distance_matrix
            else distance_matrix
        )
        self.embedding_pcoa2D = None
        self.embedding_tsne2D = None
        self.embedding_pcoa3D = None
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
            # with open(self.file, "w") as trees:
            #    with open(collection.file, "r") as file:
            #        trees.write(file.read())
            #        file.close()
            #    trees.close()

        elif len(collection) > 0:
            remove = list()
            for i, element in enumerate(collection):
                if not isinstance(element, tree_set):
                    if isinstance(element, str):
                        try:
                            file = os.path.splitext(os.path.basename(element))[0]
                            exec(f"{file} = tree_set('{element}')")
                            remove.append(i)
                        except FileNotFoundError:
                            sys.exit(f"File {element} not found")
                        except TypeError:
                            sys.exit(
                                f"Set collection can be initialized only with set_collection, tree_set, or file path elements"
                            )
                        exec(f"collection.append({file})")

                    else:
                        sys.exit(
                            f"Set collection can be initialized only with set_collection, tree_set, or file path elements"
                        )
            for i in remove[::-1]:
                collection.pop(i)

            self.collection = collection

        else:
            self.collection = collection

        self.data = dict()

        metadata_input = None
        if type(metadata) != type(None):
            try:
                # Read without forcing an index so we can clean trailing empty rows
                metadata_input = pd.read_csv(metadata)
                try:
                    metadata_input = _clean_metadata_df(metadata_input)
                except Exception:
                    pass
            except Exception:
                sys.exit(
                    "There's an error with the Metadata file - please check the correct location and name of the .csv file"
                )

        self.metadata = pd.DataFrame()

        self.n_trees = 0
        for set in self.collection:
            # Prefer the SET-ID already stored in the member's metadata
            # (tree_set will set this to the original filename base). If
            # it's not present, fall back to the member file basename.
            metadata = set.metadata
            if type(metadata) == type(None):
                metadata = pd.DataFrame()

            if "SET-ID" in metadata.columns and metadata.shape[0] > 0:
                key = str(metadata["SET-ID"].iloc[0])
            else:
                key = os.path.splitext(os.path.basename(set.file))[0]

            # Ensure the metadata exported for the collection uses the
            # canonical SET-ID (so config keys match)
            metadata["SET-ID"] = np.array([key] * set.n_trees)

            self.metadata = pd.concat([self.metadata, metadata])

            self.data[key] = {"metadata": metadata, "n_trees": set.n_trees}
            self.n_trees += set.n_trees

        self.metadata.reset_index(drop=True, inplace=True)

        if type(metadata) != type(None):
            try:
                self.metadata = pd.concat([self.metadata, metadata_input], axis=1)
            except Exception:
                sys.exit(
                    "There's an error with the Metadata file - please check that the .csv file is compatible with the set_collection structure"
                )

        # Ensure SET-ID column is string (avoid mixed types from NaN)
        try:
            if "SET-ID" in self.metadata.columns:
                self.metadata["SET-ID"] = self.metadata["SET-ID"].fillna("").astype(str)
        except Exception:
            pass

        # Compute unique sets, excluding any empty SET-ID entries
        try:
            self.sets = np.unique([s for s in self.metadata["SET-ID"].values if s != ""])
        except Exception:
            self.sets = np.unique(self.metadata["SET-ID"].astype(str).values)

    # ─── CALCULATE DISTANCES ───────────────────────────────────────────────────
    def calculate_distances(self, method):
        """Computes tree_set distance matrix with method of choice

        Args:
            method (str): method/algorithm used to compute distance matrix
        """
        methods = {
            "hashrf_RF": hashrf.hashrf,
            "hashrf_wRF": hashrf.hashrf_weighted,
            "smart_RF": maple_RF.calculate_distance_matrix,
            "tqdist_quartet": tqdist.quartet,
            "tqdist_triplet": tqdist.triplet,
            "None": None,
        }

        combined_file = self.file
        if method in (
            "hashrf_RF",
            "hashrf_wRF",
            "smart_RF",
            "tqdist_quartet",
            "tqdist_triplet",
        ):
            if not self.collection:
                sys.exit(
                    "Cannot compute distances: this collection contains no tree sets."
                )

            # This is the site of the reported bug. The old code did
            #     trees.write(file.read())
            # per member, with no separator, so a member file whose last tree had no
            # trailing newline had its last tree glued onto the next member's first
            # tree. Writing through _write_trees() terminates every tree, whatever
            # shape the member files were in.
            #
            # It also wrote into self.file with mode "w". For a collection built with
            # an explicit file= argument that truncates the user's own file, so the
            # combined input goes to a temp file instead and self.file is left alone.
            self._tmpdir = tempfile.TemporaryDirectory(prefix="pear_ebi_collection_")
            combined_file = os.path.join(
                self._tmpdir.name, os.path.basename(self.file) + ".combined"
            )

            written = 0
            per_member = []
            with open(combined_file, "w", encoding="utf-8") as out:
                for member in self.collection:
                    trees = _read_trees(member.tool_input())
                    written += _write_trees(out, trees)
                    per_member.append((os.path.basename(member.file), len(trees)))

            print("[cyan]Per-member tree counts:")
            for name, count in per_member:
                print(f"  - {name}: {count} trees")

            # The recount used to be printed and thrown away, so a mismatch between
            # the summed member counts and the combined file went to the tools
            # unnoticed. Disagreement here means the count and the data have diverged,
            # which is exactly what produces a wrong-shaped distance matrix.
            if written != self.n_trees:
                sys.exit(
                    f"Tree count mismatch: the combined input holds {written} trees but "
                    f"this collection reports {self.n_trees}.\n"
                    "Refusing to compute distances, because the resulting matrix would "
                    "not match the metadata."
                )
            print(f"[cyan]Combined input: {written} trees, one per line")

        with self.console.status("[bold green]Calculating distances...") as status:
            self.distance_matrix = methods[method](
                combined_file, self.n_trees, self.output_file
            )

        # The combined input lives in a TemporaryDirectory now, so it is cleaned up
        # even if this call raises. The old code shelled out to `rm` with an unquoted
        # path, sent the result to DEVNULL, and omitted smart_RF from the list, which
        # left that method's Set_collection_<uuid> file in the working directory.
        if combined_file != self.file:
            _exec.remove_file(combined_file)

        print(f"[bold blue]{method} | Done!")

    # the result of addition between two collections
    # is the concatenation of the two collections
    def __add__(self, other):
        """Concatenates two collectionsor collection and tree_set

        Args:
            other (tree_set ot set_colletion): tree_set ot set_colletion

        Returns:
            set_collection: concatenated set_collection
        """
        if isinstance(other, set_collection):
            return set_collection(self.collection + other.collection)
        elif isinstance(other, tree_set):
            return set_collection(self.collection + [other])
        else:
            remove = list()
            for i, element in enumerate(other):
                if not isinstance(element, tree_set):
                    if isinstance(element, str):
                        try:
                            file = os.path.splitext(os.path.basename(element))[0]
                            exec(f"{file} = tree_set('{element}')")
                            remove.append(i)
                        except FileNotFoundError:
                            sys.exit(f"File {element} not found")
                        except TypeError:
                            sys.exit(
                                "You can concatenate a set_collection \
                        only with another set_collection, a tree_set,\
                            or a list of tree_set"
                            )

                        exec(f"other.append({file})")

                    else:
                        sys.exit(
                            "You can concatenate a set_collection \
                        only with another set_collection, a tree_set,\
                            or a list of tree_set"
                        )
            for i in remove[::-1]:
                other.pop(i)

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
        if isinstance(other, set_collection):
            return set_collection(self.collection + other.collection)
        elif isinstance(other, tree_set):
            return set_collection(self.collection + [other])
        else:
            remove = list()
            for i, element in enumerate(other):
                if not isinstance(element, tree_set):
                    if isinstance(element, str):
                        try:
                            file = os.path.splitext(os.path.basename(element))[0]
                            exec(f"{file} = tree_set('{element}')")
                            remove.append(i)
                        except FileNotFoundError:
                            sys.exit(f"File {element} not found")
                        except TypeError:
                            sys.exit(
                                "You can concatenate a set_collection \
                        only with another set_collection, a tree_set,\
                            or a list of tree_set"
                            )

                        exec(f"other.append({file})")

                    else:
                        sys.exit(
                            "You can concatenate a set_collection \
                        only with another set_collection, a tree_set,\
                            or a list of tree_set"
                        )
            for i in remove[::-1]:
                other.pop(i)

            return set_collection(self.collection + other)
