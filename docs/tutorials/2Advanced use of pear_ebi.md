# 2. Advanced uses of Pear
Here we introduce some more elaborated ways we can use pear to perform in-depth analyses and produce effective representations of the embedded distance matrices.

## pear.toml
Pear will automatically look for a `pear.toml` file in the working directory, alternatively, a `.toml` file can be specified using the `--config` flag. A `.toml` file is just a convenient way of specifying many parameters. Doing so in an auxiliary file remove excessive clutter in the use of pear and promotes a more standardized way of performing a series of analyses.
<br>We will guide you through the use of this additional tool.


```python
!cat template_pear.toml
```











































Here you can see the `pear.toml` example stored in the same directory of this notebook. It is a perfect template for your future analyses!
<br>Let's go thorugh all its parts:<ul>
    <li> `[trees]` contains single file specifications. Each line associated with this key should direct pear to a file containing trees in Newick format. The nomenclature is "file$n$=filename", where $n$ is just the index of the file, whereas filename is the path to the file itself.
    <li> `[dir]` contains directory and pattern specifications. Each directory should contain only tree-containing files, and should be indicated with a . Alternatively, a `pattern` can be indicated to narrow the research of the files.
    <li> `[collection]` stores details related to the `tree_set` or `set_collection`:<ul>
        <li> `output_file` specifies an alternative name and path for the distance matrix file;
        <li> `distance_matrix` indicates the path of a precomputed distance matrix;
        <li> `metadata` indicates a `.csv` file containing metadata compatible with the collection. That means that the number of rows in the file should be equal to the number of trees in the collection. The information stored in metadata can be of any type (discrete or continuous) and can subsequently used in the representation of your data in the 3D embedding instead of the 3$^{rd}$ dimension, or to color the points (trees). </ul>
    <li> `[highlight]` allows for specifying specific trees in the set/collection which are going to highlighted in the final plots. The way one specifies this is by giving a list of indexes indicating which trees to be highlighted for a given set (either if that is part of a collection or not). You specify a list for a set by writing "file$n$"if the file has been indexed as such in the `[files]` argument. Otherwise by using the name of the file (without extension: filename.trees is just filename) if the file has been specified through `[dir]` selection.
    <li> `[distance]` specifies the `method` used to compute the distance matrix. It can be chosen among `hashrf_RF`, `hashrf_wRF`, `smart_RF`, `tqdist_quartet`, `tqdist_triplet`.
    <li> `[embedding]` specifies the `method` used to compute the embedding of the distance matrix, the `dimensions` of the embedding, and whether to display the `quality` or not. Methods are `pcoa`, `tsne`, `isomap`, `lle`.
    <li> `[plot]` defines some aspects of the plots produced by pear:<ul>
        <li>`name_plot` specifies the name of the plot produced;
        <li>`plot_meta` indicates which feature to use to color the points in the graph, default value is `SET-ID` which simply colors by `tree_set`. A `STEP` meta-variable is present and indicates the index of a tree in a `tree_set`, it can be used to color trees when the ordering is important. Other meta-variables can be specified through the `[metadata]` argument.
        <li> `z_axis`, similarly to `plot_meta`, indicates an alternative meta-variable to be used in the plots. The selected meta-variable replaces the 3$^{rd}$ dimension in the 3D graphs.
        <li>`select` indicates whether the graph should have a set of interactive buttons to display/hide specific `tree_set`s or not.
        <li>`same_scale` indicates whether the same colorscale should be applied to every `tree_set` or not.
        <li>`show` specifies whether the plot should be shown or not.
</ul>
<b>Note that</b> all these arguments are optional, and many of them can be specified otherwise using the normal functionalities of pear. In fact, should any of these arguments be specified using the `.toml` structure and the flags in pear, the arguments will be overwrited by the ones indicated on the command line. On an additional note related to this, the flag `--meta` allows to specify on the command line a metadata file, replicating the behaviour of the `metadata` argument in `[collection]`.

****

## Examples

### Example 1
We use pear to analyze 3 runs of a MCMC algorithm, called Beast, used to estimate a phylogenetic tree structure.<br>
We upload 2 files using the `[trees]` argument and 1 by specifying the `[dir]` and a `pattern`.<br>
We compute the Robinson Foulds distances using hashrf.<br>
We embed the distances in a 3-dimensional space using PCoA, and we plot the results colouring by STEP and highlighting some specific trees for each set. We also display some buttons to hide/show the sets in the plots, we use different colorscales (the default behaviour), and we show the results.


```python
!cat example_1.toml
```





























```python
!pear_ebi --config example_1.toml
```

    PEAR v0.1.85
    Looking into directory ../beast_trees/ - pattern: *run2*
    Your input:
    ─────────────────────────────
     Tree set collection containing 3003 trees;
     File: Set_collection_72c842e8-e1ac-4a90-8fbc-535c4b10ef92;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run1; Containing 1001 trees.
    beast_long; Containing 1001 trees.
    beast_run2; Containing 1001 trees.

    ⠸ Calculating distances...0m
    hashrf_RF | Done!
    ⠋ Embedding distances...0m
    pcoa | Done!

    - Leaving PEAR -


#### Example 1 Continued
We now show that we can easily modify a single parameter without touching the `.toml` file by simply overriding it on the command line. <br>
As an example, we change the embedding method to tSNE.


```python
!pear_ebi --config example_1.toml --tsne 3
```

    PEAR v0.1.85
    Looking into directory ../beast_trees/ - pattern: *run2*
    Your input:
    ─────────────────────────────
     Tree set collection containing 3003 trees;
     File: Set_collection_9a55b5c8-53b6-40c3-b4c7-55df883b3258;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run1; Containing 1001 trees.
    beast_long; Containing 1001 trees.
    beast_run2; Containing 1001 trees.

    ⠼ Calculating distances...0m
    hashrf_RF | Done!
    ⠹ Embedding distances...0m
    tsne | Done!

    - Leaving PEAR -


### Example 2
We use pear to compare 5 different algorithms used to estimate a phylogenetic tree structure with the real structure (we know it as we use simulated data).<br>
We upload 6 files using the `[trees]` argument.<br>
We compute the Robinson Foulds distances using hashrf.<br>
We embed the distances in a 3-dimensional space using PCoA, and we plot the results colouring by the Likelihood values obtained during the runs. We highlight the true tree structure.<br>
We substitute the 3$^{rd}$ dimension with the Likelihood scores, and we use the same colorscale to represent the likelihood of the proposed structures.


```python
!cat example_2.toml
```





























```python
!pear_ebi --config example_2.toml
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 138 trees;
     File: Set_collection_ec4084c4-a5a7-4878-abcc-93c634208b6e;
     Distance matrix: computed.
    ─────────────────────────────
    IQtreeStartingTree_Trees; Containing 29 trees.
    MapleStartingTree_Trees; Containing 5 trees.
    ParsimonyRAxMLStartingTree_GTRmodel_Trees; Containing 47 trees.
    RAxMLNGStartingTree_Trees; Containing 26 trees.
    UshERStartingTree_Trees; Containing 30 trees.
    TrueTreeSimulations; Containing 1 trees.

    ⠋ Embedding distances...0m
    pcoa | Done!

    - Leaving PEAR -


#### Example 2 Continued
We now show that we can easily run the same analyses in a very neat way by simply renaming our `.toml` file by simply renaming it `pear.toml`.


```python
!pear_ebi
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 138 trees;
     File: Set_collection_b87402e6-2e41-4902-8fd2-092a3f6c230a;
     Distance matrix: computed.
    ─────────────────────────────
    IQtreeStartingTree_Trees; Containing 29 trees.
    MapleStartingTree_Trees; Containing 5 trees.
    ParsimonyRAxMLStartingTree_GTRmodel_Trees; Containing 47 trees.
    RAxMLNGStartingTree_Trees; Containing 26 trees.
    UshERStartingTree_Trees; Containing 30 trees.
    TrueTreeSimulations; Containing 1 trees.

    ⠋ Embedding distances...0m
    pcoa | Done!

    - Leaving PEAR -
