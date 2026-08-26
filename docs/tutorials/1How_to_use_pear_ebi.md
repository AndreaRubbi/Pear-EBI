# 1. How to use Pear
This tutorial goes through the basic features of pear_ebi and how to use them
_____
## Setup ##
First of all, check that you have a supported python version (python 3.10 to 3.13). We strongly encourage the creation of a dedicated virtual environment in order to avoid potential conflicts with other libraries due to the mismatch of dependencies' versions. We also support the use of mamba as a more efficient version of conda.

<b>In your terminal:</b><br><br>
Download and install <b>Miniforge</b>, which provides mamba and uses conda-forge by default: <br>
<a href="https://github.com/conda-forge/miniforge#install">See documentation</a>
<br>Create new environment with one of the supported versions of python:<br>
```mamba create -n pear_env -c conda-forge python=3.11```
<br>Activate environment:<br>
```mamba activate pear_env```
<br>Install pear:<br>
```python -m pip install pear_ebi```

<font color='red'><b>Note on channels:</b></font> the <code>-c conda-forge</code> flag above is deliberate. Some institutes (EMBL-EBI among them) are no longer licensed to use Anaconda's <code>defaults</code> channel, so packages must come from conda-forge. Miniforge is configured that way out of the box; if you use Miniconda or Anaconda instead, keep passing <code>-c conda-forge</code> explicitly.

<br><font color='red'><b>Alternatively:</b></font><br>
<br>Create new environment with one of the supported versions of python:<br>
```python3.11 -m venv pear_env```
<br>Activate environment:<br>
```source pear_env/bin/activate```
<br>Install pear:<br>
```python -m pip install pear_ebi```

### Optional ###
If you're planning on performing more advanced analyses, such as the ones described in the <b><i>"Advanced Examples"</b></i>, install pear with the <code>notebook</code> extra, which adds Jupyter and JupyterLab:<br>
```python -m pip install "pear_ebi[notebook]"```
<br>and also install the new jupyter kernel:<br>
```python -m ipykernel install --user --name=pear_ebi```

### Working from a checkout ###
If you cloned the repository rather than installing from PyPI, dependencies are managed with <a href="https://python-poetry.org/">Poetry</a> and the exact versions are recorded in <code>poetry.lock</code>:<br>
```poetry install```
<br>which reproduces a known-good environment. Add <code>--with docs</code> to also install the documentation toolchain.
****

## Basic Use ##
After following the steps above to set up pear_ebi, you should be ready to use all the features of pear_ebi!
<br>This notebook is a good guide to learn how to use it and to check that your installation is successful. If you should have any problem, please contact us by filing an issue on github.<br>
To start with, simply check your installation by running:


```python
!pear_ebi
```

    PEAR v0.1.85
    No files specified (see --help for instructions)
    - Leaving PEAR -


Pear is complaining because no file was given... since it looks like you don't know how to use it, it kindly suggests to seek help using the `--help` (or simply `-h`) flag. Good idea!


```python
!pear_ebi -h
```

    usage: PEAR [-h] [-o output] [--interactive] [-d distance_matrix]
                [--meta metadata] [-m METHOD] [--pcoa PCOA] [--tsne TSNE] [--plot]
                [--config CONFIG] [--quality] [--dir DIR] [--pattern PATTERN]
                [input ...]

    PEAR-EBI v0.1.85 | Phylogeny Embedding and Approximate Representation
    Calculates Robinson-Foulds distances between large set of trees

    positional arguments:
      input                 input file : tree set in newic format

    optional arguments:
      -h, --help            show this help message and exit
      -o output             output file : storage of distance matrix
      --interactive, -i     run the program in interactive mode - only the input
                            file, distance matrix, output file, and metadata
                            arguments will be considered
      -d distance_matrix, --dM distance_matrix
                            distance matrix : specify file containing a
                            precomputed distance matrix
      --meta metadata       metadata : csv file with metadata for each tree
      -m METHOD, --method METHOD
                            calculates tree distances using specified method
                            (hashrf_RF, hashrf_wRF, smart_RF, tqdist_quartet,
                            tqdist_triplet)
      --pcoa PCOA           embedding using PCoA: select number of components
                            (int) to be calculated
      --tsne TSNE           embedding using t-SNE: select number of components
                            (int) to be calculated
      --plot, -p            plot embedding in 2 or 3 dimensions
      --config CONFIG, -c CONFIG
                            toml config file
      --quality, -q         assess quality of embedding
      --dir DIR             directory with files
      --pattern PATTERN     pattern of files in directory

    Author: Andrea Rubbi - Goldman Group | European Bioinformatics Institute


Essentially, pear_ebi does the following:<ul>
    <li>Reads phylogenetic trees in newick format from one or multiple files;
    <li>Computes the distances between the trees with one of the available metrics (refer to `--help` to see them and to the manuscript for their specific functions);
    <li>Embeds the distances in a lower dimensional space. Basic usage allows using either PCoA or tSNE;
    <li>Plots the resulting embeddings.
 </ul>
 Why should you use pear_ebi? Because it's <font color='green'><b>simple</b><font color='black'>, <font color='dark orange'><b>fast</b><font color='black'>, and <font color='purple'><b>produces nice graphs</b></font>!

### Load trees ###
Pear has two kind of structures for your trees:<ol>
    <li><font color='blue'>tree_set</font>;<li><font color='dark pink'>set_collection</font>.
</ol>
You define a <font color='blue'>tree_set</font> with a set of trees coming from a single file, whereas a <font color='dark pink'>set_collection</font> is composed of trees coming from multiple files, divided into groups depending on the file of origin. In practice, computationally and conceptually, a <font color='dark pink'>set_collection</font> is a set of multiple <font color='blue'>tree_set</font>.
<br>We'll make it clear exploring how we can load trees on pear with a few examples.


```python
# load single set of trees into a tree_set #
!pear_ebi beast_trees/beast_long.trees
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 1001 trees;
     File: beast_trees/beast_long.trees;
     Distance matrix: not computed.
    ─────────────────────────────


    - Leaving PEAR -



```python
# load two set of trees into a set_collection #
!pear_ebi beast_trees/beast_run1.trees beast_trees/beast_run2.trees
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 2002 trees;
     File: Set_collection_1bfec55f-592e-43a6-bf2b-4000219b3092;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run1; Containing 1001 trees.
    beast_run2; Containing 1001 trees.


    - Leaving PEAR -


Although you may be perplexed about the actual practical difference between these structures, we suggest you keep these questions for the next chapters and bear with us to explore some more funky ways of loading trees onto pear.
<br>In fact, you can avoid boring repetitions of the path of your files and simply indicate the directory!


```python
# load files from directory #
!pear_ebi --dir beast_trees
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 9009 trees;
     File: Set_collection_696e3bb7-7deb-46e1-97c8-3bcce630cddd;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run7; Containing 1001 trees.
    beast_run2; Containing 1001 trees.
    beast_long; Containing 1001 trees.
    beast_run5; Containing 1001 trees.
    beast_run6; Containing 1001 trees.
    beast_run3; Containing 1001 trees.
    beast_run8; Containing 1001 trees.
    beast_run4; Containing 1001 trees.
    beast_run1; Containing 1001 trees.


    - Leaving PEAR -


Whoa! Too many... thankfully we can indicate a pattern to look for


```python
# load files from directory using pattern#
!pear_ebi --dir beast_trees --pattern "*run*"
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 8008 trees;
     File: Set_collection_9781a60f-d6d7-4ffb-97c3-36d7aed9af4b;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run7; Containing 1001 trees.
    beast_run2; Containing 1001 trees.
    beast_run5; Containing 1001 trees.
    beast_run6; Containing 1001 trees.
    beast_run3; Containing 1001 trees.
    beast_run8; Containing 1001 trees.
    beast_run4; Containing 1001 trees.
    beast_run1; Containing 1001 trees.


    - Leaving PEAR -


If you are a regex wizard you can probably select any set of similarly-named files in this way. Me? I generally ask GPT to write the magic formula. <b>NB:</b> you can also use a combination of `--dir`, `--pattern`, and normal file definition if one or more files come from other directories.


```python
!pear_ebi beast_trees/beast_long.trees --dir beast_trees --pattern "*run[1,2]*"
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 3003 trees;
     File: Set_collection_21d56458-1d6b-4dab-9842-f8a65287e845;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run2; Containing 1001 trees.
    beast_run1; Containing 1001 trees.
    beast_long; Containing 1001 trees.


    - Leaving PEAR -


### Compute Distances ###
You can compute the distance matrix using different methods. Each method has a specific purpose, which is outlined in the associated paper. However, in general, the Robinson Foulds distance metric is an admittedly good choice. Additionally, pear_ebi computes this metric using the hashrf algorithm, which is the fastest way of computing such metric to date. You can indicate any method available using the `-m`, or `--method`, flag. When a metric/method is indicated, pear will use it to compute the distance matrix <b>even if a distance matrix is given</b>, and it <b>will overwrite any previous matrix saved at a file with the same standard format</b> (realistically, this happens only if the matrix was produced using pear during an interactive or advanced session - see the interactive sessions chapter).


```python
# compute Robinson Foulds distances using hashrf #
!pear_ebi beast_trees/beast_run1.trees -m hashrf_RF
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 1001 trees;
     File: beast_trees/beast_run1.trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠙ Calculating distances...0m
    hashrf_RF | Done!

    - Leaving PEAR -



```python
# compute Robinson Foulds distances using hashrf #
!pear_ebi --dir beast_trees --pattern "*run[12]*" --method hashrf_RF
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set collection containing 2002 trees;
     File: Set_collection_c9cf5472-d392-4825-b88f-e9feffdb331a;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run2; Containing 1001 trees.
    beast_run1; Containing 1001 trees.

    ⠙ Calculating distances...0m
    hashrf_RF | Done!

    - Leaving PEAR -



```python
# compute quartet distances #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m tqdist_quartet
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠴ Calculating distances...0m
    tqdist_quartet | Done!

    - Leaving PEAR -



```python
# compute modified RF distances #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m smart_RF
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠙ Calculating distances...0m
    smart_RF | Done!

    - Leaving PEAR -


If you run the exaples above you will notice the difference between the time efficiency of hashrf and any other algorithm. Please note the astonishing performance of that algorithm, especially considering the sensible difference in the number of trees analyzed!

Given that one may want to analyze the distance matrix in multiple ways,
possibly desiring to skip the hefty distance-computation step, we introduced the
convenient `-d` flag that allows specifying a previosuly computed distance matrix.
<b>Please note that</b> the order of the trees in the distance matrix is preserved and
thus one must be consistent with the input specification when reusing a distance matrix.
Should you disregard this suggestion, we suggest you also disregard your downstream analyses
(and perhaps everything else as well).


```python
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m hashrf_RF -o precomputed_distance_matrix
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠴ Calculating distances...0m
    hashrf_RF | Done!

    - Leaving PEAR -



```python
# reusing distances #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -d precomputed_distance_matrix
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: computed.
    ─────────────────────────────


    - Leaving PEAR -


Please note the difference in the input "Distance matrix" status.

### Embed distances ###
Here we show the easy way of embedding the distances and produce nice plots! First of all, one may choose any number of dimensions $M \leq N$; where $N$ is the number of trees and $M$ the chosen dimensionality of the embedded space. Note that, should $M$ be $\lt 3$, pear will produce only one 2D plot. Whereas, should $M$ be $\geq 3$ pear will produce a 2D and a 3D plot. Please note that we are currently unable to produce human-friendly representations for $\gt 3$D data... In any case, whichever dimension $M$ may be chosen, an $M$-dimensional embedding shall be produced and then saved on the machine.


```python
# compute triplet distances, embed in 2D using pcoa #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m tqdist_triplet --pcoa 2
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠹ Calculating distances...0m
    tqdist_triplet | Done!
    ⠋ Embedding distances...
    pcoa | Done!

    - Leaving PEAR -



```python
# compute RF distances, embed in 5D using tsne #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m hashrf_RF --tsne 5
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠦ Calculating distances...0m
    hashrf_RF | Done!
    ⠸ Embedding distances...0m
    tsne | Done!

    - Leaving PEAR -


Use the `--plot` flag to show the plots at the end. Note that, when an embedding method is specified, the plots are produced regardless of whether the `--plot` flag is present or not.


```python
# compute RF distances, embed in 5D using tsne #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m hashrf_RF --tsne 5 --plot
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠴ Calculating distances...0m
    hashrf_RF | Done!
    ⠸ Embedding distances...0m
    tsne | Done!

    - Leaving PEAR -


<b>NB:</b> this is something you are expected to run in a terminal, this is why the graph doesn't show up here.

The `--quality`, or `-q`, flag indicates to pear to provide some quality metrics for the embedding generated. <br>
The toolset of quality-measures provided by pear may be expanded upon request in the future.


```python
# compute RF distances, embed in 5D using tsne, computing quality #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m hashrf_RF --tsne 5 --quality
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠼ Calculating distances...0m
    hashrf_RF | Done!
    ⠼ Embedding distances...0m
    tsne | Done!
    With 5 components/dimensions, the estimated correlation with the 32-dimensional
    coordinates is -0.03

    - Leaving PEAR -



```python
# compute RF distances, embed in 2D using pcoa, computing quality #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m hashrf_RF --pcoa 2 --quality
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠼ Calculating distances...0m
    hashrf_RF | Done!
    ⠋ Embedding distances...
    pcoa | Done!
    With 2 components/dimensions, the explained variance is 96.56,
     with an estimated correlation 1.00 with the 32-dimensional coordinates

    - Leaving PEAR -



```python
# compute RF distances, embed in 10D using pcoa, computing quality #
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -m hashrf_RF --pcoa 10 --quality
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠦ Calculating distances...0m
    hashrf_RF | Done!
    ⠋ Embedding distances...
    pcoa | Done!
    With 10 components/dimensions, the explained variance is 99.98,
     with an estimated correlation 1.00 with the 32-dimensional coordinates

    - Leaving PEAR -


<i>Pro tip</i>: since in order to compute an embedded space we need the distance matrix first, pear will automatically compute the distances with `hashrf_RF` when no method is indicated (and no precomputed distance matrix is specified)


```python
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees --pcoa 2
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    ⠦ Calculating distances...0m
    hashrf_RF | Done!
    ⠋ Embedding distances...
    pcoa | Done!

    - Leaving PEAR -


## <font color='blue'>tree_set</font> vs <font color='dark pink'>set_collection</font> ##
Pear has two kind of structures for your trees:<ol>
    <li><font color='blue'>tree_set</font>;<li><font color='dark pink'>set_collection</font>.
</ol>
<font color='blue'>tree_set</font> allows for analyzing a single set of trees.
<font color='dark pink'>set_collection</font> allows for analyzing multiple set of trees and the divergences between them.
In our interactive plots, <font color='dark pink'>set_collection</font> generates different sets of points that are colored
according to the file of origin. This could be useful, for instance, if one performs multiple runs of a phylogenetic-tree estimating algorithm and wants to assess the robustness/consistency of that method.
<br>The distances computation is overall the same as a square distance matrix is generated, encompassing the whole collection of trees.

****

## Interactive Sessions ##
All the above, but staying in the loop!
<br>You can add set of trees iteratively and compute distances/embeddings as many times as you want to.


```python
!pear_ebi -i
```

    PEAR v0.1.85
    Specify file with tree set

    File: ^C

    - Leaving PEAR -



```python
!pear_ebi MAPLE_res/IQtreeStartingTree_slower_Trees -i
```

    PEAR v0.1.85
    Your input:
    ─────────────────────────────
     Tree set containing 32 trees;
     File: MAPLE_res/IQtreeStartingTree_slower_Trees;
     Distance matrix: not computed.
    ─────────────────────────────

    PEAR | Interactive Mode
    ⣿⣿⣿⣿⣿⣿⣿⣿⠿⣟⠉⡿⠿⣿⣿⣿⣿⣿⣿⣿ -- Controls --
    ⣿⣿⣿⣿⡿⣿⢉⢳⠴⣞⠉⡷⢥⡏⡙⡿⢿⣿⣿⣿ 1 --> see status
    ⣿⣿⡋⢻⡤⣼⠉⢯⡆⣞⠙⣧⣢⠏⠪⣣⢦⡛⠹⣿ 2 --> calculate distances
    ⣿⣿⠓⢻⣄⣼⠋⢷⡠⡽⠚⣉⣤⡞⢚⢦⢢⠟⠹⣿ 3 --> embed distances
    ⣿⣿⢓⢻⡄⡼⠗⢃⣂⡒⠻⣧⣂⡿⠚⣨⢨⡓⠻⣿ 4 --> plot embeddings
    ⣿⣿⠗⢎⢄⠂⠾⣯⣂⡽⢓⢆⢔⠐⠿⣇⢅⡗⠻⣿ 5 --> add set to collection
    ⣿⣿⠖⢯⡡⣹⠗⣤⣉⠛⠶⣏⢌⡿⠲⡌⢌⢞⠼⣿ 6 --> get subset
    ⣿⣿⣮⣾⡉⣹⠦⣞⡉⡽⠶⡌⠌⡮⠲⣏⢩⣷⣼⣿
    ⣿⣿⣿⣿⣿⣿⣤⣞⢉⣳⠥⣏⠍⣧⣵⣿⣿⣿⣿⣿ 7 --> exit
    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣤⣿⣿⣿⣿⣿⣿⣿⣿⣿ 8 --> see list of controls
    Command:

## There is more...
the flags `--meta` and `--config` introduce another layer of flexibility into our analyses.
<br> Please, refer to the <b><i>"Advanced Examples"</b></i> folder to get a gist of the full potential of pear!
