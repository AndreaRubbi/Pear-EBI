# 3. How to use Pear on python
<h1></h1>

Here we report the code used to generate the examples in the manuscript, providing a general outline of the use of the python library

`tree_set` is the main module of `pear_ebi` - it contains both the <font color="blue"><b>tree_set</b></font> and the <font color="red"><b>set_collection</b></font> objects.


```python
from pear_ebi import tree_set
import numpy as np
import pandas as pd
import time
import os
```

## Examples
Here we do not delve into the specifics or the rationale of the experiments - for that we redirect you to the related manuscript - but rather focus on the use of the library to analyze sets of trees in different settings.

***
# BEAST
### <font color='purple'> Bayesian - Markov Chain Monte Carlo </font>

This is a relatively simple example: we have a few files containing trees in Newick format, where each set of trees is produced by a program that produced them sequentially. We want to represent the distribution of these trees in order to analyze the single trajectories and to compare the different trjectories coming from different runs. <br>
With this idea in mind, we simply compute the square distance matrix using the `Robinson Foulds` metric, which generally represents the relations between trees in an effective way, using the most efficient algorithm: `hashrf`.<br>
We then compute and plot the `PCoA` embeddings in 2 and 3 dimensions.

#### Step 1: load trees

<font color="blue"><b>tree_set</b></font> is meant to store and manage the analysis of a single set of trees stored in the same file. These trees have to be in the (standard) Newick format.<br>
`tree_set = tree_set.tree_set(file, output_file=None, distance_matrix=None, metadata=None)`
<br>The first argument indicates the file where the trees are stored. `output_file` can be empty and may be used to indicate a specific name and path for the output distance matrix. `distance_matrix` can be used to indicate a precomputed distance matrix in the `.csv` format. `metadata` may be used to indicate a `.csv` file with additional information relative to the trees. It must be of the same dimensionality of the set of trees, hence $size = (|trees|, |meta$-$variables|)$.
 ***
<font color="red"><b>set_collection</b></font> performs the same tasks of <font color="blue"><b>tree_set</b></font>, but stores multiple set of trees. <b>NB:</b> for each file, a set of trees is defined within the set collection. In practice, a <font color="red"><b>set_collection</b></font> is a collection of different <font color="blue"><b>tree_set</b></font>, each one containing a set of trees coming from a different file. It should not come as a surprise, then, that the input `collection=list()`, potentially empty (one may initialize an empty collection), is a list of files containing trees in the Newick format.
`set_collection=tree_set.set_collection(collection=list(), file="Set_collection_", output_file=None, distance_matrix=None, metadata=None,)`
<br>`file` may be used to indicate an alternative name for the file containing the collection information. This file will be used by pear to store information relative to the sets of trees. `output_file` may be used to indicate a specific name for the distance matrix file.
<br>`distance_matrix` can be used to indicate a precomputed distance matrix in the `.csv` format. `metadata` may be used to indicate a `.csv` file with additional information relative to the trees. It must be of the same dimensionality of the set of trees, hence $size = (|trees|, |meta$-$variables|)$.


```python
beast_dir = '../beast_trees'
beast_runs = [os.path.join(beast_dir, f"beast_run{i}.trees") for i in range(1,9)] + [os.path.join(beast_dir, f"beast_long.trees")]
beast_collection = tree_set.set_collection(beast_runs)
print(beast_collection)
```

    ─────────────────────────────
     Tree set collection containing 9009 trees;
     File: Set_collection_d67b3fc3-bac4-4e5f-934c-0c8dc668f012;
     Distance matrix: not computed.
    ─────────────────────────────
    beast_run1; Containing 1001 trees.
    beast_run2; Containing 1001 trees.
    beast_run3; Containing 1001 trees.
    beast_run4; Containing 1001 trees.
    beast_run5; Containing 1001 trees.
    beast_run6; Containing 1001 trees.
    beast_run7; Containing 1001 trees.
    beast_run8; Containing 1001 trees.
    beast_long; Containing 1001 trees.



#### Step 2: compute distances


```python
beast_collection.calculate_distances(method="hashrf_RF")
# where method can be chosen among hashrf_RF,
# hashrf_wRF, smart_RF, tqdist_quartet, tqdist_triplet
```


    Output()



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080; font-weight: bold">hashrf_RF | Done!</span>
</pre>



#### Step 3: compute embeddings


```python
beast_collection.embed(method="pcoa", dimensions=3, quality=False)
# where method can be chosen among pcoa, tsne, isomap, lle
# pro tip: if the distance matrix has not been computed prior to calling this function,
# it will automatically be computed using hashrf_RF
```


    Output()



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080; font-weight: bold">pcoa | Done!</span>
</pre>



#### Step 4: Plot embeddings


```python
beast_collection.plot_2D(
        method='pcoa', # method used for the embedding - if not previously computed, it will be computed when calling this function
        save=False, # save the plot in a pdf format? The .html will be saved anyway!
        name_plot=None, # specific name for the plot
        static=False, # create a static element rather than a dynamic widget
        plot_meta="SET-ID", # meta-variable used to colour the points (trees)
        plot_set=None, # select specific set of sets of trees in the collection for the plot
        select=False, # create widgets to hide/display tree_sets in the graph
        same_scale=False,) # use the same colorscale for each tree_set (useful when the same metric is compared among sets)
```


```python
beast_collection.plot_3D(
        method='pcoa',
        save=False,
        name_plot=None,
        static=False,
        plot_meta="SET-ID",
        plot_set=None,
        select=False,
        same_scale=False,
        z_axis=None,) # substitute the 3rd axis (Z-axis/3rd PCoA) with a meta-variable of choice
```

***
# RAxML
### <font color='purple'> Bootstrap Analysis </font>

We do not delve into the biological interpretation of this specific example, as we discussed that in our manuscript. <br> Instead, we point your attention to some specific features of pear's interactive representations: <ul>
    <li> First: you can directly plot the `PCoA` embedding of the Robinson Foulds distances computed with `hashrf` without explicitly going through these steps. In fact, running a plotting function bypasses steps 2 and 3, performing them "under the hood";
    <li> In the 2D representation:<ul>
        <li> From the dropdown menu, one can choose the meta-vatiable used to color the points in the graphs;
        <li> The red button "Save plot as pdf" saves the plot in a `.pdf` file;
        <li> The $#show `name of the tree set`$ buttons allow showing and hiding specific sets of trees;
        <li> One may choose between the scatter representation or a contour plot;
        <li> While visualizing a contour plot, one may choose to hide or show the points (trees) with the light blue button "Show points on a Contour plot";
        <li> One can choose between visualizing only the points or also the connections between them (indicating the progress through an hypotetical sequential order);
        <li> One can use the native tools of the plot so zoom in to focus on a specific area of the graph;
        <li> If one wants to focus their attention on specific sets of trees, they just have to click them: this will highlight the selected set and make the other sets transparent. One may continue and highlight on other sets by simply clicking on them. In order to go back to the original representation, one just needs to click on one of the previously highlighted sets for a second time.
    </ul>
   <li> In the 3D representation:<ul>
       <li> Same features of the 2D representation. However, there is no contour plot in this case.


```python
input_dir = '../bootstrap_mammals/'
files = ([f'bootstrap_{i}' for i in [6086, 5261, 5092, 281, 11289, 10409]] +
    [f'bestTree_{i}' for i in [6086, 5261, 5092, 281, 11289, 10409]])
collection = tree_set.set_collection([os.path.join(input_dir, f) for f in files])
print(collection)
```

    ─────────────────────────────
     Tree set collection containing 3656 trees;
     File: Set_collection_4cbb9500-1e98-40cd-84ec-a9c8bd511698;
     Distance matrix: not computed.
    ─────────────────────────────
    bootstrap_6086; Containing 500 trees.
    bootstrap_5261; Containing 600 trees.
    bootstrap_5092; Containing 800 trees.
    bootstrap_281; Containing 600 trees.
    bootstrap_11289; Containing 600 trees.
    bootstrap_10409; Containing 550 trees.
    bestTree_6086; Containing 1 trees.
    bestTree_5261; Containing 1 trees.
    bestTree_5092; Containing 1 trees.
    bestTree_281; Containing 1 trees.
    bestTree_11289; Containing 1 trees.
    bestTree_10409; Containing 1 trees.




```python
collection.plot_2D('pcoa', select=True)
```


```python
collection.plot_3D('pcoa', select = True)
```


<h1> MAPLE </h1>
<h3> <font color='purple'> Maximum Likelihood Search </font></h3>

In this specific example, we focus the readers' attention on the possibility of specifying a metadata file to add meta-variable and use them to color the points in the plots or to replace the $Z_{axis}$ in the 3D plot.<br>
In the cell below we specify a list of files containing trees and a list of files containing the likelihood of each tree.


```python
maple_tree = [
    "IQtreeStartingTree_Trees",
    "MapleStartingTree_Trees",
    "ParsimonyRAxMLStartingTree_GTRmodel_Trees",
    "RAxMLNGStartingTree_Trees",
    "UshERStartingTree_Trees",
    "TrueTreeSimulations",
    ]

maple_LK = [
    "IQtreeStartingTree_LK",
    "MapleStartingTree_LK",
    "ParsimonyRAxMLStartingTree_GTRmodel_LK",
    "RAxMLNGStartingTree_LK",
    "UshERStartingTree_LK",
]

maple_dir = '../MAPLE_res/'
```

In the cell below we check that "DisMat_Maple.csv", the precomputed distance matrix, is present in the folder. We precomputed the distances as the `smart_RF` algorithm, which computes a modified version of the Robinson Foulds metric, takes much longer than `hashrf_RF`.


```python
set_list = [os.path.join(maple_dir, tree) for tree in maple_tree]
try:collection_maple = tree_set.set_collection(set_list, distance_matrix = "DistMat_Maple.csv")
except:collection_maple = tree_set.set_collection(set_list)
finally:print(collection_maple)
```

    ─────────────────────────────
     Tree set collection containing 138 trees;
     File: Set_collection_17e30756-0754-4fbb-9ed7-cb862716a7f8;
     Distance matrix: computed.
    ─────────────────────────────
    IQtreeStartingTree_Trees; Containing 29 trees.
    MapleStartingTree_Trees; Containing 5 trees.
    ParsimonyRAxMLStartingTree_GTRmodel_Trees; Containing 47 trees.
    RAxMLNGStartingTree_Trees; Containing 26 trees.
    UshERStartingTree_Trees; Containing 30 trees.
    TrueTreeSimulations; Containing 1 trees.



This is how we extracted the information from multiple files - one may simply skip this step as we then saved the information in the "Likelihoods.csv" file.


```python
LKs = dict()
for lk_file in maple_LK:
    file = os.path.join(maple_dir, lk_file)
    with open(file, 'r') as f:
        LKs[lk_file] = np.array(f.readlines())
        f.close()

LK = list()
for lk in LKs.values(): LK.extend(lk)

LK = np.array([tree.replace('\n', '') for tree in LK], dtype=np.float64)
LK = np.concatenate([LK,np.array([-257513])], axis = 0)
LK = np.interp(LK, (LK.min(), LK.max()), (0, +1)) # scale LKs between 0 and 1

df_LK = pd.DataFrame({'LK':LK})
df_LK.to_csv("Likelihoods.csv")
df_LK
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>LK</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.983439</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.983822</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.984404</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.985124</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.985631</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th>133</th>
      <td>0.999844</td>
    </tr>
    <tr>
      <th>134</th>
      <td>0.999850</td>
    </tr>
    <tr>
      <th>135</th>
      <td>0.999976</td>
    </tr>
    <tr>
      <th>136</th>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>137</th>
      <td>0.706890</td>
    </tr>
  </tbody>
</table>
<p>138 rows × 1 columns</p>
</div>



There are two ways to introduce new variables in a <font color="red"><b>set_collection</b></font>'s or <font color="blue"><b>tree_set</b></font>'s metadata:


```python
# Method 1:
# Given that metadata is a pandas dataframe,
# one may simply add columns to it!
df_LK = pd.read_csv("Likelihoods.csv")
collection_maple.metadata['LK'] = df_LK['LK']
collection_maple.metadata
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SET-ID</th>
      <th>STEP</th>
      <th>LK</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>0</td>
      <td>0.983439</td>
    </tr>
    <tr>
      <th>1</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>1</td>
      <td>0.983822</td>
    </tr>
    <tr>
      <th>2</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>2</td>
      <td>0.984404</td>
    </tr>
    <tr>
      <th>3</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>3</td>
      <td>0.985124</td>
    </tr>
    <tr>
      <th>4</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>4</td>
      <td>0.985631</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>133</th>
      <td>UshERStartingTree_Trees</td>
      <td>26</td>
      <td>0.999844</td>
    </tr>
    <tr>
      <th>134</th>
      <td>UshERStartingTree_Trees</td>
      <td>27</td>
      <td>0.999850</td>
    </tr>
    <tr>
      <th>135</th>
      <td>UshERStartingTree_Trees</td>
      <td>28</td>
      <td>0.999976</td>
    </tr>
    <tr>
      <th>136</th>
      <td>UshERStartingTree_Trees</td>
      <td>29</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>137</th>
      <td>TrueTreeSimulations</td>
      <td>0</td>
      <td>0.706890</td>
    </tr>
  </tbody>
</table>
<p>138 rows × 3 columns</p>
</div>




```python
# Method 2:
# One may specify a file with additional
# metadata when the object is created
collection_maple = tree_set.set_collection(set_list, distance_matrix = "DistMat_Maple.csv", metadata="Likelihoods.csv")
collection_maple.metadata
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SET-ID</th>
      <th>STEP</th>
      <th>LK</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>0</td>
      <td>0.983439</td>
    </tr>
    <tr>
      <th>1</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>1</td>
      <td>0.983822</td>
    </tr>
    <tr>
      <th>2</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>2</td>
      <td>0.984404</td>
    </tr>
    <tr>
      <th>3</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>3</td>
      <td>0.985124</td>
    </tr>
    <tr>
      <th>4</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>4</td>
      <td>0.985631</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>133</th>
      <td>UshERStartingTree_Trees</td>
      <td>26</td>
      <td>0.999844</td>
    </tr>
    <tr>
      <th>134</th>
      <td>UshERStartingTree_Trees</td>
      <td>27</td>
      <td>0.999850</td>
    </tr>
    <tr>
      <th>135</th>
      <td>UshERStartingTree_Trees</td>
      <td>28</td>
      <td>0.999976</td>
    </tr>
    <tr>
      <th>136</th>
      <td>UshERStartingTree_Trees</td>
      <td>29</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>137</th>
      <td>TrueTreeSimulations</td>
      <td>0</td>
      <td>0.706890</td>
    </tr>
  </tbody>
</table>
<p>138 rows × 3 columns</p>
</div>



If, for any reason, the "DisMat_Maple.csv" is not present in the folder, we need to compute the distance matrix. We chose to upload the precomputed distance matrix as the `smart_RF` method takes much longer to compute than `hashrf_RF`, and we wanted this notebook to run smoothly and rapidly.


```python
if collection_maple.distance_matrix is None:
    start = time.time()
    collection_maple.calculate_distances('days_RF')
    np.savetxt("DistMat_Maple.csv", collection_maple.distance_matrix, delimiter=",")
    print(time.time() - start)
collection_maple.distance_matrix.shape
```




    (138, 138)



We colour the trees based on their likelihood in the 2D plot, while we replace the $Z_{axis}$ with the likelihood in the 3D plot. <br>
`same_scale` makes sure that points with the same value have the same colour in the graphs.


```python
collection_maple.plot_3D('pcoa', plot_meta = 'LK', same_scale = True, select = True)
```


```python
collection_maple.plot_3D('pcoa', plot_meta = 'LK', same_scale = True, select = True, z_axis = 'LK')
```

What if I wanted to make one of these trees to "pop out"? Well, we can add a coulumn called `"highlight"` that will automatically be read by the plotting function to "highloght" the points specified. The columns has to be binary (0s and 1s), where the 1s indicate that a tree should be highlighted. It is crucial here to know the order and size of the sets, as the column is shared by the whole collection. <b>Note that</b>: the vector `"highlight"` <b>must be integer (dtype = int)!</b>.

<h4> Easy example:</h4>
we want to highight the true tree, which is the last in the collection.


```python
highlight_mask = np.zeros(collection_maple.metadata.shape[0], dtype=int) # vector of 0s with size = n_trees
highlight_mask[-1] = 1 # last element set to 1
collection_maple.metadata['highlight'] = highlight_mask
collection_maple.metadata
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SET-ID</th>
      <th>STEP</th>
      <th>LK</th>
      <th>highlight</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>0</td>
      <td>0.983439</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>1</td>
      <td>0.983822</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>2</td>
      <td>0.984404</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>3</td>
      <td>0.985124</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>4</td>
      <td>0.985631</td>
      <td>0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>133</th>
      <td>UshERStartingTree_Trees</td>
      <td>26</td>
      <td>0.999844</td>
      <td>0</td>
    </tr>
    <tr>
      <th>134</th>
      <td>UshERStartingTree_Trees</td>
      <td>27</td>
      <td>0.999850</td>
      <td>0</td>
    </tr>
    <tr>
      <th>135</th>
      <td>UshERStartingTree_Trees</td>
      <td>28</td>
      <td>0.999976</td>
      <td>0</td>
    </tr>
    <tr>
      <th>136</th>
      <td>UshERStartingTree_Trees</td>
      <td>29</td>
      <td>1.000000</td>
      <td>0</td>
    </tr>
    <tr>
      <th>137</th>
      <td>TrueTreeSimulations</td>
      <td>0</td>
      <td>0.706890</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
<p>138 rows × 4 columns</p>
</div>




```python
collection_maple.plot_2D('pcoa', plot_meta = 'LK', same_scale = True, select = False, static = True)
```

<h4> Hard (not really) example:</h4>
we want to highight the true tree, and the last tree of each set.
In order to perform this task, we can exploit an additional tool: the `set_collection.data`.
This object contains some useful information regarding the structure of each set composing the collection. In particular, we may be interested in the `n_trees` variable which can be used to assess the size of each set!


```python
# this is the object --> collection_maple.data # Note that this is a nested dictionary
for Tset, info in collection_maple.data.items():
    print(Tset, f"has {info['n_trees']} trees")
```

    IQtreeStartingTree_Trees has 29 trees
    MapleStartingTree_Trees has 5 trees
    ParsimonyRAxMLStartingTree_GTRmodel_Trees has 47 trees
    RAxMLNGStartingTree_Trees has 26 trees
    UshERStartingTree_Trees has 30 trees
    TrueTreeSimulations has 1 trees



```python
highlight_mask = np.zeros(collection_maple.metadata.shape[0], dtype=int)
last_tree = 0
for Tset, info in collection_maple.data.items():
    last_tree += info["n_trees"]
    highlight_mask[last_tree - 1] = 1
collection_maple.metadata['highlight'] = highlight_mask
collection_maple.metadata
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SET-ID</th>
      <th>STEP</th>
      <th>LK</th>
      <th>highlight</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>0</td>
      <td>0.983439</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>1</td>
      <td>0.983822</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>2</td>
      <td>0.984404</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>3</td>
      <td>0.985124</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>IQtreeStartingTree_Trees</td>
      <td>4</td>
      <td>0.985631</td>
      <td>0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>133</th>
      <td>UshERStartingTree_Trees</td>
      <td>26</td>
      <td>0.999844</td>
      <td>0</td>
    </tr>
    <tr>
      <th>134</th>
      <td>UshERStartingTree_Trees</td>
      <td>27</td>
      <td>0.999850</td>
      <td>0</td>
    </tr>
    <tr>
      <th>135</th>
      <td>UshERStartingTree_Trees</td>
      <td>28</td>
      <td>0.999976</td>
      <td>0</td>
    </tr>
    <tr>
      <th>136</th>
      <td>UshERStartingTree_Trees</td>
      <td>29</td>
      <td>1.000000</td>
      <td>1</td>
    </tr>
    <tr>
      <th>137</th>
      <td>TrueTreeSimulations</td>
      <td>0</td>
      <td>0.706890</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
<p>138 rows × 4 columns</p>
</div>




```python
collection_maple.plot_2D('pcoa', plot_meta = 'LK', same_scale = True, select = False)
```
