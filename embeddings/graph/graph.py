__author__ = 'Andrea Rubbi'

import plotly.graph_objects as go
import pandas as pd
import numpy as np
from ipywidgets import widgets
import builtins
import matplotlib.colors as mcolors
from pylab import cm
import matplotlib
import random
import warnings


def plot_embedding(data, metadata, dimensions, save=False, name_plot='Tree_embedding', static=False, plot_meta = 'SET-ID', plot_set = None, select = False):
    """Plot embedding of distance matrix - in 2D or 3D

    Args:
        data (pandas.DataFrame): embedding of distance matrix 
        metadata (pandas.DataFrame): metadata of tree_set or set_collection
        dimensions (int): number of dimensions - either 2 or 3
        save (bool, optional): save plot in pdf format. Defaults to False.
        name_plot (str, optional): name of plot. Defaults to 'Tree_embedding'.
        static (bool, optional): if True, returns a less interactive format of plot. Defaults to False.
        plot_meta (str, optional): defines the meta-feature used to color the points. Defaults to 'SET-ID'.
        plot_set (list or str, optional): sets in set_collection to be plotted. Defaults to None.
        select (bool, optional): if True, generates widgets that allow to show or hide uo to 16 set traces. Defaults to False.

    Returns:
       image: plot with related widgets - interactive or static format
    """
    
    # ─── Format Controls ──────────────────────────────────────────────────
    assert dimensions <= 3 and dimensions > 1, 'Please select either 2 or 3 dimensions for the plot'  
    assert plot_meta in metadata.columns, f'Could not find {plot_meta} in metadata'
    assert 'SET-ID' in metadata.columns, 'Could not find SET-ID in metadata - please redefine SET-ID or reset the tree_set/set_collection and subsequently define additional metadata'
    
    # ─── plot_set Selects The Tree Set To Be Displayed ────────────────────
    if type(plot_set) == type(None): plot_set = np.unique(metadata['SET-ID'])
    else: plot_set = list(plot_set)

    idxs_plot_set = [Tset in plot_set for Tset in metadata['SET-ID']]
    data, metadata = data[idxs_plot_set], metadata[idxs_plot_set]
    
    # ─── Define Color Map And Mapping Array For Each Metadata ─────────────
    metadata_colors = dict()
    for meta in metadata.columns: # for each metadata column
        elements = np.unique(metadata[meta]) # unique elements
        color_plot = np.array(metadata[meta].values) # get meta column
        
        cmap = cm.get_cmap('jet', max(20,len(elements)))
        # here we generate a color map --> if there are more than 1000 unique elements,
        # we assume that the variable is some sort of common parameter. e.g. likelihood, step, parsimony...
        if len(elements) > 100: col_list = "jet" # jet color scale
        # else, we presume this is some kind of difference between sets - e.g. different tree_set
        # we generate random colors to maximize the visual division
        else: col_list = [mcolors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]
        
        # for each single element in meta, we generate a common color_plot array
        # where each element is encoded as distinctive number
        for i, elem in enumerate(elements):
            idx = color_plot == elem
            color_plot[idx] = i
        
        metadata_colors[f'{meta}_color_plot'] = list()
        metadata_colors[f'{meta}_color_map'] = list()
        
        Sets = np.unique(metadata['SET-ID'])
        for i, SetID in enumerate(Sets):
            idx = metadata['SET-ID'] == SetID
            color_plot_set = color_plot[idx]
            color = col_list
            if len(np.unique(color_plot_set)) == 1: color = random.sample([mcolors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)], 10)
            metadata_colors[f'{meta}_color_plot'].append(color_plot_set.tolist())
            metadata_colors[f'{meta}_color_map'].append(color)
            
        
    # initialize go.Figure()
    fig = go.Figure()
    
    # widget for the modification of metadata coloring
    meta_widget = widgets.Dropdown(
        options=list(metadata.columns),
        value=plot_meta,
        description='Metadata:',
        )
    
    # widget to save plot current view as pdf
    save_pdf = widgets.Button(
            description='Save plot as PDF',
            button_style='danger',
            layout= widgets.Layout(width="150px"),)
    
    # define response_meta to change color of points with different metadata selections
    def response_meta(change):
        traces_start_at = 0
        if dimensions == 2: traces_start_at = 1
        with fig.batch_update():
            for i in range(nUnique):
                fig.data[i + traces_start_at].marker['color'] = metadata_colors[f'{meta_widget.value}_color_plot'][i]
                fig.data[i + traces_start_at].marker['colorscale'] = metadata_colors[f'{meta_widget.value}_color_map'][i]

    # bind widget to function response_metadata
    meta_widget.observe(response_meta, names="value")
    
    # define save_pdf_funct to save pdf
    def save_pdf_func(b):
        fig.write_image(name_plot + '.pdf')
    
    # bind widget to function save_pdf_func
    save_pdf.on_click(save_pdf_func)
    
    
    # ─── 3D Plot ──────────────────────────────────────────────────────────
    if dimensions == 3: # if input dimension is 3
        assert data.shape[1] == 3, 'Embed distance_matrix in 3D before requesting a 3D plot'
        
        # add a scatter3d trace for each tree_set - i.e. SET-ID unique value
        Sets, nSetID = np.unique(metadata['SET-ID'], return_counts=True)
        # number of sets
        nUnique = len(Sets)
        for i, SetID in enumerate(Sets):
            idx = metadata['SET-ID'] == SetID
            fig.add_trace(go.Scatter3d(x=data[idx,0],
                                y=data[idx,1],
                                z=data[idx,2], mode= 'markers', 
                                showlegend=False,
                                visible=True,
                                marker_color=metadata_colors[f'{plot_meta}_color_plot'][i],
                                text=metadata['SET-ID'].values[idx],
                                opacity = 0.7,
                                marker=dict(
                                colorscale=metadata_colors[f'{plot_meta}_color_map'][i],
                                size = 10,
                                line=dict(
                                        #color='MediumPurple',
                                        width=10,
                                        color = "rgba( 30, 30, 30, 0.3)",
                                    ),
                                ),
                                ))
        
        # defines the layout of the plot    
        fig.update_layout(
            width=900,
            height=900 + 10*nUnique,
            autosize=True,
            margin=dict(t=0, b=100, l=0, r=0),
            template="seaborn",
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            #hovermode="x unified",
        )

        # updates 3D scene options
        fig.update_scenes(
            aspectratio=dict(x=1, y=1, z=1),
            aspectmode="manual"
        )
        
        
        buttons_ID = list()
        # buttons_ID are switch on-off buttons that allow to show or hyde specific tree_sets | activated by select
        # TODO : transform double buttons in single toggle button --> specify the hyde-button behaviours in arg2
        if select:
            if nUnique > 16:
                warnings.warn('Selection not available for number of sets > 16 - please select a smaller range of sets using plot_set')
            else:
                for i,ID in enumerate(Sets):
                    if len(ID) > 15: ID = ID[:15] + '...'
                        
                    buttons_ID.append(dict(type="buttons",
                                        direction="left",
                                        buttons=list([
                                                dict(
                                                    args=[dict(visible=True), [i]],
                                                    label=f"Show {ID}",
                                                    method="restyle"
                                                    ),
                                                dict(
                                                    args=[dict(visible=False ), [i]],
                                                    label=f"Hide {ID}",
                                                    method="restyle"
                                                    ),
                                            ]),
                                            pad={"r": 0, "t": 0},
                                            showactive=True,
                                            x=0 + 0.4 * int(round((i+1)/(nUnique))),
                                            xanchor="left",
                                            y=1.5 - 0.06*(i - nUnique//2 * int(round((i+1)/(nUnique)))),
                                            yanchor="top"),
                                    )
        
        
        
        # buttons_meta define the metadata variable used to color the traces
        # TODO : this is currently not working for some reson - coming back to it asap
        buttons_meta = [dict(args=[{"marker" : {"color" : f"metadata_colors['{plot_meta}_color_plot']",
                                    #"colorscale" : f"metadata_colors['{plot_meta}_color_map']"
                                    }},
                                    ],
                                    label=f"Meta {plot_meta}",
                                    method="update"
                                    )]
        
        for i, meta in enumerate(metadata.columns):
            if meta == plot_meta: continue
            buttons_meta.append(dict(args=[{"marker" : {"color" : f"metadata_colors['{meta}_color_plot']",
                                    #"colorscale" : f"metadata_colors['{meta}_color_map']"
                                    }},
                                    ],
                                    label=f"Meta {meta}",
                                    method="update"
                                    ))
            
        buttons_meta = [dict(type="dropdown", 
                                direction="down",
                                buttons=buttons_meta,
                                pad={"r": 0, "t": 0},
                                showactive=True,
                                x=0.5,
                                xanchor="left",
                                y=1.074,
                                yanchor="top")]
                              

        # add buttons to figure
        # updatemenus is a list of distinct dictionaries defining 
        # widgets controlling different parameters in the figure 
        fig.update_layout(
            updatemenus=[
                # Markers or Markers + Lines
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        dict(
                            args=[{"mode":['markers' for i in range(nUnique)]}],
                            label="Markers",
                            method="restyle"
                        ),
                        dict(
                            args=[{"mode":['lines+markers' for i in range(nUnique)]}],
                            label="Markers & Lines",
                            method="restyle"
                        ),
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.1,
                    xanchor="left",
                    y= 1.09,
                    yanchor="top"
                ),
            ] + buttons_ID #+ buttons_meta

            )
            
        # add label
        fig.update_layout(
            annotations=[
                dict(text="<b>Points:", showarrow=False,
                                    x=0, y=1.06, yref="paper", xref="paper", align="left"),
            ]
        )   
        
        
        # define container for widgets
        container = widgets.HBox([meta_widget, save_pdf]) 
        
        # "static" plot
        no_widget_fig = fig
        
        # "non-static" plot 
        fig = go.FigureWidget(fig)
        
    
    # ─── 2d Plot ──────────────────────────────────────────────────────────
    else: # dimension can be either 2 or 3 - if dim == 3: ... else:
        
        # Getting the margins for the plot axes 
        max_x, min_x = max(data[:,0]), min(data[:,0])
        max_y, min_y = max(data[:,1]), min(data[:,1])
        
        # add scatter trace of whole collection
        fig.add_trace(go.Scatter(x=data[:,0],
                                   y=data[:,1],
                                   mode= 'markers',
                                   xaxis='x',
                                   yaxis='y',
                                   visible=False,
                                   showlegend=False, 
                                   marker_color=color_plot,
                                   opacity=0.7,
                                   text=metadata[plot_meta],
                                   marker=dict(
                                   colorscale=col_list,)))
        
        # for each tree_set add a scatter trace to the figure
        Sets, nSetID = np.unique(metadata['SET-ID'], return_counts=True)
        # number of unique tree_sets
        nUnique = len(Sets)
        for i, SetID in enumerate(Sets):
            idx = metadata['SET-ID'] == SetID
            fig.add_trace(go.Scatter(x=data[idx,0],
                                y=data[idx,1],
                                mode= 'markers', 
                                showlegend=False,
                                visible=True,
                                marker_color=metadata_colors[f'{plot_meta}_color_plot'][i],
                                text=metadata['SET-ID'].values[idx],
                                opacity = 0.7,
                                marker=dict(
                                    colorscale=metadata_colors[f'{plot_meta}_color_map'][i],
                                    size = 10,    
                                    line=dict(
                                        #color='MediumPurple',
                                        width=2,
                                    ),
                                ),
                                ))
        
        # ─── Marginal Distribution Of Points ──────────────────────────
        fig.add_trace(go.Histogram(
                y = data[:,1],
                xaxis = 'x2',
                showlegend=False,
                visible=True,
                marker = dict(
                    color = 'rgba(200,200,250,1)'
                )
            ))
        
        fig.add_trace(go.Histogram(
                x=data[:,0],
                yaxis = 'y2',
                showlegend=False,
                visible=True,
                marker = dict(
                    color = 'rgba(200,200,250,1)'
                )
            ))
        
        # ─────────────────────────────────────────────────────────────
            
        # defines the layout of the plot  
        fig.update_layout(
            width=800,
            height=1200,
            autosize=False,
            margin=dict(t=0, b=0, l=0, r=0),
            template="seaborn",
            xaxis= dict(showgrid=True,
                        showline=True,
                        zeroline=True,
                        showticklabels = True,
                        domain = [0,0.85],
                        range = [min_x - abs(0.1 * min_x), max_x + abs(0.1 * max_x)]
                        ),
            yaxis= dict(showgrid=True,
                        showline=True,
                        zeroline=True,
                        showticklabels = True,
                        domain = [0,0.85],
                        range = [min_y - abs(0.1 * min_y), max_y + abs(0.1 * max_y)]
                        ),
            plot_bgcolor= 'rgb(240,240,240)',
            xaxis2 = dict(
                zeroline = True,
                showticklabels = True,
                domain = [0.85,1],
                showgrid = True
            ),
            yaxis2 = dict(
                zeroline = True,
                showticklabels = True,
                domain = [0.85,1],
                showgrid = True
            ),
            paper_bgcolor='rgba(0,0,0,0)',
            #hovermode="x unified",
        )
        
        
        buttons_ID = list()
        # buttons_ID are switch on-off buttons that allow to show or hyde specific tree_sets | activated by select
        # TODO : transform double buttons in single toggle button --> specify the hyde-button behaviours in arg2
        if select:
            if nUnique > 16:
                warnings.warn('Selection not available for number of sets > 16 - please select a smaller range of sets using plot_set')
            else:
                for i,ID in enumerate(Sets):
                    if len(ID) > 15: ID = ID[:15] + '...'
                        
                    buttons_ID.append(dict(type="buttons",
                                        direction="left",
                                        buttons=list([
                                                dict(
                                                    args=[dict(visible=True), [i + 1]],
                                                    label=f"Show {ID}",
                                                    method="restyle"
                                                    ),
                                                dict(
                                                    args=[dict(visible=False ), [i + 1]],
                                                    label=f"Hide {ID}",
                                                    method="restyle"
                                                    ),
                                            ]),
                                            pad={"r": 0, "t": 0},
                                            showactive=True,
                                            x=0 + 0.48 * int(round((i+1)/(nUnique))),
                                            xanchor="left",
                                            y=1.8 - 0.06*(i - nUnique//2 * int(round((i+1)/(nUnique)))),
                                            yanchor="top"),
                                    )
        
        # add updatemenus
        fig.update_layout(
            updatemenus=[
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        # change graph from scatterplot to histogram2dconotour plot
                        dict(
                            args=[{"type":  ["scatter"] + ["scatter" for i in range(nUnique)] + ["histogram", "histogram"],
                                   "visible": [False] + [True for i in range(nUnique+2)],},
                                  {"xaxis":{"showgrid":True,
                                                "showline":True,
                                                "zeroline":True,
                                                "showticklabels":True,
                                                "domain":[0,0.85],},
                                   "yaxis":{"showgrid":True,
                                                "showline":True,
                                                "zeroline":True,
                                                "showticklabels":True,
                                                "domain":[0,0.85],},
                                   # xaxis 2 and yaxis2 refer to the second set of axes 
                                   # that are needed to plot the marginal distribution
                                   # the range indicates the portion of total graph space
                                   # allocated to that specific axis
                                   "xaxis2":{"domain":[0.85,1],
                                             "showgrid":True,
                                             "showline":True,
                                             "zeroline":True,
                                             "showticklabels":True,},
                                   "yaxis2":{"domain":[0.85,1],
                                             "showgrid":True,
                                             "showline":True,
                                             "zeroline":True,
                                             "showticklabels":True,},
                                   
                                   "plot_bgcolor":'rgb(240,240,240)',
                                   },
                                  {"skip_invalid" : True,},
                                  ],
                            label="Scatter2D",
                            method="update"
                        ),
                        dict(
                            args=[{"type": ["histogram2dcontour"] + ["scatter" for i in range(nUnique)] + ["histogram", "histogram"], 
                                   "visible": [True] + [False for i in range(nUnique+2)],}, #"showscale": True, "contours":{"coloring":"fill"}
                                  {"xaxis":{"showgrid":True,
                                                "showline":True,
                                                "zeroline":True,
                                                "showticklabels":True,
                                                "domain":[0,1],
                                                },
                                   "yaxis":{"showgrid":True,
                                                "showline":True,
                                                "zeroline":True,
                                                "showticklabels":True,
                                                "domain":[0,1],
                                                },
                                   "xaxis2":{"domain":[0.999,1],
                                             "showgrid":False,
                                                "showline":False,
                                                "zeroline":False,
                                                "showticklabels":False,},
                                   "yaxis2":{"domain":[0.999,1],
                                             "showgrid":False,
                                             "showline":False,
                                             "zeroline":False,
                                             "showticklabels":False,},
                                   "plot_bgcolor":'rgb(240,240,240)',
                                   },
                                  {"skip_invalid" : True,}
                                  ],
                            label="Contour",
                            method="update"
                        )
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.15,
                    xanchor="left",
                    y=1.3,
                    yanchor="top",
                ),
                # Markers or Markers + Lines
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        dict(
                            args=[{"mode":"markers"},list(range(1, nUnique+1))],
                            label="Markers",
                            method="restyle"
                        ),
                        dict(
                            args=[{"mode":"lines+markers"},list(range(1, nUnique+1))], 
                            label="Markers & Lines",
                            method="restyle"
                        ),
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.15,
                    xanchor="left",
                    y=1.2,
                    yanchor="top"
                ),
            ] + buttons_ID
        )

        # Add annotation
        fig.update_layout(
            annotations=[
                dict(text="<b>Plot:", showarrow=False,
                                    x=0, y=1.275, yref="paper", xref="paper", align="left"),
                dict(text="<b>Points:", showarrow=False,
                                    x=0, y=1.175, yref="paper", xref="paper", align="left"),
            ]
        )   
        
        # shows points on contour plot
        def response_points(b):
            if fig.data[0].type == "histogram2dcontour":
                if not np.array([fig.data[i].visible for i in range(1, nUnique+1)]).all():
                    for i in range(1, nUnique+1): fig.data[i].visible = True
                else:
                    for i in range(1, nUnique+1): fig.data[i].visible = False

        # widget button -> response_points 
        points_density = widgets.Button(
            description='Show points on Contour plot',
            button_style='info',
            layout= widgets.Layout(width="200px"),)    
        
        # observer meta values (dropdown menu)
        meta_widget.observe(response_meta, names="value")
        # oberver points_density button
        points_density.on_click(response_points)
        
        container = widgets.HBox([meta_widget,points_density, save_pdf]) 
        
        no_widget_fig = fig
        
        fig = go.FigureWidget(fig)
    
    
    # ─────────────────────────────────────────────────────────────────────
    
    # update_point hihglights specific traces on click
    def update_point(trace, points, selector):  
        if points.point_inds:
            traces_start_at = 0
            if dimensions == 2:
                if points.trace_index == 0: return
                traces_start_at = 1
                
            idx_trace = points.trace_index - traces_start_at
            other_idxs = list(range(traces_start_at, nUnique + traces_start_at))
            other_idxs.pop(idx_trace)

            if trace.opacity == 1.0: 
                trace.opacity = 0.7
                for idx in other_idxs: fig.data[idx].opacity = 0.7
            elif trace.opacity == 0.05: trace.opacity = 1.0
            else: 
                trace.opacity = 1.0
                for idx in other_idxs: fig.data[idx].opacity = 0.05
    
    # binds update_point function on every trace
    for trace in fig.data:
            trace.on_click(update_point)
        
    # VBox contains both widgets and the figure
    image = widgets.VBox([container, fig])   
    
    # if save is requested, then a pdf file is returned 
    if save: no_widget_fig.write_html(name_plot + '.html')
    
    # if static is requested, then a less interactive plot is returned
    #! NB: this can be the only way to obtain a plot in cases where 
    #! plotly - ipywidgets - jupyter versions are incompatible
    if static: return no_widget_fig
    # check if we are working in an interactive environment
    # hence ok using widgets
    if hasattr(builtins,'__IPYTHON__'): return image
    else: return no_widget_fig
