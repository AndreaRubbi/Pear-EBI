import plotly.graph_objects as go
import pandas as pd
import numpy as np
from ipywidgets import widgets
import builtins


def plot_embedding(data, metadata, dimensions, save=False, name_plot='Tree_embedding', static=False, plot_meta = 'SET-ID'):
    assert dimensions <= 3 and dimensions > 1, 'Please select either 2 or 3 dimensions for the plot'  
    assert plot_meta in metadata.columns, f'Could not find {plot_meta} in metadata'
    assert 'SET-ID' in metadata.columns, 'Could not find SET-ID in metadata - please redefine SET-ID or reset the tree_set/set_collection and subsequently define additional metadata'
    
    sets = np.unique(metadata[plot_meta])
    color_plot = np.array(metadata[plot_meta].values)
    col_list = ["red", "green", "cyan", "blue", "yellow", 'black', "magenta"]
    for i, elem in enumerate(sets):
        idx = metadata[plot_meta] == elem
        color_plot[idx] = i
    
    if len(col_list) < len(sets): col_list = 'Jet'
    
    fig = go.Figure()
    
    meta_widget = widgets.Dropdown(
        options=list(metadata.columns),
        value=plot_meta,
        description='Metadata:',
        )
    save_pdf = widgets.Button(
            description='Save plot as PDF',
            button_style='danger',
            layout= widgets.Layout(width="150px"),)
    
    def save_pdf_func(b):
        fig.write_image(name_plot + '.pdf')
    
    save_pdf.on_click(save_pdf_func)
    
    if dimensions == 3:
        assert data.shape[1] == 3, 'Embed distance_matrix in 3D before requesting a 3D plot'
        fig.add_trace(go.Scatter3d(x=data[:,0],
                                   y=data[:,1],
                                   z=data[:,2], mode= 'markers', 
                                   showlegend=False,
                                   marker_color=color_plot,
                                   text=metadata[plot_meta],
                                   opacity = 0.6,
                                   marker=dict(
                                   colorscale=col_list,
                                   size = 10)))
        
        SetID, nSetID = np.unique(metadata['SET-ID'], return_counts=True)
        nUnique = len(SetID)
        for i, uSet in enumerate(SetID):
            if i == 0: start = 0
            else: start = sum(nSetID[:i])
            end = start + int(nSetID[i])
            if end - start == 1: continue
            fig.add_trace(go.Scatter3d(x=data[start:end,0],
                                   y=data[start:end,1],
                                   z=data[start:end,2], mode= 'lines', 
                                   showlegend=False,
                                   visible=False,
                                   marker_color=color_plot[start:end],
                                   text=metadata[plot_meta].values[start:end],
                                   opacity = 0.6,
                                   marker=dict(
                                   colorscale=col_list,
                                   size = 10)))
        
        
            
        fig.update_layout(
            width=800,
            height=900,
            autosize=False,
            margin=dict(t=0, b=0, l=0, r=0),
            template="seaborn",
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            #hovermode="x unified",
        )

        # Update 3D scene options
        fig.update_scenes(
            aspectratio=dict(x=1, y=1, z=1),
            aspectmode="manual"
        )

        # Add dropdown
        fig.update_layout(
            updatemenus=[
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        dict(
                            args=[{"visible":False},list(range(1,1+nUnique))],
                            label="Markers",
                            method="restyle"
                        ),
                        dict(
                            args=[{"visible":True},list(range(1,1+nUnique))],
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
            ]
        )

        # Add annotation
        fig.update_layout(
            annotations=[
                dict(text="<b>Points:", showarrow=False,
                                    x=0, y=1.175, yref="paper", xref="paper", align="left"),
            ]
        )   
        
        def response_meta(change):
            sets = np.unique(metadata[meta_widget.value])
            color_plot = np.array(metadata[meta_widget.value].values)
            col_list = ["red", "green", "cyan", "blue", "yellow", 'black', "magenta"]
            for i, elem in enumerate(sets):
                idx = metadata[meta_widget.value] == elem
                color_plot[idx] = i
            
            if len(col_list) < len(sets): col_list = 'Jet'

            with fig.batch_update():
                fig.data[0].marker['color'] = color_plot
                fig.data[0].marker['colorscale'] = col_list

            
        
        meta_widget.observe(response_meta, names="value")
        
        container = widgets.HBox([meta_widget, save_pdf]) 
        
        no_widget_fig = fig
        
        fig = go.FigureWidget(fig)
        
        image = widgets.VBox([container, fig])  
    
    else:
        fig.add_trace(go.Scatter(x=data[:,0],
                                   y=data[:,1],
                                   mode= 'markers', 
                                   showlegend=False,
                                   marker_color=color_plot,
                                   text=metadata[plot_meta],
                                   opacity = 0.6,
                                   marker=dict(
                                   colorscale=col_list,
                                   size = 10)))
        
        fig.add_trace(go.Scatter(x=data[:,0],
                                   y=data[:,1],
                                   mode= 'markers',
                                   xaxis='x',
                                   yaxis='y',
                                   visible=False,
                                   showlegend=False, 
                                   marker_color=color_plot,
                                   text=metadata[plot_meta],
                                   marker=dict(
                                   colorscale=col_list,)))
        
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
        
        SetID, nSetID = np.unique(metadata['SET-ID'], return_counts=True)
        nUnique = len(SetID)
        for i, uSet in enumerate(SetID):
            if i == 0: start = 0
            else: start = sum(nSetID[:i])
            end = start + int(nSetID[i])
            if end - start == 1: continue
            fig.add_trace(go.Scatter(x=data[start:end,0],
                                   y=data[start:end,1],
                                   mode= 'lines', 
                                   showlegend=False,
                                   visible=False,
                                   marker_color=color_plot[start:end],
                                   text=metadata[plot_meta].values[start:end],
                                   opacity = 0.6,
                                   marker=dict(
                                   colorscale=col_list,
                                   size = 10)))
            
            
        fig.update_layout(
            width=800,
            height=900,
            autosize=False,
            margin=dict(t=0, b=0, l=0, r=0),
            template="seaborn",
            xaxis= dict(showgrid=True,
                        showline=True,
                        zeroline=True,
                        showticklabels = True,
                        domain = [0,0.85],
                        ),
            yaxis= dict(showgrid=True,
                        showline=True,
                        zeroline=True,
                        showticklabels = True,
                        domain = [0,0.85],
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

        # Add dropdown
        fig.update_layout(
            updatemenus=[
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        dict(
                            args=[{"type": ["scatter", "scatter", "histogram", "histogram"] + ["scatter" for i in range(nUnique)],
                                   "visible":[True, False, True, True] + [False for i in range(nUnique)]},
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
                                  {"skip_invalid" : [True, True, True, True],},
                                  ],
                            label="Scatter2D",
                            method="update"
                        ),
                        dict(
                            args=[{"type": ["histogram2dcontour", "scatter",  "histogram", "histogram"] + ["scatter" for i in range(nUnique)],
                                   "visible":[True, False, False, False] + [False for i in range(nUnique)]},
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
                                  {"skip_invalid" : [True, True, True, True],}
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
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        dict(
                            args=[{"visible":False},list(range(4,4+nUnique))],
                            label="Markers",
                            method="restyle"
                        ),
                        dict(
                            args=[{"visible":True},list(range(4,4+nUnique))], 
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
            ]
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
        
        def response_meta(change):
            sets = np.unique(metadata[meta_widget.value])
            color_plot = np.array(metadata[meta_widget.value].values)
            col_list = ["red", "green", "cyan", "blue", "yellow", 'black', "magenta"]
            for i, elem in enumerate(sets):
                idx = metadata[meta_widget.value] == elem
                color_plot[idx] = i
            
            if len(col_list) < len(sets): col_list = 'Jet'

            with fig.batch_update():
                fig.data[0].marker['color'] = color_plot
                fig.data[0].marker['colorscale'] = col_list
                
        def response_points(b):
            if fig.data[0].type == "histogram2dcontour" and not fig.data[1].visible:
                fig.data[1].visible = True
            else:fig.data[1].visible = False

        points_density = widgets.Button(
            description='Show points on Contour plot',
            button_style='info',
            layout= widgets.Layout(width="200px"),)    
        
        meta_widget.observe(response_meta, names="value")
        points_density.on_click(response_points)
        
        container = widgets.HBox([meta_widget,points_density, save_pdf]) 
        
        no_widget_fig = fig
        
        fig = go.FigureWidget(fig)
        
        image = widgets.VBox([container, fig])     
    
    if save: no_widget_fig.write_html(name_plot + '.html')
    
    if static: return no_widget_fig
    if hasattr(builtins,'__IPYTHON__'): return image
    else: return no_widget_fig
