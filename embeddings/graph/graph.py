import plotly.graph_objects as go
import pandas as pd
import numpy as np
from ipywidgets import widgets
from contextlib import suppress


def plot_embedding(data, metadata, dimensions, name_plot='Tree_embedding'):
    assert dimensions <= 3 and dimensions > 1, 'Please select either 2 or 3 dimensions for the plot'  
    
    fig = go.Figure()
    
    meta_widget = widgets.Dropdown(
        options=list(metadata.columns),
        value='SET-ID',
        description='Metadata:',
        )
    
    points_density = widgets.Button(
    description='Show points on Contour plot',
    button_style='info',
    layout= widgets.Layout(width="200px"),)
    
    with suppress(AttributeError):
        data_reformatted = {"D1":data[:,0], "D2":data[:,1], "D3":data[2]}
        sets = np.unique(metadata['SET-ID'])
        color_plot = np.array(metadata['SET-ID'].values)
        col_list = ["red", "green", "cyan", "blue", "yellow", 'black', "magenta"]
        for i, elem in enumerate(sets):
            idx = metadata['SET-ID'] == elem
            color_plot[idx] = i
        
        if len(col_list) < len(sets): col_list = 'Jet'
        
        
        fig.add_trace(go.Scatter3d(x=data[:,0],
                                   y=data[:,1],
                                   z=data[:,2], mode= 'markers', 
                                   showlegend=False,
                                   marker_color=color_plot,
                                   text=metadata['SET-ID'],
                                   marker=dict(
                                   colorscale=col_list,)))
        
        fig.add_trace(go.Scatter(x=data[:,0],
                                   y=data[:,1],
                                   mode= 'markers',
                                   xaxis='x',
                                   yaxis='y',
                                   visible=False,
                                   showlegend=False, 
                                   marker_color=color_plot,
                                   text=metadata['SET-ID'],
                                   marker=dict(
                                   colorscale=col_list,)))
        
        fig.add_trace(go.Histogram(
                y = data[:,1],
                xaxis = 'x2',
                showlegend=False,
                visible=False,
                marker = dict(
                    color = 'rgba(200,200,250,1)'
                )
            ))
        
        fig.add_trace(go.Histogram(
                x=data[:,0],
                yaxis = 'y2',
                showlegend=False,
                visible=False,
                marker = dict(
                    color = 'rgba(200,200,250,1)'
                )
            ))
            
            
        fig.update_layout(
            width=800,
            height=900,
            autosize=False,
            margin=dict(t=0, b=0, l=0, r=0),
            template="seaborn",
            xaxis= dict(showgrid=False,
                        showline=False,
                        zeroline=False,
                        showticklabels = False,
                        domain = [0,1],
                        ),
            yaxis= dict(showgrid=False,
                        showline=False,
                        zeroline=False,
                        showticklabels = False,
                        domain = [0,1],
                        ),
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis2 = dict(
                zeroline = False,
                showticklabels = False,
                domain = [0,0],
                showgrid = False
            ),
            yaxis2 = dict(
                zeroline = False,
                showticklabels = False,
                domain = [0,0],
                showgrid = False
            ),
            #paper_bgcolor='rgba(0,0,0,0)',
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
                            args=[{"type":["scatter3d", "scatter", "histogram", "histogram"],
                                   "visible":[True, False, False, False]},
                                  { "scene":{"aspectratio":dict(x=1,y=1,z=1),
                                             "aspectmode":"manual",
                                             "xaxis":{"showgrid":True,
                                                      "showline":True,
                                                      "zeroline":True,
                                                      "showticklabels":True,
                                                      "domain":[0,1],},
                                             "yaxis":{"showgrid":True,
                                                      "showline":True,
                                                      "zeroline":True,
                                                      "showticklabels":True,
                                                      "domain":[0,1],},
                                             "zaxis":{"showgrid":True,
                                                      "showline":True,
                                                      "zeroline":True,
                                                      "showticklabels":True,
                                                      "domain":[0,1],},},
                                      "xaxis":{"showgrid":False,
                                                "showline":False,
                                                "zeroline":False,
                                                "showticklabels":False,
                                                "domain":[0,1],},
                                   "yaxis":{"showgrid":False,
                                                "showline":False,
                                                "zeroline":False,
                                                "showticklabels":False,
                                                "domain":[0,1],},
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
                                   "plot_bgcolor":'rgba(0,0,0,0)',
                                   },
                                  {"skip_invalid" : [True, True, True, True],},
                                  ],
                            label="Scatter3D",
                            method="update",
                            
                        ),
                        dict(
                            args=[{"type": ["scatter", "scatter", "histogram", "histogram"],
                                   "visible":[True, False, True, True]},
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
                            args=[{"type": ["histogram2dcontour", "scatter",  "histogram", "histogram"],
                                   "visible":[True, False, False, False]},
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
                            args=[{"mode":"markers"},[0,1]],
                            label="Markers",
                            method="restyle"
                        ),
                        dict(
                            args=[{"mode":"lines+markers"},[0,1]], #[0, 1] indicates the traces that have to be modified!
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

            
        
        meta_widget.observe(response_meta, names="value")
        points_density.on_click(response_points)
        
        container = widgets.HBox([meta_widget,points_density]) 
        
        fig = go.FigureWidget(fig)
        
        image = widgets.VBox([container, fig])     
    
    return image
