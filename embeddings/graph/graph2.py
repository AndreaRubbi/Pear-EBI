# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, html, dcc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from ipywidgets import widgets
import base64

app = Dash(__name__)

df = pd.DataFrame({
    "Fruit": ["Apples", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    "Amount": [4, 1, 2, 2, 4, 5],
    "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
})

fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")
import os
print()
image_filename = os.path.join(os.getcwd(), 'TreeEmbedding/goldman_logo.png')# replace with your own image
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

app.layout = html.Div(children=[
    html.H1(children='TreeEmbedding'),
    html.H4(children='EMBL-EBI Goldman Group'),
    
    html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
             style={"width":"150px",
                    "position":"absolute",
                    "top":"8px",
                    "left":"350px",}),
    
    html.Div(children='''
        Inreactive plotting interface
    '''),

    dcc.Graph(
        id='example-graph',
        figure=fig,
        style={"position":"relative",
            "top":"100px",}
    ),
    ])


if __name__ == "__main__":
    app.run_server(debug=True)

def plot_embedding(data, metadata, dimensions, name_plot='Tree_embedding'):
    assert dimensions <= 3 and dimensions > 1, 'Please select either 2 or 3 dimensions for the plot'  
    
    fig = go.FigureWidget()
    
    meta_widget = widgets.Dropdown(
        options=list(metadata.columns),
        value='SET-ID',
        description='<b> Metadata:',
        )
    
    if dimensions == 3:
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
                                   marker_color=color_plot,
                                   text=metadata['SET-ID'],
                                   marker=dict(
                                   colorscale=col_list,)))
        fig.update_layout(
            width=800,
            height=900,
            autosize=False,
            margin=dict(t=0, b=0, l=0, r=0),
            template="seaborn",
        )

        # Update 3D scene options
        fig.update_scenes(
            aspectratio=dict(x=1, y=1, z=1.0),
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
                            args=["type", "scatter3d"],
                            label="Scatter3D",
                            method="restyle"
                        ),
                        dict(
                            args=["type", "scatter"],
                            label="Scatter2D",
                            method="restyle"
                        ),
                        dict(
                            args=[{"type": "histogram2dcontour",
                                   "autobinx" : True, 
                                   "autobiny" : True}],
                            label="Contour",
                            method="update"
                        )
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.15,
                    xanchor="left",
                    y=1.3,
                    yanchor="top"
                ),
                dict(
                    type = "buttons",
                    direction = "left",
                    buttons=list([
                        dict(
                            args=["mode", "markers"],
                            label="Markers",
                            method="restyle"
                        ),
                        dict(
                            args=["mode", "lines+markers"],
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
        
        
    else: pass
    
    def response(change):
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
    
    meta_widget.observe(response, names="value")
    
    container = widgets.HBox([meta_widget])
    
    image = widgets.VBox([container, fig])
    
    return image