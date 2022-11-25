import plotly.graph_objects as go
import pandas as pd



def plot_embedding(data, metadata, dimensions, name_plot='Tree_embedding'):
    assert dimensions <= 3 and dimensions > 1, 'Please select either 2 or 3 dimensions for the plot'  
    
    fig = go.Figure()
    
    if dimensions == 3:
        fig.add_trace(go.Scatter3d(data, x=0, y=1, z=2))

        # Update plot sizing
        fig.update_layout(
            width=800,
            height=900,
            autosize=False,
            margin=dict(t=0, b=0, l=0, r=0),
            template="plotly_white",
        )

        # Update 3D scene options
        fig.update_scenes(
            aspectratio=dict(x=1, y=1, z=0.7),
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
                            args=["type", "Scatter3D"],
                            label="3D",
                            method="restyle"
                        ),
                        dict(
                            args=["type", "Scatter"],
                            label="2D",
                            method="restyle"
                        )
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.11,
                    xanchor="left",
                    y=1.1,
                    yanchor="top"
                ),
            ]
        )

        # Add annotation
        fig.update_layout(
            annotations=[
                dict(text="Trace type:", showarrow=False,
                                    x=0, y=1.08, yref="paper", align="left")
            ]
        )

        fig.write_html(name_plot + '.html')
        
        
        
        
    else: pass
        
        
        
        
            
        fig = px.scatter_3d(
            components, x=0, y=1, z=2,
            title=f'Total Explained Variance: {total_var:.2f}%',
            labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
        )
        

    fig = px.scatter(components, x=0, y=1, 
            title=f'Total Explained Variance: {total_var:.2f}%') #color=Distances[WHAT]
    fig.write_html("./graph_PCA2.html")

if n_dimensions >= 3:    
        fig = px.scatter_3d(Distances_embedded_ND, x=0, y=1, z=2,
                 labels={'0': 'D1', '1': 'D2', '2': 'D3'},
                 title=f't-SNE 3D Embedding') #color=Distances[WHAT]
        fig.write_html("./graph_t-SNE3.html")
          
    fig = px.scatter(Distances_embedded_ND, x=0, y=1, 
                    title=f't-SNE 2D Embedding') #color=Distances[WHAT]
    fig.write_html("./graph_t-SNE2.html")
    





# create figure


# Add surface trace
