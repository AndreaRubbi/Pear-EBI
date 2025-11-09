import plotly.graph_objects as go

ax_style = dict(
    showbackground=False,
    backgroundcolor="rgb(240, 240, 240)",
    showgrid=False,
    zeroline=False,
)


fig = go.Figure(
    go.Scatter3d(
        x=[0, 1], y=[-1, 0.2], z=[1, 2], mode="lines", line_color="red", line_width=2
    )
)
fig.update_layout(
    template="none",
    width=600,
    height=600,
    font_size=11,
    scene=dict(
        xaxis=ax_style,
        yaxis=ax_style,
        zaxis=ax_style,
        camera_eye=dict(x=1.85, y=1.85, z=1),
    ),
)

fig.show()
