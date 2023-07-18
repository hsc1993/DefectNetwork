
import plotly.graph_objects as go
import numpy as np
import ast
from plotly.graph_objs import *
import plotly
from plotly import express


if __name__ == "__main__":

    color_file = open('orthogonal_merged20_voro_color.txt', "r")
    color_data_read = color_file.read()
    color_data_lit = list(ast.literal_eval(color_data_read))
    color_data = np.array(color_data_lit)

    xlow = 10
    xhigh = 100
    ylow = 10
    yhigh = 100
    zlow = 10
    zhigh = 100

    X, Y, Z = np.mgrid[-1:1:20j, -1:1:20j, -1:1:10j]
    X, Y, Z = np.mgrid[xlow:xhigh:20j, ylow:yhigh:20j, zlow:zhigh:10j]
    X, Y, Z = np.mgrid[xlow:xhigh:30j, ylow:yhigh:30j, zlow:zhigh:20j]

    values_old = np.sin(np.pi * X) * np.cos(np.pi * Z) * np.sin(np.pi * Y)
    from scipy import ndimage

    color_data_reshape = color_data.reshape(20, 20, 10)
    margin_thickness = 5
    color_data_withmargin = np.zeros((20 + 2 * margin_thickness, 20 + 2 * margin_thickness, 10 + 2 * margin_thickness))

    for i in range(len(color_data_withmargin)):
        for j in range(len(color_data_withmargin[0])):
            for k in range(len(color_data_withmargin[0, 0])):
                if i < margin_thickness or j < margin_thickness or k < margin_thickness:
                    color_data_withmargin[i, j, k] = 0
                elif i > 19 + margin_thickness or j > 19 + margin_thickness or k > 9 + margin_thickness:
                    color_data_withmargin[i, j, k] = 0
                else:
                    color_data_withmargin[i, j, k] = color_data_reshape[
                        i - margin_thickness, j - margin_thickness, k - margin_thickness]

    color_data_withmargin = color_data_withmargin.flatten()
    print(color_data_withmargin.shape)
    print(X.shape)
    vol = ndimage.gaussian_filter(color_data_withmargin, 0.8)
    vol = ndimage.gaussian_filter(vol, 0.9)
    print(vol.shape)
    print(vol.min())

    # vol /= vol.max()

    layout1 = Layout(
        # paper_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='black',
        plot_bgcolor='rgba(0,0,0,0)',
    )

    print(max(vol))

    colorscale = plotly.express.colors.get_colorscale('jet')
    newjet = plotly.express.colors.colorscale_to_colors(colorscale)
    print('newjet=', newjet)

    opacity_value = 1
    fig = go.Figure(data=[

        go.Volume(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            # value=values_old.flatten(),
            # value=color_data,
            value=vol,
            # flatshading=True,
            # isomin=0.1,
            # isomax=0.3,
            isomin=3e16,
            isomax=3e17,
            opacity=0.1,  # needs to be small to see through all surfaces
            surface_count=21,  # needs to be a large number for good volume rendering
            colorscale=newjet,
            # colorscale='jet',
            # colorscale=["red", "green", "blue"],

            cmin=3e16,
            cmax=3e17,
            colorbar=dict(
                # tickvals=np.array([0.5e17,1e17,1.5e17,2e17]),
                thickness=30,
                tickfont=dict(color='white', size=30),
                title=dict(text='GND'),
                tickwidth=3,
                tickcolor='white',
                bordercolor='black',
                # tickformatstops=dict(enabled=True)

            ),

        ),
        # GB surface
        go.Surface(
            x=[25, 88, 88, 27, 25],
            y=[26, 85, 85, 27, 26],
            z=[
                [28, 28],
                [68.5, 68.5],
                [68, 68],
                [28,28],
                [28, 28]
            ],
            opacity=opacity_value,
            colorscale='cividis',
            cmin=0,
            cmax=1,
            surfacecolor=[[0.4, 0.4], [0.4, 0.4], [0.4, 0.4], [0.4, 0.4], [0.4, 0.4]],
            showscale=False
        ),
        # GB side surface
        go.Surface(
            x=[88, 88, 88, 88],
            y=[26, 27, 85, 85],
            z=[
                [28, 28],
                [28, 28.5],
                [68, 68.5],
                [68, 68.5]
            ],
            opacity=opacity_value,
            colorscale='cividis',
            cmin=0,
            cmax=1,
            surfacecolor=[[0.4, 0.4], [0.4, 0.4], [0.4, 0.4], [0.4, 0.4]],
            showscale=False
        ),
    ],
        layout=layout1
    )

    fig.update_layout(scene=dict(
        xaxis=dict(
            # backgroundcolor="rgb(200, 200, 230)",
            backgroundcolor="black",
            gridcolor="rgb(200, 200, 230)",
            showbackground=True,
            zerolinecolor="white",
            showgrid=False,
            visible=False,

        ),
        yaxis=dict(
            backgroundcolor="rgb(230, 200,230)",
            gridcolor="rgb(200, 200, 230)",
            showbackground=False,
            zerolinecolor="white",
            showgrid=False,
            visible=False

        ),

        zaxis=dict(
            backgroundcolor="rgb(230, 230,200)",
            gridcolor="rgb(200, 200, 230)",
            showbackground=False,
            zerolinecolor="white",
            showgrid=False,
            visible=False
        ), ),
        # width=700,
        margin=dict(
            r=10, l=10,
            b=10, t=10),

    )


    fig.show()

main()

