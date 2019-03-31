# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

from starkit.gridkit import load_grid


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

grid = load_grid('/home/nimish/Documents/phoenix_t4000_10000_w3000_9000_r3000.h5')

teff_low = grid.teff.bounds[0]
teff_high = grid.teff.bounds[1]



app.layout = html.Div([
    html.Div([

    html.Div([
        html.H4(
        "Paramters for 1st Spectra"),
        html.Label("Teff"),
        dcc.Slider(
            id='teff-slider',
            min=4000.0,
            max = 10000.0,
            step = 1,
            value=grid.teff.value
        ),
        html.Label("logg"),
        dcc.Slider(
            id='logg-slider',
            min=3.0,
            max=5.0,
            step = 0.01,
            value= grid.logg.value
        ),
        html.Label("mh"),
        dcc.Slider(
            id='mh-slider',
            min=-1.0,
            max=0.0,
            step = 0.01,
            value=-0.5,
        )
    ],style={'marginRight':100}),
    html.Div([
        html.H4(
        "Paramters for 2nd Spectra"),
        html.Label("Teff"),
        dcc.Slider(
            id='teff-slider2',
            min=4000.0,
            max = 10000.0,
            step = 1,
            value=grid.teff.value
        ),
        html.Label("logg"),
        dcc.Slider(
            id='logg-slider2',
            min=3.0,
            max=5.0,
            step = 0.01,
            value= grid.logg.value
        ),
        html.Label("mh"),
        dcc.Slider(
            id='mh-slider2',
            min=-1.0,
            max=0.0,
            step = 0.01,
            value=-0.5
        )
    ],style={'marginRight':100})


    ], style={'columnCount': 2}),

    dcc.Graph(id='spectra-graphic'),


],style={'margin':100})



@app.callback(
    Output('spectra-graphic', 'figure'),
    [Input('logg-slider', 'value'),
     Input('mh-slider', 'value'),
     Input('teff-slider', 'value'),
     Input('logg-slider2', 'value'),
     Input('mh-slider2', 'value'),
     Input('teff-slider2', 'value'),])
def update_graph(logg, mh,teff,logg2, mh2,teff2):

    grid.teff = teff
    grid.logg = logg
    grid.mh = mh
    wave,flux = grid()
    grid.teff = teff2
    grid.logg = logg2
    grid.mh = mh2
    wave2,flux2 = grid()
    wave2

    return {
        'data': [go.Scatter(
            x=wave,
            y=flux,
            # text=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'],
            mode='lines',
            name = 'teff = {}  logg = {}  mh = {}'.format(teff, logg, mh),

        ),
        go.Scatter(
            x=wave2,
            y=flux2,
            # text=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'],
            mode='lines',
            name = 'teff = {}  logg = {}  mh = {}'.format(teff2, logg2, mh2),

        )
        ],
        'layout': go.Layout(

            xaxis={
                'title': 'Wavelength(Armstrongs)',
                'type': 'linear'
            },
            yaxis={
                'title': 'Logarithmic Brightness',
                'type': 'linear'
            },
            margin={'l': 100, 'b': 40, 't': 200, 'r': 50},
            hovermode='closest',
            legend=dict(
                x=0,
                y=1.5,
                traceorder='normal',
                font=dict(
                    family='sans-serif',
                    size=14,
                    color='#000'
                ),
                bgcolor='#E2E2E2',
                bordercolor='#FFFFFF',
                borderwidth=2
    )
        )
    }
if __name__ == '__main__':
    app.run_server(debug=True)
