import dash
from dash import dcc, html, dash_table, callback
from dash.dependencies import Input, Output, State


dash.register_page(__name__)
layout = html.Div(children=[
    html.H1(children='This is our Prediction page'),

    html.Div(children='''
        This is our Prediction page content.
    '''),

])
