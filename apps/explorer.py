# -*- coding: utf-8 -*-


import os
import sys
import json
import math
import base64
import datetime
import io

import pandas as pd
import numpy as np
import flask

import dash
# import dash_cytoscape as cyto
import dash_table
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from plotly import graph_objs as go
import plotly.plotly as py

from app import app, indicator, millify, df_to_table, parse_df_upload, parse_contents

normalizeLen = lambda x: np.ceil(
    (x.length-x.length.min()) / (x.length.max()-x.length.min())
) * 2 + 4

# def modal():
#     return html.Div(
#         html.Div(
#             [
#                 html.Div(
#                     [
#                         # modal header
#                         html.Div(
#                             [
#                                 html.Span(
#                                     "Export Selections",
#                                     style={
#                                         "color": "#c5040d",
#                                         "fontWeight": "bold",
#                                         "fontSize": "20",
#                                     },
#                                 ),
#                                 html.Span(
#                                     "×",
#                                     id="analysis_modal_close",
#                                     n_clicks=0,
#                                     style={
#                                         "float": "right",
#                                         "cursor": "pointer",
#                                         "marginTop": "0",
#                                         "marginBottom": "17",
#                                     },
#                                 ),
#                             ],
#                             className="row",
#                             style={"borderBottom": "1px solid #C8D4E3"},
#                         ),
#                         # modal form
#                         html.Div(
#                             [
#                                 html.P(
#                                     [
#                                         "Length Cutoff",
#                                     ],
#                                     style={
#                                         "float": "left",
#                                         "marginTop": "4",
#                                         "marginBottom": "2",
#                                     },
#                                     className="row",
#                                 ),
#                                 dcc.Input(
#                                     id="length_cutoff",
#                                     placeholder="Length Cutoff (default 3000bp)",
#                                     type="text",
#                                     value="3000",
#                                     style={"width": "100%"},
#                                 ),
#                                 html.P(
#                                     style={
#                                         "textAlign": "left",
#                                         "marginBottom": "2",
#                                         "marginTop": "4",
#                                     },
#                                     id="completeness_display"
#                                 ),
#                                 dcc.Slider(
#                                     id="completeness_cutoff",
#                                     min=5.0,
#                                     max=100.0,
#                                     value=20.0,
#                                     updatemode="drag",
#                                 ),
#                                 html.P(
#                                     "Select Metagenome Assembly",
#                                     style={
#                                         "textAlign": "left",
#                                         "marginBottom": "2",
#                                         "marginTop": "4",
#                                     },
#                                 ),
#                                 dcc.Upload(
#                                     id='upload-data',
#                                     children=['Drag and Drop or ', html.A('Select a File')],
#                                     style={
#                                         'width': '100%',
#                                         'height': '60px',
#                                         'lineHeight':' 60px',
#                                         'borderWidth': 'dashed',
#                                         'borderRadius': '5px',
#                                         'textAlign': 'center',
#                                     },
#
#                                 )
#                             ],
#                             className="row",
#                             style={"padding": "2% 8%"},
#                         ),
#                         # submit button
#                         html.Span(
#                             "Submit",
#                             id="submit_new_lead",
#                             n_clicks=0,
#                             className="button button--primary add"
#                         ),
#                     ],
#                     className="modal-content",
#                     style={"textAlign": "center"},
#                 )
#             ],
#             className="modal",
#         ),
#         id="analysis_modal",
#         style={"display": "none"},
#     )


layout = [
    # 2D-scatter plot row div
    html.Div(
            [
                html.Div(
                    [
                        html.P("2D Binning Overview"),
                        dcc.Graph(
                            id='scatter2d_graphic',
                            style={"height": "90%", "width": "98%"},
                            clear_on_unhover=True,
                        ),
                    ],
                    className="ten columns chart_div"
                ),
                html.Div([
                    html.P("Color By:"),
                    dcc.Dropdown(
                        id="cluster_col",
                        options=[],
                        value="cluster",
                        clearable=False,
                    ),
                    html.P("X-Axis:"),
                    dcc.Dropdown(
                        id="2d_xaxis",
                        options=[
                            {'label':'bh-tsne-x','value':'bh_tsne_x'},
                            {'label':'Coverage','value':'cov'},
                            {'label':'GC%','value':'gc'},
                            {'label':'Length','value':'length'},
                        ],
                        value="bh_tsne_x",
                        clearable=False,
                    ),
                    html.P("Y-Axis:"),
                    dcc.Dropdown(
                        id="2d_yaxis",
                        options=[
                            {'label':'bh-tsne-y','value':'bh_tsne_y'},
                            {'label':'Coverage','value':'cov'},
                            {'label':'GC%','value':'gc'},
                            {'label':'Length','value':'length'},
                        ],
                        value="bh_tsne_y",
                        clearable=False,
                    ),
                    # html.Pre('Legend:\nCircle Colors:'),
                    daq.ToggleSwitch(id='2d-legend', label='Legend', labelPosition='right', size=40),
                    daq.ToggleSwitch(id='2d-circles', label='Circle bins', labelPosition='right', size=40),
                    html.Pre(
                        style={
                            "textAlign": "middle",
                        },
                        id="selection_summary"
                    ),
                ],
                className='two columns',
                ),
            ],
        ),
    # 3D-scatter plot row div
    html.Div(
        [
            html.Div(
                [
                    html.P("3D Binning Overview"),
                    dcc.Graph(
                        id='scatter3d_graphic',
                        clear_on_unhover=True,
                        style={"height": "90%", "width": "98%"},
                        config={
                            'toImageButtonOptions':dict(
                                format='svg',
                                filename='scatter3dPlot.autometa.binning',
                            ),
                            # 'displayModeBar':False,
                        }
                    ),
                ],
                className="ten columns chart_div"
            ),
            html.Div([
                # add button
                # html.Span(
                #     html.A(id="save_link", children="Save Selected"),
                #     id="save_button",
                #     n_clicks=0,
                #     className="button button--primary",
                #     style={
                #         "height": "34",
                #         "background": "#c5040d",
                #         "border": "1px solid #c5040d",
                #         "color": "white",
                #     },
                # ),
                html.P("Z-Axis:"),
                dcc.Dropdown(
                    id="zaxis_column",
                    options=[
                        {"label": "Coverage", "value": "cov"},
                        {"label": "GC%", "value": "gc"},
                        {"label": "Length", "value": "length"},
                    ],
                    value="cov",
                    clearable=False,
                ),
                daq.ToggleSwitch(
                    id='3d-legend',
                    label='Legend',
                    labelPosition='top',
                ),
            ],
            className="two columns",
            ),
        ],
        className="row",
        style={"marginTop": "0", "marginBottom": "2"},
    ),
    # Taxa Distribution plot row div
    html.Div(
        [
        html.Div(
            [
                html.P("Taxa Distribution"),
                dcc.Graph(
                    id='taxa_piechart',
                    style={"height": "90%", "width": "98%"},
                    config=dict(displayModeBar=False),
                ),
            ],
            className="ten columns chart_div",
        ),
        html.Div(
            [
                html.P("Rank:"),
                dcc.Dropdown(
                    id="rank_dropdown",
                    options=[
                        {'label':'Kingdom','value':'kingdom'},
                        {'label':'Phylum','value':'phylum'},
                        {'label':'Class','value':'class'},
                        {'label':'Order','value':'order'},
                        {'label':'Family','value':'family'},
                        {'label':'Genus','value':'genus'},
                        {'label':'Species','value':'species'},
                    ],
                    value="phylum",
                    clearable=False,
                ),
            ],
            className="row two columns",
        ),
    ]),
    # table div
    html.Div(
        className="row twelve columns",
        style={
            "maxHeight": "350px",
            "overflowY": "scroll",
            "padding": "8",
            "marginTop": "5",
            "backgroundColor":"white",
            "border": "1px solid #C8D4E3",
            "borderRadius": "3px"
        },
        id='binning_table'
    ),
    html.Div(id='uploaded-data'),

    # modal(),
]

# # TODO: See https://plot.ly/python/custom-buttons/ and relayout to remove
# legend and draw circles (or convex hulls) over clusters
# @app.callback(
#     Output("2d-figure", ""),
#     [Input("2d-legend", "value")]
# )
# def draw_2d_circles(value):
#     return circles

# # TODO: DEBUG: Legend Hide/Show See line : 495
# Hide legend
@app.callback(
    Output('scatter2d_graphic','layout'),
    [Input('2d-legend','value')]
)
def hide_2d_legend(value):
    return dict(showlegend=value)

@app.callback(
    Output("cluster_col", "options"),
    [Input("binning_df", "children")]
)
def get_color_by_cols(df):
    df = pd.read_json(df, orient='split')
    df.drop(axis=1,labels=["single_copy_PFAMs"],inplace=True)
    options = [{"label":col.title().replace('_',' '), "value":col}
        for col in df.columns if df[col].dtype.name not in {'float64', 'int64'}]
    return options

@app.callback(
    Output("selection_summary", "children"),
    [Input("scatter2d_graphic", "selectedData"),
    Input("cluster_col", "value"),
    Input("binning_df", "children")]
)
def display_selection_summary(selectedData, clusterCol, df):
    if not selectedData:
        lines =[
            'Select points to update indicators',
            'Marker Set: N/A',
            'Completeness: N/A',
            'Purity: N/A',
            'Selected: 0',
        ]
        return '\n'.join(lines)
    # Update to return total number of bins
    df = pd.read_json(df, orient='split')
    ctg_list = {point['text'] for point in selectedData['points']}
    n_selected = len(selectedData['points'])
    pfams = df[df.contig.isin(ctg_list)].single_copy_PFAMs.dropna().tolist()
    all_pfams = [p for pfam in pfams for p in pfam.split(',')]
    total = len(all_pfams)
    markers = 139
    n_marker_sets = round(float(total)/markers, 2)
    marker_sets = '{} ({} markers)'.format(n_marker_sets, total)
    nunique = len(set(all_pfams))
    completeness = round(float(nunique)/markers * 100, 2)
    purity = '-' if total == 0 else round(float(nunique)/total * 100, 2)
    indicators = ['Marker Set(s)','Completeness','Purity','Selected']
    lines = map(str,[
        marker_sets,
        completeness,
        purity,
        n_selected,
    ])
    p = []
    for i,l in zip(indicators, lines):
        p.append(': '.join([i,l]))
    return '\n'.join(p)

@app.callback(
    Output("binning_table", "children"),
    [Input("binning_df", "children")],
)
def bin_table_callback(df):
    df = pd.read_json(df, orient="split")
    child = dash_table.DataTable(
        id='datatable',
        data=df.to_dict('records'),
        columns=[{'name':col, 'id':col} for col in df.columns],
        # sorting=True,
        # filtering=True,
        virtualization=True,
        # pagination_mode=None,
    ),
    return child

@app.callback(
    Output("taxa_piechart", "figure"),
    [Input("scatter2d_graphic", "selectedData"),
    Input("rank_dropdown", "value"),
    Input("binning_df", "children")]
)
def taxa_piechart_callback(selectedData, selectedRank, df):
    df = pd.read_json(df, orient="split")
    layout = dict(margin=dict(l=15, r=10, t=35, b=45))
    if not selectedData:
        n_ctgs = len(df.index)
        labels = df[selectedRank].unique().tolist()
        values = [len(df[df[selectedRank] == label])/float(n_ctgs) for label in labels]
        trace = go.Pie(
            labels=labels,
            values=values,
            textinfo='label',
            hoverinfo='label+percent',
            showlegend=False,
        )
        return dict(data=[trace], layout=layout)

    ctg_list = {point['text'] for point in selectedData['points']}
    dff = df[df.contig.isin(ctg_list)]
    n_ctgs = len(dff.index)
    labels = dff[selectedRank].unique().tolist()
    values = [len(dff[dff[selectedRank] == label])/float(n_ctgs) for label in labels]
    trace = go.Pie(
        labels=labels,
        values=values,
        hoverinfo='label+percent',
        textinfo='label',
        showlegend=False,
    )
    return dict(data=[trace], layout=layout)

@app.callback(
    Output('scatter3d_graphic', 'figure'),
    [Input('zaxis_column', 'value'),
    Input("cluster_col", "value"),
    Input('scatter2d_graphic','selectedData'),
    Input("binning_df", "children")])
def update_zaxis(zaxis_column, cluster_col, selectedData, df):
    df = pd.read_json(df, orient="split")
    titles = {
        'cov':'Coverage',
        'gc':'GC%',
        'length':'Length',
    }
    zaxis_title = titles[zaxis_column]
    if selectedData is None:
        selected = df.contig.tolist()
    else:
        selected = {point['text'] for point in selectedData['points']}

    dff = df[df.contig.isin(selected)]


    return {
        'data': [
            go.Scatter3d(
                x=dff[dff[cluster_col] == i]['bh_tsne_x'],
                y=dff[dff[cluster_col] == i]['bh_tsne_y'],
                z=dff[dff[cluster_col] == i][zaxis_column],
                text=dff[dff[cluster_col] == i]['contig'],
                mode='markers',
                textposition='top center',
                opacity=0.45,
                hoverinfo='all',
                marker={
                    'size': dff.assign(normLen = normalizeLen)['normLen'],
                    'line': {'width': 0.1, 'color': 'black'},
                },
                name=i
            ) for i in dff[cluster_col].unique()
        ],
        'layout': go.Layout(
            scene = dict(
                xaxis=dict(title='bh-tsne-x'),
                yaxis=dict(title='bh-tsne-y'),
                zaxis=dict(title=zaxis_title),
            ),
            legend={'x': 0, 'y': 1},
            showlegend=False,
            autosize=True,
            margin=dict(r=0, b=0, l=0, t=25),
            # title='3D Clustering Visualization',
            hovermode='closest',
        ),
    }
# OPTIMIZE: # DEBUG: # TODO: Figure should be defined in the div section above... Not in callbacks
# This needs to be refactored as an update layout approach
@app.callback(
    Output('scatter2d_graphic', 'figure'),
    [Input('2d_xaxis', 'value'),
    Input('2d_yaxis', 'value'),
    Input('cluster_col', 'value'),
    Input("binning_df", "children"),
    # Input('datatable', 'data')]
])
def update_axes(xaxis_column, yaxis_column, cluster_col, df):
    df = pd.read_json(df, orient="split")
    titles = {
        'bh_tsne_x':'bh-tsne-x',
        'bh_tsne_y':'bh-tsne-y',
        'cov':'Coverage',
        'gc':'GC%',
        'length':'Length',
    }
    xaxis_title = titles[xaxis_column]
    yaxis_title = titles[yaxis_column]

    return {
        'data': [
            go.Scattergl(
                x=df[df[cluster_col] == i][xaxis_column],
                y=df[df[cluster_col] == i][yaxis_column],
                text=df[df[cluster_col] == i]['contig'],
                mode='markers',
                opacity=0.45,
                marker={
                    'size': df.assign(normLen = normalizeLen)['normLen'],
                    'line': {'width': 0.1, 'color': 'black'},
                },
                name=i
            ) for i in df[cluster_col].unique()
        ],
        'layout': go.Layout(
            scene = dict(
                xaxis=dict(title=xaxis_title),
                yaxis=dict(title=yaxis_title),
            ),
            legend={'x': 1, 'y': 1},
            showlegend=False,
            margin=dict(r=50, b=50, l=50, t=50),
            # title='2D Clustering Visualization',
            hovermode='closest',
        ),
    }

@app.callback(
    Output('datatable', 'data'),
    [Input('scatter2d_graphic','selectedData'),
    Input("binning_df", "children")]
)
def update_table(selectedData, df):
    df = pd.read_json(df, orient="split")
    if selectedData is None:
        return df.to_dict('records')
    selected = {point['text'] for point in selectedData['points']}
    return df[df.contig.isin(selected)].to_dict('records')
