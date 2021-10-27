# from typing import Any, Dict

import dash
import dash_bootstrap_components as dbc
import dash_daq as daq
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import dash_table, dcc, html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash_extensions import Download
from dash_extensions.snippets import send_data_frame
from plotly import graph_objects as go

from app import app
from functions import df_to_table, get_assembly_stats, load_markers, marker_size_scaler

# JSONDict = Dict[str, Any]
colors = {"background": "#F3F6FA", "background_div": "white"}

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SKETCHY])
app.config.suppress_callback_exceptions = True

# tab 2 mag refinements elements
hidden_div_refinements_clusters = html.Div(
    id="refinements-clusters", style={"display": "none"}
)
binning_overview_2D = dbc.Row(
    [
        dbc.Col(
            [
                html.H5("Figure 1: 2D Binning Overview"),
                dcc.Graph(
                    id="scatterplot-2d",
                    style={"height": "90%", "width": "98%"},
                    clear_on_unhover=True,
                ),
            ],
            width=9,
        ),
        dbc.Col(
            [
                dbc.Button(
                    "Download Refinements",
                    id="refinements-download-button",
                    n_clicks=0,
                    color="primary",
                ),
                Download(id="refinements-download"),
                html.H3("Color contigs by:"),
                dcc.Dropdown(
                    id="color-by-column",
                    options=[],
                    value="cluster",
                    clearable=False,
                ),
                html.H3("X-Axis:"),
                dcc.Dropdown(
                    id="x-axis-2d",
                    options=[
                        {"label": "X_1", "value": "x_1"},
                        {"label": "Coverage", "value": "coverage"},
                        {"label": "GC%", "value": "gc_content"},
                        {"label": "Length", "value": "length"},
                    ],
                    value="x_1",
                    clearable=False,
                ),
                html.H3("Y-Axis:"),
                dcc.Dropdown(
                    id="y-axis-2d",
                    options=[
                        {"label": "X_2", "value": "x_2"},
                        {"label": "Coverage", "value": "coverage"},
                        {"label": "GC%", "value": "gc_content"},
                        {"label": "Length", "value": "length"},
                    ],
                    value="x_2",
                    clearable=False,
                ),
                html.Pre(
                    style={
                        "textAlign": "middle",
                    },
                    id="selection-binning-metrics",
                ),
            ],
            width=3,
        ),
    ]
)
binning_overview_3D_pie_chart = dbc.Row(
    [
        dbc.Col(
            [
                html.Label("Fig. 2: Change Z-axis"),
                dcc.Dropdown(
                    id="scatterplot-3d-zaxis-dropdown",
                    options=[
                        {"label": "Coverage", "value": "coverage"},
                        {"label": "GC%", "value": "gc_content"},
                        {"label": "Length", "value": "length"},
                    ],
                    value="coverage",
                    clearable=False,
                ),
                dcc.Graph(
                    id="scatterplot-3d",
                    clear_on_unhover=True,
                    style={"height": "90%", "width": "98%"},
                    config={
                        "toImageButtonOptions": dict(
                            format="svg",
                            filename="scatter3dPlot.autometa.binning",
                        ),
                    },
                ),
            ],
            width=6,
        ),
        dbc.Col(
            [
                html.Label("Fig. 3: Distribute Taxa by Rank"),
                dcc.Dropdown(
                    id="taxonomy-piechart-dropdown",
                    options=[
                        {"label": "Kingdom", "value": "superkingdom"},
                        {"label": "Phylum", "value": "phylum"},
                        {"label": "Class", "value": "class"},
                        {"label": "Order", "value": "order"},
                        {"label": "Family", "value": "family"},
                        {"label": "Genus", "value": "genus"},
                        {"label": "Species", "value": "species"},
                    ],
                    value="superkingdom",
                    clearable=False,
                ),
                html.Label("Figure 3: Taxonomic Distribution"),
                dcc.Graph(
                    id="taxonomy-piechart",
                    style={"height": "90%", "width": "98%"},
                    config=dict(displayModeBar=False),
                ),
            ],
            width=6,
        ),
    ]
)
toggles = dbc.Row(
    [
        dbc.Col(
            # add show-legend-toggle
            daq.ToggleSwitch(
                id="show-legend-toggle",
                size=40,
                color="#c5040d",
                label="Show/Hide 2D Legend",
                labelPosition="top",
                vertical=False,
                value=True,
            ),
            width=3,
            class_name="mt-3 mb-3 border",
        ),
        dbc.Col(
            # add hide selection toggle
            daq.ToggleSwitch(
                id="hide-selections-toggle",
                size=40,
                color="#c5040d",
                label="Hide Selections",
                labelPosition="top",
                vertical=False,
                value=False,
            ),
            width=3,
        ),
        dbc.Col(
            # add save selection toggle
            daq.ToggleSwitch(
                id="save-selections-toggle",
                size=40,
                color="#c5040d",
                label="Store Selections",
                labelPosition="top",
                vertical=False,
                value=False,
            ),
            width=3,
        ),
        dbc.Col(
            # add scatterplot-3d-legend-toggle
            daq.ToggleSwitch(
                id="scatterplot-3d-legend-toggle",
                size=40,
                color="#c5040d",
                label="Show/Hide 3D scatterplot Legend",
                labelPosition="top",
                vertical=False,
                value=True,
            ),
            width=3,
        ),
    ]
)
note = dbc.Row(
    dbc.Col(
        html.P(
            "Toggling save with contigs selected will save them as a refinement group."
        )
    )
)
table = dbc.Row(dbc.Col(dbc.Table(id="refinements-table"), width=12))

mag_file_upload = dcc.Upload(
    id="upload-mag-file",
    children=dbc.Button("upload a MAG file", color="primary"),
)
marker_file_upload = dcc.Upload(
    id="upload-markers-file",
    children=dbc.Button("upload a markers file", color="secondary"),
)
fasta_file_upload = dcc.Upload(
    id="upload-fasta-file",
    children=dbc.Button("upload a fasta file", color="secondary"),
)
# tab 1 upload files
tab1 = dbc.Tab(
    label="Upload MAG and Summary Files",
    id="tab-1-upload-files",
    children=[
        dbc.Row(
            [
                dbc.Col(mag_file_upload, width=3),
                dbc.Col(marker_file_upload, width=3),
                dbc.Col(fasta_file_upload, width=3),
                # TODO: do we need to add a button to press "compute" after uploading
                # the required files? that way if the user does upload the optional
                # files, we take them into account. It would also hold off before firing
                # all the callbacks as soon as    
                # something is uploaded.
                # dbc.Col(dbc.Button("Compute", color='success', id='ok_to_compute'))
            ],
            align="center",
            class_name="mt-3 mb-3",
        )
    ],
)

# mag refinement
tab2 = dbc.Tab(
    label="MAG Refinement",
    tab_id="tab-2-mag-refinement",
    children=[
        hidden_div_refinements_clusters,
        toggles,
        binning_overview_2D,
        binning_overview_3D_pie_chart,
        table,
    ],
)

# summary
# TODO: add children after declaring them as separate objects similar to tab2 above
tab3 = dbc.Tab(label="MAG Summary", tab_id="tab-3-mag-summary", children=[])

app.layout = dbc.Container(dbc.Tabs([tab1, tab2, tab3]), fluid=True)


def indicator(color, text, id_value):
    return html.Div(
        [
            html.P(text, className="twelve columns indicator_text"),
            html.Pre(id=id_value, className="indicator_value"),
        ],
        className="two columns indicator",
    )


def bin_dropdown(df, column):
    options = [{"label": bin, "value": bin} for bin in df[column].unique()]
    return options


# TODO:
# def plot_pie_chart(taxonomy: JSONDict, rank: str) -> Dict:
# removing the typing since it's not consistent all through and doesn't add much IMO(?)
#
# can this function be taken out of here? Additionally, is it possible to return a figure directly instead of a dict?
# example:
# fig = go.Figure()
# fig.add_trace(go.Pie(...))
# fig.update_layout(go.Layout(...))
# return fig
def plot_pie_chart(taxonomy, rank):
    df = pd.read_json(taxonomy, orient="split")
    total_contigs = df.shape[0]
    values = [
        contig / total_contigs for contig in df.groupby(rank)[rank].count().tolist()
    ]
    labels = df.groupby(rank)[rank].count().index.tolist()
    layout = go.Layout(
        margin=dict(l=0, r=0, b=0, t=4, pad=8),
        legend=dict(orientation="h"),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    trace = go.Pie(
        labels=labels,
        values=values,
        hoverinfo="label+percent",
        textinfo="label",
        showlegend=False,
    )
    return {"data": [trace], "layout": layout}


def taxa_by_rank(df, column, rank):
    clusters = dict(list(df.groupby(column)))
    clusters = df[column].unique().tolist()
    clusters.pop(clusters.index("unclustered"))
    nuniques = [df[df[column] == cluster][rank].nunique() for cluster in clusters]
    data = [
        go.Bar(y=clusters, x=nuniques, orientation="h")
    ]  # x could be any column value since its a count

    layout = go.Layout(
        barmode="stack",
        margin=dict(l=210, r=25, b=20, t=0, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return {"data": data, "layout": layout}


@app.callback(
    Output("color-by-column", "options"), [Input("metagenome-annotations", "children")]
)
def get_color_by_cols(annotations):
    df = pd.read_json(annotations, orient="split")
    return [
        {"label": col.title().replace("_", " "), "value": col}
        for col in df.columns
        if df[col].dtype.name not in {"float64", "int64"} and col != "contig"
    ]


@app.callback(
    Output("selection-binning-metrics", "children"),
    [
        Input("kingdom-markers", "children"),
        Input("scatterplot-2d", "selectedData"),
    ],
)
def display_selection_summary(markers, selected_contigs):
    if not selected_contigs:
        num_expected_markers = pd.NA
        num_markers_present = pd.NA
        total_markers = pd.NA
        n_marker_sets = pd.NA
        completeness = pd.NA
        purity = pd.NA
        n_selected = 0
    else:
        df = pd.read_json(markers, orient="split").set_index("contig")
        contigs = {point["text"] for point in selected_contigs["points"]}
        n_selected = len(selected_contigs["points"])
        num_expected_markers = df.shape[1]
        pfam_counts = df.loc[df.index.isin(contigs)].sum()
        if pfam_counts.empty:
            total_markers = 0
            num_single_copy_markers = 0
            num_markers_present = 0
            completeness = pd.NA
            purity = pd.NA
            n_marker_sets = pd.NA
        else:
            total_markers = pfam_counts.sum()
            num_single_copy_markers = pfam_counts.eq(1).count()
            num_markers_present = pfam_counts.ge(1).count()
            completeness = num_markers_present / num_expected_markers * 100
            purity = num_single_copy_markers / num_markers_present * 100
            n_marker_sets = total_markers / num_expected_markers
    # TODO: Create cleaner table for Sam to read completeness/purity etc.
    return f"""
    Selection Binning Metrics:
    -----------------------
|    Markers Expected:\t{num_expected_markers}\t|
|    Unique Markers:\t{num_markers_present}\t|
|    Total Markers:\t{total_markers:,}\t|
|    Marker Set(s):\t{n_marker_sets:.02f}\t|
|    Completeness:\t{completeness:.02f}\t|
|    Purity:\t\t{purity:.02f}\t|
|    Contigs Selected:\t{n_selected:,}\t|
    -----------------------
    """


@app.callback(
    Output("taxonomy-piechart", "figure"),
    [
        Input("metagenome-annotations", "children"),
        Input("scatterplot-2d", "selectedData"),
        Input("taxonomy-piechart-dropdown", "value"),
    ],
)
def taxa_piechart(annotations, selected_contigs, selected_rank):
    df = pd.read_json(annotations, orient="split")
    layout = dict(margin=dict(l=15, r=10, t=35, b=45))
    if not selected_contigs:
        n_ctgs = len(df.index)
        # Get taxa and their respective count at selected canonical-rank
        labels = df[selected_rank].unique().tolist()
        values = [
            len(df[df[selected_rank] == label]) / float(n_ctgs) for label in labels
        ]
        # Add in to pie chart... Sankey (go.Sankey) diagram or Parallel Categories plot (go.ParCat)
        trace = go.Pie(
            labels=labels,
            values=values,
            textinfo="label",
            hoverinfo="label+percent",
            showlegend=False,
        )
        return dict(data=[trace], layout=layout)

    ctg_list = {point["text"] for point in selected_contigs["points"]}
    dff = df[df.contig.isin(ctg_list)]
    n_ctgs = len(dff.index)
    labels = dff[selected_rank].unique().tolist()
    values = [len(dff[dff[selected_rank] == label]) / float(n_ctgs) for label in labels]
    trace = go.Pie(
        labels=labels,
        values=values,
        hoverinfo="label+percent",
        textinfo="label",
        showlegend=False,
    )
    return dict(data=[trace], layout=layout)


@app.callback(
    Output("scatterplot-3d", "figure"),
    [
        Input("metagenome-annotations", "children"),
        Input("scatterplot-3d-zaxis-dropdown", "value"),
        Input("scatterplot-3d-legend-toggle", "value"),
        Input("color-by-column", "value"),
        Input("scatterplot-2d", "selectedData"),
    ],
)
def update_zaxis(annotations, zaxis, show_legend, colorby_col, selected_contigs):
    df = pd.read_json(annotations, orient="split")
    titles = {
        "coverage": "Coverage",
        "gc_content": "GC Content",
        "length": "Length",
    }
    zaxis_title = titles[zaxis]
    if not selected_contigs:
        contigs = df.contig.tolist()
    else:
        contigs = {point["text"] for point in selected_contigs["points"]}
    # Subset DataFrame by selected contigs
    df = df[df.contig.isin(contigs)]
    if colorby_col == "cluster":
        # Categoricals for binning
        df[colorby_col] = df[colorby_col].fillna("unclustered")
    else:
        # Other possible categorical columns all relate to taxonomy
        df[colorby_col] = df[colorby_col].fillna("unclassified")
    return {
        "data": [
            go.Scatter3d(
                x=dff.x_1,
                y=dff.x_2,
                z=dff[zaxis],
                text=dff.contig,
                mode="markers",
                textposition="top center",
                opacity=0.45,
                hoverinfo="all",
                marker={
                    "size": dff.assign(normLen=marker_size_scaler)["normLen"].fillna(1),
                    "line": {"width": 0.1, "color": "black"},
                },
                name=colorby_col_value,
            )
            for colorby_col_value, dff in df.groupby(colorby_col)
        ],
        "layout": go.Layout(
            scene=dict(
                xaxis=dict(title="X_1"),
                yaxis=dict(title="X_2"),
                zaxis=dict(title=zaxis_title),
            ),
            legend={"x": 0, "y": 1},
            showlegend=show_legend,
            autosize=True,
            margin=dict(r=0, b=0, l=0, t=25),
            hovermode="closest",
        ),
    }


@app.callback(
    Output("scatterplot-2d", "figure"),
    [
        Input("metagenome-annotations", "children"),
        Input("x-axis-2d", "value"),
        Input("y-axis-2d", "value"),
        Input("show-legend-toggle", "value"),
        Input("color-by-column", "value"),
        Input("refinements-clusters", "children"),
        Input("hide-selections-toggle", "value"),
    ],
)
def update_axes(
    annotations,
    xaxis_column,
    yaxis_column,
    show_legend,
    cluster_col,
    refinement,
    hide_selection_toggle,
):
    df = pd.read_json(annotations, orient="split").set_index("contig")
    titles = {
        "x_1": "X_1",
        "x_2": "X_2",
        "coverage": "Coverage",
        "gc_content": "GC Content",
        "length": "Length",
    }
    xaxis_title = titles[xaxis_column]
    yaxis_title = titles[yaxis_column]
    # Subset metagenome-annotations by selections iff selections have been made
    df[cluster_col] = df[cluster_col].fillna("unclustered")
    if hide_selection_toggle:
        refine_df = pd.read_json(refinement, orient="split").set_index("contig")
        refine_cols = [col for col in refine_df.columns if "refinement" in col]
        if refine_cols:
            refine_col = refine_cols.pop()
            # Retrieve only contigs that have already been refined...
            refined_contigs_index = refine_df[
                refine_df[refine_col].str.contains("refinement")
            ].index
            # print(f"refined df shape: {refine_df.shape}, df shape: {df.shape}")
            df.drop(refined_contigs_index, axis="index", inplace=True, errors="ignore")
            # print(f"new df shape: {df.shape}")
    return {
        "data": [
            go.Scattergl(
                x=df[df[cluster_col] == cluster][xaxis_column],
                y=df[df[cluster_col] == cluster][yaxis_column],
                text=df[df[cluster_col] == cluster].index,
                mode="markers",
                opacity=0.45,
                marker={
                    "size": df.assign(normLen=marker_size_scaler)["normLen"],
                    "line": {"width": 0.1, "color": "black"},
                },
                name=cluster,
            )
            for cluster in df[cluster_col].unique()
        ],
        "layout": go.Layout(
            scene=dict(
                xaxis=dict(title=xaxis_title),
                yaxis=dict(title=yaxis_title),
            ),
            legend={"x": 1, "y": 1},
            showlegend=show_legend,
            margin=dict(r=50, b=50, l=50, t=50),
            # title='2D Clustering Visualization',
            hovermode="closest",
        ),
    }


@app.callback(
    Output("datatable", "data"),
    [
        Input("scatterplot-2d", "selectedData"),
        Input("refinements-clusters", "children"),
    ],
)
def update_table(selected_data, refinements):
    df = pd.read_json(refinements, orient="split")
    if not selected_data:
        return df.to_dict("records")
    contigs = {point["text"] for point in selected_data["points"]}
    return df[df.contig.isin(contigs)].to_dict("records")


@app.callback(
    Output("refinements-table", "children"),
    [Input("refinements-clusters", "children")],
)
def bin_table(df):
    df = pd.read_json(df, orient="split")
    return dash_table.DataTable(
        id="datatable",
        data=df.to_dict("records"),
        columns=[{"name": col, "id": col} for col in df.columns],
        style_cell={"textAlign": "center"},
        style_cell_conditional=[{"if": {"column_id": "contig"}, "textAlign": "right"}],
        virtualization=True,
    )


@app.callback(
    Output("refinements-download", "data"),
    [
        Input("refinements-download-button", "n_clicks"),
        Input("refinements-clusters", "children"),
    ],
)
def download_refinements(n_clicks, intermediate_selections):
    if not n_clicks:
        raise PreventUpdate
    df = pd.read_json(intermediate_selections, orient="split")
    return send_data_frame(df.to_csv, "refinements.csv", index=False)


@app.callback(
    Output("refinements-clusters", "children"),
    [
        Input("scatterplot-2d", "selectedData"),
        Input("refinement-data", "children"),
        Input("save-selections-toggle", "value"),
    ],
    [State("refinements-clusters", "children")],
)
def store_binning_refinement_selections(
    selected_data, refinement_data, save_toggle, intermediate_selections
):
    if not selected_data and not intermediate_selections:
        # We first load in our binning information for refinement
        # Note: this callback should trigger on initial load
        # TODO: Could also remove and construct dataframes from selected contigs
        # Then perform merge when intermediate selections are downloaded.
        bin_df = pd.read_json(refinement_data, orient="split")
        bin_df["cluster"].fillna("unclustered", inplace=True)
        return bin_df.to_json(orient="split")
    if not save_toggle or not selected_data:
        raise PreventUpdate
    contigs = {point["text"] for point in selected_data["points"]}
    pdf = pd.read_json(intermediate_selections, orient="split").set_index("contig")
    refinement_cols = [col for col in pdf.columns if "refinement" in col]
    refinement_num = len(refinement_cols) + 1
    group_name = f"refinement_{refinement_num}"
    pdf.loc[contigs, group_name] = group_name
    pdf = pdf.fillna(axis="columns", method="ffill")
    pdf.reset_index(inplace=True)
    return pdf.to_json(orient="split")


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', debug=True)