import dash
from dash import html, dcc, callback, Output, Input, State
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.express as px
import dash
from dash import dcc, html, dash_table, callback
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import umap.umap_ as umap
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.manifold import Isomap, LocallyLinearEmbedding
import plotly.graph_objs as go
import os
import diskcache

# Set up diskcache for caching
cache = diskcache.Cache("./cache")

# Caching duration (e.g., 60 seconds = 1 minute)
CACHE_DURATION = 180

dash.register_page(__name__, path='/', suppress_callback_exceptions=True)
data = pd.read_csv("METABRIC_RNA_Mutation3.csv")
available_variables = list(data.columns)
# Extract the numerical variables for correlation
numerical_data = data.select_dtypes(include='number')

EPSILON = 1e-9 

def linear_layout_algorithm(data, iterations=1000, learning_rate=0.1, Vmax=5, Smax=10):
    N, D = data.shape

    # Initialize positions and velocities randomly
    positions = np.random.random((N, 2))
    velocities = np.zeros((N, 2))

    for iteration in range(iterations):
        # Update neighbor sets
        neighbor_sets = []
        for i in range(N):
            distances = np.linalg.norm(data - data[i], axis=1)
            indices = np.argsort(distances)[1:Vmax+1]  
            neighbor_sets.append(indices)

        # Update positions and velocities
        for i in range(N):
            forces = np.zeros(2)
            for j in neighbor_sets[i]:
                dist_ij = np.linalg.norm(positions[j] - positions[i])
                dist_ij_high = np.linalg.norm(data[j] - data[i])
                force_ij = dist_ij - dist_ij_high
                if np.isfinite(force_ij):
                    forces += force_ij * (positions[j] - positions[i]) * np.divide(1, (dist_ij + EPSILON), where=(dist_ij > EPSILON))

            subset = np.random.choice(N, Smax, replace=False)
            for j in subset:
                dist_ij = np.linalg.norm(positions[j] - positions[i])
                dist_ij_high = np.linalg.norm(data[j] - data[i])
                force_ij = dist_ij - dist_ij_high
                if np.isfinite(force_ij):
                    forces += force_ij * (positions[j] - positions[i]) * np.divide(1, (dist_ij + EPSILON), where=(dist_ij > EPSILON))

            max_force_magnitude = np.max(np.linalg.norm(forces))  # Calculate norm without specifying axis
            normalized_forces = forces / max_force_magnitude  # Normalize forces
            velocities[i] += learning_rate * normalized_forces  # Update velocity with learning_rate
            positions[i] += velocities[i]  # Update position

        positions -= np.mean(positions, axis=0)  # Normalize positions
        std = np.std(positions, axis=0)
        positions = positions / std[np.newaxis, :] if not np.allclose(std, 0) else positions

    return positions



layout = html.Div(children=[

    html.Div(className='row', children=[
        html.Div( style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)'}, children=[
            html.H3('Statistics for each variable'),

            # Dropdown to select the variable for statistics table
            html.Label("Choose variable to filter:"),
            dcc.Dropdown(
                id='statistics-variable-dropdown',
                options=[{'label': var, 'value': var} for var in available_variables],
                value='',  # Set default value
            ),

            # Scrollable table container for statistics table
            html.Div(
                style={'height': '400px', 'overflowY': 'scroll', 'border': '1px solid #ccc',  # Add a border to the scrollable div
                        'borderRadius': '5px',      # Add border radius for a nicer look
                        'boxShadow': '0px 0px 5px rgba(0, 0, 0, 0.1)',  # Add a subtle shadow
                        'padding': '10px',          # Add some padding for spacing
                        'background': '#fff',  },
                children=[
                    html.Table(id='statistics-table', children=[
                        html.Thead(
                            html.Tr([
                                html.Th('Variable'),
                                html.Th('Number of Values'),
                                html.Th('Number of Gaps'),
                                html.Th('Number of Unique Values'),
                                html.Th('Number of 0 Values'),
                                html.Th('Most Common Value'),
                                html.Th('Min Value'),
                                html.Th('Max Value'),
                                html.Th('Median'),
                                html.Th('Mean'),
                                html.Th('Std Value'),
                                html.Th('Number of Outliers (2 σ)'),
                                html.Th('Number of Outliers (3 σ)'),
                                html.Th('Number of Outliers (4 σ)'),
                            ]), style={'position': 'sticky', 'top': '0', 'background-color': '#f0f0f0'},
                        ),
                        html.Tbody(id='statistics-table-body'),  
                    ]),
                ],
            ),
        ]),

    ], style={'margin-bottom': '20px'}),

    html.Div(className='row', children=[
        html.Div(className='six columns', style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)'}, children=[
            dcc.Graph(id='correlation-heatmap'),
        ]),

        html.Div(className='six columns', style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)'}, children=[
            # Dropdown to select the variable for correlation table
            html.H3('Pairwise correlation'),
            html.Label("Choose variable:"),
            dcc.Dropdown(
                id='variable-dropdown',
                options=[{'label': var, 'value': var} for var in numerical_data.columns],
                value=numerical_data.columns[5],  # Set default value
            ),

            # Scrollable table container for correlation table
            html.Div(
                style={'height': '570px', 'overflowY': 'scroll'},
                children=[
                    html.Table(id='correlation-table', children=[
                        html.Thead(
                            html.Tr([
                                html.Th('Variable'),
                                html.Th('Correlation'),
                                html.Th('Diagram'),
                            ])
                        ),
                        html.Tbody(id='correlation-table-body'),  # Include the correlation-table-body 
                    ]),
                ],
            ),
        ]),
    ], style={'margin-bottom': '20px'}),

    # Frequency Graph by Label section
    html.Div(style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)', 'margin-bottom': '20px'}, children=[
        html.Div(className='row', children=[
            # Dropdown container
            html.Div( children=[
                html.Label("Choose variable for frequency graph:"),
                dcc.Dropdown(
                    id='frequency-variable-dropdown',
                    options=[{'label': var, 'value': var} for var in data.columns[:31]],
                    value=data.columns[7],  # Set default value
                ),
            ]), 
            # Graph container
            html.Div( children=[
                html.H3('Frequency Graph by Label'),
                dcc.Graph(id='frequency-graph'),
            ]),
        ]),
    ]),

    # Scatter Plot and Scatter Matrix section
    html.Div(className='row', children=[
        html.Div( style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)'}, children=[
            html.H3('Scatter Plot'),

            # Dropdowns for X and Y axes
            html.Label("Choose variable for X-axis:"),
            dcc.Dropdown(
                id='scatter-x-variable-dropdown',
                options=[{'label': var, 'value': var} for var in numerical_data.columns],
                value=numerical_data.columns[9],  # Set default value
            ),
            html.Label("Choose variable for Y-axis:"),
            dcc.Dropdown(
                id='scatter-y-variable-dropdown',
                options=[{'label': var, 'value': var} for var in numerical_data.columns],
                value=numerical_data.columns[12],  # Set default value
            ),
            html.Div(style={'margin-top': '20px'}),
            html.Button("Show All Data", id="show-all-data-button"),

            html.Div(style={'display': 'inline'}, children=[  # Add this inline style to the button container
                html.Div(style={'margin': '20px'}),
                # html.Button("Show All Data", id="show-all-data-button", n_clicks=0),
            ]),
            html.Div(
        dcc.Loading(
            type="circle",
            children=[
                dcc.Graph(id='scatter-plot', selectedData=None),
            ],
        )
    ),
        ]),

    ], style={'margin-bottom': '20px'}),

    # Scatter Matrix section 
    html.Div(className='row', children=[
        html.Div(className='twelve columns', style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)'}, children=[
            html.H3('Scatter Matrix'),

            # Dropdown for scatter matrix
            html.Label("Choose variables for scatter matrix:"),
            dcc.Dropdown(
                id='scatter-matrix-variables',
                options=[{'label': var, 'value': var} for var in numerical_data.columns],
                value=[numerical_data.columns[11], numerical_data.columns[12]],  # Set default value
                multi=True,  # Allow multiple variable selection
            ),

            html.Div(
                dcc.Loading(
                    type="circle",
                    children=[
                        dcc.Graph(id='scatter-matrix-plot', selectedData=None),
                    ],
                )
    ),
        ]),
    ], style={'margin-bottom': '20px'}),

    # Second layout starts here
    html.Div(
        children=[
            html.Div(
                dcc.Loading(
                    type="circle",
                    children=[html.Div("Loading...", style={"font-size": "20px"})],
                    fullscreen=True,
                ),
                id="loading-spinner",
            ),
            html.Div(
                className="row",
                children=[
                    html.Div(
                        className="eight columns", style={'background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'},
                        children=[
                            html.Div(
                                children=[

                                    html.Div(
                                        dcc.Loading(
                                            type="circle",
                                            children=[
                                                dcc.Graph(
                                        id="cluster-graph",
                                        hoverData={"points": [{"hovertemplate": ""}]},
                                    ),
                                            ],
                                        )
                                    ),
                                ]
                            ),
                        ],
                    ),
                    html.Div(
                        className="four columns",
                        children=[
                            html.Div(
                                children=[
                                    html.Label("Method:"),
                                    dcc.Dropdown(
                                        id="method-dropdown",
                                        options=[
                                            {"label": "UMAP", "value": "umap"},
                                            {"label": "t-SNE", "value": "tsne"},
                                            {"label": "PCA", "value": "pca"},
                                            {"label": "Isomap", "value": "isomap"},
                                            {"label": "LLE", "value": "lle"},
                                            {
                                                "label": "Linear Layout Algorithm",
                                                "value": "linear_layout_algorithm",
                                            },
                                        ],
                                        value="umap",
                                        style={
                                            "width": "150px",
                                        },
                                    ),
                                ]
                            ),
                            html.Div(
                                children=[
                                    html.Label("Number of Components:"),
                                    dcc.Input(
                                        id="n-components-input",
                                        type="number",
                                        min=2,
                                        max=len(available_variables),
                                        value=3,
                                        step=1,
                                        style={
                                            "width": "100px",
                                            "opacity": 1,  # Full opacity by default
                                        },
                                    ),
                                ],
                                id="n-components-container",
                            ),
                            html.Div(
                                children=[
                                    html.Label("Number of Neighbours:"),
                                    dcc.Input(
                                        id="n-neighbours-input",
                                        type="number",
                                        min=2,
                                        max=len(available_variables),
                                        value=5,
                                        step=1,
                                        style={
                                            "width": "100px",
                                            "display": "inline-block",
                                            "opacity": 1,  # Full opacity by default
                                        },
                                    ),
                                ],
                                id="n-neighbours-container",
                            ),
                            html.Div(
                                children=[
                                    html.Label("Select Color Variable:"),
                                    dcc.Dropdown(
                                        id="color-variable-input",
                                        options=[
                                            {
                                                "label": variable,
                                                "value": variable,
                                            }
                                            for variable in available_variables[:32]
                                        ],
                                        value="chemotherapy",
                                        style={
                                            "width": "150px",
                                        },
                                    ),
                                ]
                            ),
                            html.Div(
                                className="row",
                                children=[
                                    html.Div(
                                        className="six columns",
                                        children=[
                                            html.Label(
                                                "Select Variables for Clustering Graph:"
                                            ),
                                            html.Div(
                                                className="scrollable",
                                                style={
                                                    "height": "230px",
                                                    "overflowY": "scroll",
                                                    "border": "1px solid #ccc",
                                                    "border-radius": "5px",
                                                    "padding": "5px",
                                                },
                                                children=[
                                                    dcc.Checklist(
                                                        id="variable-checkboxes",
                                                        options=[
                                                            {
                                                                "label": variable,
                                                                "value": variable,
                                                            }
                                                            for variable in numerical_data
                                                        ],
                                                        value=[
                                                            "brca1",
                                                            "brca2",
                                                            "palb2",
                                                            "pten",
                                                            "tp53",
                                                            "atm",
                                                            "cdh1",
                                                            "chek2",
                                                            "nbn",
                                                            "nf1",
                                                            "stk11",
                                                            "bard1",
                                                            "mlh1",
                                                            "msh2",
                                                            "msh6",
                                                            "pms2",
                                                            "epcam",
                                                            "rad51c",
                                                            "rad51d",
                                                            "rad50",
                                                            "rb1",
                                                            "rbl1",
                                                            "rbl2",
                                                        ],
                                                    ),
                                                ],
                                            ),
                                        ],
                                    ),
                                    html.Div(
                                        className="six columns",
                                        children=[
                                            html.Label(
                                                "Select Variables for Parallel Coordinates:"
                                            ),
                                            html.Div(
                                                className="scrollable",
                                                style={
                                                    "height": "230px",
                                                    "overflowY": "scroll",
                                                    "border": "1px solid #ccc",
                                                    "border-radius": "5px",
                                                    "padding": "5px",
                                                },
                                                children=[
                                                    dcc.Checklist(
                                                        id="show-variables-checkboxes",
                                                        options=[
                                                            {
                                                                "label": variable,
                                                                "value": variable,
                                                            }
                                                            for variable in numerical_data
                                                        ],
                                                        value=[
                                                            "brca1",
                                                            "brca2",
                                                            "palb2",
                                                            "pten",
                                                            "tp53",
                                                            "atm",
                                                            "cdh1",
                                                            "chek2",
                                                        ],
                                                    ),
                                                ],
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
            ),
            html.Div(
                className="row", style={'margin-top': '20px', 'background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'},
                children=[
                    html.Div(
                        className="twelve columns",
                        children=[
                            html.Div(
                                children=[
                                    html.Label("Parallel Coordinates graph:"),
                                    html.Div(
                                        dcc.Loading(
                                            type="circle",
                                            children=[
                                                dcc.Graph(id="parallel-coordinates-graph"),
                                            ],
                                        )
                                    ),                                    
                                ]
                            ),
                        ],
                    ),
                ],
            ),
            html.Div(
                className="row", style={'margin-top': '20px','background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'},
                children=[
                    html.Div(
                        style={"width": "100%"},
                        children=[
                            html.H4("Data Table"),
                            dash_table.DataTable(
                                id="data-table",
                                columns=[
                                    {"name": col, "id": col} for col in data.columns
                                ],
                                data=data.to_dict("records"),
                                style_table={
                                    "height": "300px",
                                    "overflowY": "scroll",
                                },
                                style_cell={"padding": "3px"},
                                style_header={
                                    "backgroundColor": "rgb(230, 230, 230)",
                                    "fontWeight": "bold",
                                },
                                row_selectable="multi",
                                selected_rows=[],
                                active_cell={
                                    "row": 1,
                                    "column": 0,
                                }, 
                            ), 
                        ], 
                    )
                ],
            ),
            html.Div(className="row", style={'margin-top': '20px'}, children=[
        
    ]),

    # Second row: Patient information table, spider chart, and mutation count pie chart
    html.Div(className="row", children=[
        html.Div(className="six columns", style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)', 'margin-right': '10px'}, children=[
            html.H2(children='Patient Information:'),
            html.Table(id='patient-info-table',
                       children=[
                           html.Tr([html.Th("Variable Name"), html.Th("Value")])
                       ]),
        ]),
        html.Div(className="six columns", style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)', 'margin-left': '10px'}, children=[
            dcc.Graph(
                id='spider-chart',
                figure={'layout': {'width': 370, 'height': 285}},  
            ),
            dcc.Graph(
                id='mutation-pie-chart',
                figure={'layout': {'width': 370, 'height': 285}}, 
            ),
            
        ]),
        
    ]),

    # Third row: Heatmap and barchart
    html.Div(className="row", children=[

        html.Div( style={'padding': '10px', 'border': '1px solid #ccc', 'border-radius': '5px', 'box-shadow': '2px 2px 5px rgba(0, 0, 0, 0.1)', 'margin-left': '10px'}, children=[
            dcc.Graph(id='gene-zscore-bar-chart'),
        ]),
    ]),
        ],
        style={
            "font-family": "Arial, sans-serif",
            "margin": "2px",
            "padding": "2px",
        },
    ),

], className='app-container')


# Define a callback to update the correlation heatmap plot
@callback(
    dash.dependencies.Output('correlation-heatmap', 'figure'),
    [dash.dependencies.Input('statistics-variable-dropdown', 'value')]
)
def update_correlation_heatmap(selected_variable):
    # Calculate the correlation matrix between all numerical variables
    correlation_matrix = numerical_data.iloc[:,:503].corr()

    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=correlation_matrix.values,
        x=correlation_matrix.columns,
        y=correlation_matrix.columns,
        colorscale='RdBu',  
        zmin=-1,
        zmax=1,
        colorbar=dict(title='Correlation'),
    ))

    layout = go.Layout(
        title='Correlation Heatmap',
        xaxis=dict(tickangle=45),
        yaxis=dict(tickangle=0),
        height=700,
    )
    fig.update_layout(layout)

    return fig

# Define a callback to update the selected variable based on click data
@callback(
    Output('variable-dropdown', 'value'),
    [Input('correlation-heatmap', 'clickData')],
    [State('variable-dropdown', 'value')]
)
def update_selected_variable_from_click_data(click_data, selected_variable):
    if click_data is not None:
        selected_variable = click_data['points'][0]['x']
    return selected_variable

# Define a callback to update the correlation table when a new variable is selected in the dropdown
@callback(
    Output('correlation-table-body', 'children'),
    [Input('variable-dropdown', 'value')]
)
def update_correlation_table(selected_variable):
    if selected_variable is not None:
        # Calculate pairwise correlations
        correlations = numerical_data.corr()[selected_variable].reset_index()
        correlations.columns = ['Variable', 'Correlation']

        # Sort correlations by absolute values in descending order
        correlations = correlations.iloc[correlations['Correlation'].abs().sort_values(ascending=False).index]

        # Generate the table rows with correlation bars
        table_rows = [
            html.Tr([
                html.Td(correlations.iloc[i]['Variable']),
                html.Td(f"{correlations.iloc[i]['Correlation']:.2f}"),
                html.Td(html.Div(style={
                    'width': f"{abs(correlations.iloc[i]['Correlation']) * 50}px",
                    'height': '20px',
                    'background-color': 'orange' if correlations.iloc[i]['Correlation'] > 0 else 'blue',
                    'margin-right': '10px',  # Add some space between the correlation value and the bar
                    'position': 'relative',  # Set position to relative for proper alignment
                })),
            ])
            for i in range(len(correlations))
        ]

        return table_rows

    # If no variable is selected, return an empty table
    return []

# Define a dictionary of variable descriptions
variable_descriptions = {
    'patient_id': 'Unique identifier for each patient',
    'age_at_diagnosis': 'Age of the patient at diagnosis time',
    'type_of_breast_surgery': 'Breast cancer surgery type: 1- MASTECTOMY, which refers to a surgery to remove all breast tissue from a breast as a way to treat or prevent breast cancer. 2- BREAST CONSERVING, which refers to a urgery where only the part of the breast that has cancer is removed',
    'cancer_type': 'Breast cancer types: 1- Breast Cancer or 2- Breast Sarcoma', 
    'cancer_type_detailed': 'Detailed Breast cancer types: 1- Breast Invasive Ductal Carcinoma 2- Breast Mixed Ductal and Lobular Carcinoma 3- Breast Invasive Lobular Carcinoma 4- Breast Invasive Mixed Mucinous Carcinoma 5- Metaplastic Breast Cancer', 
    'cellularity': 'Cancer cellularity post chemotherapy, which refers to the amount of tumor cells in the specimen and their arrangement into clusters', 
    'chemotherapy': 'Whether or not the patient had chemotherapy as a treatment (yes/no)', 
    'pam50_+_claudin-low_subtype': 'Pam 50: is a tumor profiling test that helps show whether some estrogen receptor-positive (ER-positive), HER2-negative breast cancers are likely to metastasize (when breast cancer spreads to other organs). The claudin-low breast cancer subtype is defined by gene expression characteristics, most prominently: Low expression of cell–cell adhesion genes, high expression of epithelial–mesenchymal transition (EMT) genes, and stem cell-like/less differentiated gene expression patterns', 
    'cohort': 'Cohort is a group of subjects who share a defining characteristic (It takes a value from 1 to 5)', 
    'er_status_measured_by_ihc': 'To assess if estrogen receptors are expressed on cancer cells by using immune-histochemistry (a dye used in pathology that targets specific antigen, if it is there, it will give a color, it is not there, the tissue on the slide will be colored) (positive/negative)', 
    'er_status': 'Cancer cells are positive or negative for estrogen receptors', 
    'neoplasm_histologic_grade': 'Determined by pathology by looking the nature of the cells, do they look aggressive or not (It takes a value from 1 to 3)', 
    'her2_status_measured_by_snp6': 'To assess if the cancer positive for HER2 or not by using advance molecular techniques (Type of next generation sequencing)', 
    'her2_status': 'Whether the cancer is positive or negative for HER2', 
    'tumor_other_histologic_subtype': "Type of the cancer based on microscopic examination of the cancer tissue (It takes a value of 'Ductal/NST', 'Mixed', 'Lobular', 'Tubular/ cribriform', 'Mucinous', 'Medullary', 'Other', 'Metaplastic' )", 
    'hormone_therapy': 'Whether or not the patient had hormonal as a treatment (yes/no)', 
    'inferred_menopausal_state': 'Whether the patient is is post menopausal or not (post/pre)', 
    'integrative_cluster': "Molecular subtype of the cancer based on some gene expression (It takes a value from '4ER+', '3', '9', '7', '4ER-', '5', '8', '10', '1', '2', '6')", 
    'primary_tumor_laterality': 'Whether it is involving the right breast or the left breast', 
    'lymph_nodes_examined_positive': 'To take samples of the lymph node during the surgery and see if there were involved by the cancer', 
    'mutation_count': 'Number of gene that has relevant mutations', 
    'nottingham_prognostic_index': 'It is used to determine prognosis following surgery for breast cancer. Its value is calculated using three pathological criteria: the size of the tumour; the number of involved lymph nodes; and the grade of the tumour.', 
    'oncotree_code': 'The OncoTree is an open-source ontology that was developed at Memorial Sloan Kettering Cancer Center (MSK) for standardizing cancer type diagnosis from a clinical perspective by assigning each diagnosis a unique OncoTree code.', 
    'overall_survival_months': 'Duration from the time of the intervention to death', 
    'overall_survival': 'Target variable wether the patient is alive of dead.', 
    'pr_status': 'Cancer cells are positive or negative for progesterone receptors', 
    'radio_therapy': 'Whether or not the patient had radio as a treatment (yes/no)', 
    '3-gene_classifier_subtype': "Three Gene classifier subtype It takes a value from 'ER-/HER2-', 'ER+/HER2- High Prolif', nan, 'ER+/HER2- Low Prolif','HER2+'", 
    'tumor_size': 'Tumor size measured by imaging techniques', 
    'tumor_stage': 'Stage of the cancer based on the involvement of surrounding structures, lymph nodes and distant spread', 
    'death_from_cancer': "Wether the patient's death was due to cancer or not", 
}

# Define a callback to update the statistics table
@callback(
    dash.dependencies.Output('statistics-table-body', 'children'),
    [dash.dependencies.Input('statistics-variable-dropdown', 'value')]
)
def update_statistics_table(selected_variable):
    # Calculate statistics for all variables
    statistics = []
    for column in data.columns:
        num_values = data[column].count()
        num_gaps = data[column].isnull().sum()
        num_unique_values = data[column].nunique()
        most_common_value = data[column].mode().values[0]

        if pd.api.types.is_numeric_dtype(data[column]):
            # For numerical variables, calculate additional statistics
            num_zeros = (data[column] == 0).sum()
            min_value = data[column].min()
            max_value = data[column].max()
            median = round(data[column].median(), 2)
            mean = round(data[column].mean(), 2)
            std_value = round(data[column].std(), 2)
            num_outliers_2sigma = data[(data[column] - mean).abs() > 2 * std_value].shape[0]
            num_outliers_3sigma = data[(data[column] - mean).abs() > 3 * std_value].shape[0]
            num_outliers_4sigma = data[(data[column] - mean).abs() > 4 * std_value].shape[0]
            num_text_values = '-'
        else:
            num_zeros = '-'
            min_value = '-'
            max_value = '-'
            median = '-'
            mean = '-'
            std_value = '-'
            num_outliers_2sigma = '-'
            num_outliers_3sigma = '-'
            num_outliers_4sigma = '-'

        statistics.append([
            column, num_values, num_gaps, num_unique_values, num_zeros,
            most_common_value, min_value, max_value, median, mean, std_value,
            num_outliers_2sigma, num_outliers_3sigma, num_outliers_4sigma,
        ])

    # Filter rows based on the selected variable
    if selected_variable:
        statistics = [row for row in statistics if row[0] == selected_variable]

    table_rows = [
        html.Tr(
            title=variable_descriptions.get(stat_row[0], ''),  # Get the description from the dictionary
            children=[html.Td(stat) for stat in stat_row]
        )
        for stat_row in statistics
    ]

    return table_rows


# Define a callback to update the frequency graph
@callback(
    dash.dependencies.Output('frequency-graph', 'figure'),
    [dash.dependencies.Input('frequency-variable-dropdown', 'value')]
)
def update_frequency_graph(selected_variable):
    # Calculate the frequency distribution by label
    frequency_data = data.groupby(selected_variable).size().reset_index(name='Frequency')

    # Get a list of unique labels for the selected variable
    labels = frequency_data[selected_variable]

    # Create a color scale for the bars based on the number of labels
    colors = [f"hsl({hue}, 50%, 50%)" for hue in range(0, 360, int(360 / len(labels)))]

    # Create the bar chart
    trace = go.Bar(
        x=frequency_data[selected_variable],
        y=frequency_data['Frequency'],
        marker=dict(color=colors)
    )

    layout = go.Layout(
        title=f'Frequency Distribution of {selected_variable}',
        xaxis=dict(title=selected_variable),
        yaxis=dict(title='Frequency'),
    )

    figure = go.Figure(data=[trace], layout=layout)

    return figure


def find_matching_variable(data, selected_value):
    for column in data.columns:
        if selected_value in data[column].values:
            return column
    return None

@callback(Output("loading-spinner", "style"), Input("data-table", "data"))
def hide_loading_spinner(data):
    if data is not None:
        return {"display": "none"}
    else:
        return {"display": "block"}

# Callback to adjust the opacity/transparency based on the selected method
@callback(
    [
        dash.dependencies.Output("n-components-container", "style"),
        dash.dependencies.Output("n-neighbours-container", "style"),
    ],
    [dash.dependencies.Input("method-dropdown", "value")]
)
def adjust_opacity_based_on_method(method):
    if method in ["umap", "isomap", "lle"]:
        return {"opacity": 1}, {"opacity": 1}
    elif method == "tsne":
        return {"opacity": 1}, {"opacity": 0.5}
    else:  # For PCA and Linear Layout Algorithm
        return {"opacity": 0.5}, {"opacity": 0.5}


# Function to generate Spider Chart
def generate_spider_chart(patient_id):
    # Filter data for the selected patient ID
    selected_patient_data = data[data['patient_id'] == patient_id]

    # Select the variables for the spider chart
    spider_variables = ['neoplasm_histologic_grade',
                        'overall_survival_months', 'integrative_cluster',
                        'nottingham_prognostic_index']

    # Create the spider chart traces
    spider_chart_trace = go.Scatterpolar(
        r=selected_patient_data[spider_variables].values.tolist()[0],
        theta=spider_variables,
        fill='toself',
        name='Patient ' + str(patient_id)
    )

    return {
        'data': [spider_chart_trace],
        'layout': go.Layout(
            title=f'Spider Chart for Patient {patient_id}',
            polar=dict(
                radialaxis=dict(
                    visible=True,
                ),
            ),
        )
    }

def generate_patient_info_table(patient_id):
    # Filter data for the selected patient ID
    selected_patient_data = data[data['patient_id'] == patient_id]

    # Select the variables for the patient information table
    patient_info_variables = ['age_at_diagnosis', 'cancer_type_detailed',
                              'type_of_breast_surgery', 'chemotherapy',
                              'hormone_therapy', 'inferred_menopausal_state',
                              'lymph_nodes_examined_positive', 'tumor_size',
                              'tumor_stage']

    # Create the patient information table as a list of HTML table rows
    patient_info_rows = []
    for col in patient_info_variables:
        variable_name = col.replace("_", " ").title()
        variable_value = selected_patient_data[col].values[0]

        # If the column is Chemotherapy or Hormone Therapy, replace the value
        # 1 with "Yes" and 0 with "No"
        if col == 'chemotherapy' or col == 'hormone_therapy':
            variable_value = "Yes" if variable_value == 1 else "No"

        patient_info_rows.append(html.Tr([html.Td(variable_name), html.Td(variable_value)]))

    # Add Overall Survival Status row
    overall_survival_status = selected_patient_data['death_from_cancer'].values[0]
    patient_info_rows.append(html.Tr([html.Td("Overall Survival Status"), html.Td(overall_survival_status)]))

    return patient_info_rows


# Function to generate the pie chart of mutation_count
def generate_mutation_pie_chart(patient_id):
    # Filter data for the selected patient ID
    selected_patient_data = data[data['patient_id'] == patient_id]

    # Get mutation_count value
    mutation_count = selected_patient_data['mutation_count'].values[0]

    # If mutation_count is null or 0 or empty, display "-" in the pie chart
    if pd.isnull(mutation_count) or mutation_count == 0 or mutation_count == "":
        mutation_count = '-'
        chart_title = f'No data available about mutation count for Patient {patient_id}'
    else:
        chart_title = f'Mutation Count Pie Chart for Patient {patient_id}'

    # Create the pie chart trace
    pie_chart_trace = go.Pie(
        values=[mutation_count],
        labels=['Mutation Count'],
        hole=0.6,
        textinfo='label+value',  # Show label (Mutation Count) and value
        hoverinfo='label+value',  # Show label (Mutation Count) and value on hover
        text=[str(mutation_count)],
        showlegend=False, 
    )

    return {
        'data': [pie_chart_trace],
        'layout': go.Layout(
            title=chart_title,
        )
    }


# Updated function to generate the gene z-score bar chart
def generate_gene_zscore_bar_chart(patient_id):
    # Filter data for the selected patient ID
    selected_patient_data = data[data['patient_id'] == patient_id]

    # Select the gene z-score values for columns 520 to the end
    gene_mut_zscore_data = selected_patient_data.iloc[:, 520:]

    columns_with_non_zero_values = []

    # Loop through the columns and check for non-zero values
    for col in gene_mut_zscore_data.columns:
        if selected_patient_data[col].values[0] != 0 and selected_patient_data[col].values[0] != '0':
            columns_with_non_zero_values.append(selected_patient_data[col].name.replace('_mut', ''))

    # Select the gene z-score values for columns 31 to 520
    gene_zscore_data = selected_patient_data.iloc[:, 31:326].values.tolist()[0]

    # Create the bar chart trace
    bar_chart_trace = go.Bar(
        x=data.columns[31:326],  # Gene names
        y=gene_zscore_data,
        marker=dict(
            color=['red' if gene_name in columns_with_non_zero_values else 'blue' for gene_name in data.columns[31:326]],
        ),
    )


    # Create the figure for the bar chart
    bar_chart_figure = {
        'data': [bar_chart_trace],
        'layout': go.Layout(
            title=f'Gene Z-Score Mutation Bar Chart for Patient {patient_id}',
            xaxis=dict(title='Genes'),
            yaxis=dict(title='Z-Score'),
        )
    }

    return bar_chart_figure

# Callback to update all elements when the patient ID is changed
@callback(
    [Output('patient-info-table', 'children'),
     Output('spider-chart', 'figure'),
     Output('mutation-pie-chart', 'figure'),
     Output('gene-zscore-bar-chart', 'figure')],
     Input("data-table", "active_cell"),
     Input("data-table", "data")
)
def update_all_elements(active_cell, data):
    if active_cell is None:
        patient_id = 0
    else:
        row = active_cell["row"]
        col = active_cell["column"]
        patient_id = data[row]["patient_id"]
    # Generate patient information table
    patient_info_table_rows = generate_patient_info_table(patient_id)

    # Generate spider chart
    spider_chart_figure = generate_spider_chart(patient_id)

    # Generate gene z-score bar chart
    gene_zscore_bar_chart_figure = generate_gene_zscore_bar_chart(patient_id)

    return patient_info_table_rows, spider_chart_figure, gene_zscore_bar_chart_figure


@callback(
    Output("parallel-coordinates-graph", "figure"),
    [Input("show-variables-checkboxes", "value"),
     Input("data-table", "selected_rows"),
     Input('data-table', 'data')],
)
def update_parallel_coordinates(selected_show_variables, selected_rows, data_table):
    
    datatable = pd.DataFrame(data_table)

    df_selected_variables = datatable[selected_show_variables]

    # Create a copy of the data for plotting
    plot_data = df_selected_variables.copy()

    # Initialize the color column for all rows to 0 (unselected)
    plot_data["color"] = 0

    if selected_rows:
        valid_selected_indices = [row for row in selected_rows if row < len(datatable)]
        selected_indices = [datatable.iloc[row].name for row in valid_selected_indices]
        plot_data.loc[selected_indices, "color"] = 1

    fig = go.Figure(data=go.Parcoords(
        line=dict(color=plot_data["color"],  # Color by the 'color' column
                  colorscale=[[0, '#EAEDED'], [1, 'red']],  # Colorscale: gray for unselected, red for selected
                  showscale=False),  # Show color scale bar
        dimensions=[
            dict(range=[min(plot_data[col]), max(plot_data[col])],
                 label=col,
                 values=plot_data[col])
            for col in selected_show_variables]
    ))

    if selected_rows:
        valid_selected_indices = [row for row in selected_rows if row < len(datatable)]
        selected_indices = [datatable.iloc[row].name for row in valid_selected_indices]
        
        # Ensure that the selected indices are in the datatable
        valid_selected_indices = [
            idx for idx in selected_indices if idx in datatable.index
        ]
        
        fig.update_traces(
            line=dict(color="red"), selector={"line.color": "gray"}  # Highlight selected rows in red
        )


    return fig


@callback(
    Output('data-table', 'data'),
    Input('frequency-graph', 'clickData'),
    State('frequency-variable-dropdown', 'value'),
    Input("show-all-data-button", "n_clicks")
)
def update_data_table(clickData, var, n_clicks):
    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'] if ctx.triggered else None
    
    if triggered_id == "show-all-data-button.n_clicks":
        return data.to_dict("records")
    elif clickData:
        clicked_value = clickData["points"][0]["x"]
        return data[data[var] == clicked_value].to_dict("records")
    else:
        return data.to_dict("records")

@callback(
    [Output('scatter-plot', 'figure'),
    Output('scatter-matrix-plot', 'figure'),
     Output('cluster-graph', 'figure'),
     Output("data-table", "selected_rows")],
    [Input('scatter-x-variable-dropdown', 'value'),
     Input('scatter-y-variable-dropdown', 'value'),
     Input('scatter-matrix-variables', 'value'),
     Input('scatter-plot', 'clickData'),
     Input('scatter-matrix-plot', 'clickData'),
     Input('frequency-graph', 'clickData'),
     Input('data-table', 'selected_rows'),
     Input('data-table', 'data'),
     Input("variable-checkboxes", "value"),
     Input("method-dropdown", "value"),
     Input("n-components-input", "value"),
     Input("n-neighbours-input", "value"),
     Input("color-variable-input", "value"),
     Input("cluster-graph", "clickData")],
    [State("cluster-graph", "figure"),
     State("data-table", "selected_rows")]
)
def update_plots(x_variable, y_variable, selected_variables, scatter_click_data, matrix_click_data, frequency_click_data,
                 selected_rows_data_table, data_table_data, 
                 selected_variables_cluster, method, n_components, n_neighbours, color_variable, click_data_cluster_graph,
                current_figure_cluster_graph, current_selected_rows_cluster_graph):

    # Initialize variables
    cache_key = None

    # Determine which component triggered the callback
    triggered_component = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

    # Handle scatter matrix plot updates
    if triggered_component == "scatter-matrix-plot":
        if matrix_click_data and matrix_click_data.get("points"):
            clicked_point = matrix_click_data["points"][0]
            if "pointIndex" in clicked_point:
                point_index = clicked_point["pointIndex"]
                if point_index not in current_selected_rows_cluster_graph:
                    current_selected_rows_cluster_graph.append(point_index)
                else:
                    current_selected_rows_cluster_graph.remove(point_index)
    
    if triggered_component == "scatter-plot":
        if scatter_click_data and scatter_click_data.get("points"):
            clicked_point = scatter_click_data["points"][0]
            if "pointIndex" in clicked_point:
                point_index = clicked_point["pointIndex"]
                if point_index not in current_selected_rows_cluster_graph:
                    current_selected_rows_cluster_graph.append(point_index)
                else:
                    current_selected_rows_cluster_graph.remove(point_index)

    # Handle data table updates
    if triggered_component == "data-table":
        # Convert selected_rows_data_table to a set for efficient handling
        selected_rows_set = set(selected_rows_data_table)
        # Update the scatter trace to highlight selected points in the cluster graph
        if current_figure_cluster_graph is not None:
            current_figure_cluster_graph["data"][0]["selectedpoints"] = [
                data[data["patient_id"].eq(data_table_data[row]["patient_id"])].index[0]
                for row in selected_rows_set
            ]
            current_figure_cluster_graph["data"][0]["selected"] = dict(marker=dict(color="red", size=10))
        # Update current_selected_rows_cluster_graph with the latest selections from the data table
        current_selected_rows_cluster_graph = list(selected_rows_set)

    # Handle cluster graph updates
    if triggered_component == "cluster-graph":
        if click_data_cluster_graph and click_data_cluster_graph.get("points"):
            clicked_point = click_data_cluster_graph["points"][0]
            if "pointIndex" in clicked_point:
                point_index = clicked_point["pointIndex"]
                if point_index not in current_selected_rows_cluster_graph:
                    current_selected_rows_cluster_graph.append(point_index)
                else:
                    current_selected_rows_cluster_graph.remove(point_index)

    datatable = pd.DataFrame(data_table_data)
    # Determine whether to show all data or filtered data

    # Create the scatter matrix plot
    num_variables = len(selected_variables)
    scatter_matrix_fig = make_subplots(rows=num_variables, cols=num_variables, shared_xaxes=True, shared_yaxes=True)

    # Loop through each pair of selected variables and create the scatter plot or bar plot accordingly
    for i in range(num_variables):
        for j in range(num_variables):
            x_var = selected_variables[i]
            y_var = selected_variables[j]

            if i == j:
                # Create bar plot for diagonal elements
                frequency_data = datatable[x_var].value_counts().reset_index()
                frequency_data.columns = [x_var, 'Frequency']
                trace = go.Bar(
                    x=frequency_data[x_var],
                    y=frequency_data['Frequency'],
                    marker=dict(
                        color= 'blue', 
                    ),
                    text=frequency_data[x_var],  # Set patient IDs as hover text
                    hovertemplate="Patient ID: %{text}<br>"  # Define the hover template
                                  f"{x_var}: %{{x}}<br>"
                                  f"Frequency: %{{y}}<extra></extra>",  # Include the variable values in hover
                )
                scatter_matrix_fig.add_trace(trace, row=i + 1, col=j + 1)
            else:
                # Create scatter plot for off-diagonal elements
                trace = go.Scatter(
                    x=datatable[x_var],
                    y=datatable[y_var],
                    mode='markers',
                    marker=dict(
                        size=8,
                        color= 'blue',
                        opacity=0.7,
                        line=dict(width=0)
                    ),
                    text=datatable['patient_id'],  # Set patient IDs as hover text
                    hovertemplate="Patient ID: %{text}<br>"  # Define the hover template
                                  f"{x_var}: %{{x}}<br>"
                                  f"{y_var}: %{{y}}<extra></extra>",  # Include the variable values in hover
                )
                scatter_matrix_fig.add_trace(trace, row=i + 1, col=j + 1)

    # Update layout of the scatter matrix plot
    layout = go.Layout(
        title='Scatter Matrix',
        showlegend=False,
        height=800,
        width=1400,
    )
    scatter_matrix_fig.update_layout(layout)

    if current_selected_rows_cluster_graph:
        valid_selected_indices = [
            row for row in current_selected_rows_cluster_graph if row < len(data_table_data)
        ]
        
        for trace in scatter_matrix_fig.data:
            if isinstance(trace, go.Scatter):
                selected_indices = [
                    datatable[datatable["patient_id"].eq(data_table_data[row]["patient_id"])].index[0]
                    for row in valid_selected_indices if row < len(data_table_data)
                ]
                
                # Ensure that the selected indices are in the datatable.index
                valid_selected_indices = [
                    idx for idx in selected_indices if idx in datatable.index
                ]
                
                trace.selectedpoints = valid_selected_indices
                trace.selected = dict(marker=dict(color="red", size=10))
            elif isinstance(trace, go.Bar):
                trace.selected = dict(marker=dict(color="red"))


    trace = go.Scatter(
        x=datatable[x_variable],
        y=datatable[y_variable],
        mode='markers',
        marker=dict(
            size=8,
            color='blue',
            opacity=0.7,
            line=dict(width=0)
        ),
        text=datatable['patient_id'],  # Set patient IDs as hover text
        hovertemplate="Patient ID: %{text}<br>"  # Define the hover template
                     f"{x_variable}: %{{x}}<br>"
                     f"{y_variable}: %{{y}}<extra></extra>",  # Include the variable values in hover
    )

    layout = go.Layout(
        title='Scatter Plot',
        xaxis=dict(title=x_variable),
        yaxis=dict(title=y_variable),
    )

    figure = go.Figure(data=[trace], layout=layout)

    if current_selected_rows_cluster_graph:
        valid_selected_indices = [
            row for row in current_selected_rows_cluster_graph if row < len(data_table_data)
        ]
        
        selected_indices = [
            datatable[datatable["patient_id"].eq(data_table_data[row]["patient_id"])].index[0]
            for row in valid_selected_indices if row < len(data_table_data)
        ]
        
        # Ensure that the selected indices are in the datatable.index
        valid_selected_indices = [
            idx for idx in selected_indices if idx in datatable.index
        ]
        
        if valid_selected_indices:
            figure["data"][0]["selectedpoints"] = valid_selected_indices
            figure["data"][0]["selected"] = dict(marker=dict(color="red", size=10))


    # Perform dimensionality reduction based on the selected method
    if data_table_data is not None:
        data_subset = pd.DataFrame(data_table_data)
        selected_data = data_subset[selected_variables_cluster]
    else:
        selected_data = datatable[selected_variables_cluster]  # Use filtered_data instead of full data


    if method == "umap":
        reducer = umap.UMAP(
            n_neighbors=n_neighbours, min_dist=0.3, n_components=n_components, random_state=42
        )
    elif method == "tsne":
        reducer = TSNE(n_components=n_components, random_state=42)
    elif method == "pca":
        reducer =  PCA(
            n_components=min(n_components, len(selected_variables_cluster)),
            svd_solver="randomized",
            random_state=42,
        )
    elif method == "isomap":
        reducer =  Isomap(n_neighbors=n_neighbours, n_components=n_components)
    elif method == "lle":
        reducer =  LocallyLinearEmbedding(n_neighbors=n_neighbours, n_components=n_components)
    elif method == "linear_layout_algorithm":
        reducer =  None 

    if reducer:
        embedding = reducer.fit_transform(selected_data)
        # Plot the scatter plot with different colors for each variable
        cluster_fig = px.scatter(
            data_frame=datatable,
            x=embedding[:, 0],
            y=embedding[:, 1],
            color=color_variable,
            title="Clustering",
            hover_data={"patient_id": True},
        )
    else:
        # Apply linear layout algorithm
        positions = linear_layout_algorithm(selected_data.values)
        df = pd.DataFrame(selected_data, columns=selected_variables_cluster)
        df["x"] = positions[:, 0]
        df["y"] = positions[:, 1]
        df['color_variable'] = datatable[color_variable]  # Assign the color variable to the DataFrame
        df['patient_id'] = datatable['patient_id']  # Include the 'patient_id' column in df

        # Plot the scatter plot with different colors for each variable
        cluster_fig = go.Figure()
        cluster_fig.add_trace(
            go.Scatter(
                x=df["x"],
                y=df["y"],
                mode='markers',
                marker=dict(
                    size=8,
                    color=df['color_variable'],  # Use the assigned color variable
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title=color_variable),
                ),
                text=df['patient_id'],  # Display patient ID when hovering
                hovertemplate="<b>Patient ID:</b> %{text}<extra></extra>",
            )
        )
        cluster_fig.update_traces(marker=dict(color=df['color_variable']))  # Update marker color based on selected variable

    cluster_fig.update_layout(transition_duration=200, width=850, height=600)

    if current_selected_rows_cluster_graph:
        valid_selected_indices = [
            row for row in current_selected_rows_cluster_graph if row < len(data_table_data)
        ]
        
        selected_indices = [
            datatable[datatable["patient_id"].eq(data_table_data[row]["patient_id"])].index[0]
            for row in valid_selected_indices if row < len(data_table_data)
        ]
        
        # Ensure that the selected indices are in the datatable.index
        valid_selected_indices = [
            idx for idx in selected_indices if idx in datatable.index
        ]
        
        if valid_selected_indices:
            cluster_fig["data"][0]["selectedpoints"] = valid_selected_indices
            cluster_fig["data"][0]["selected"] = dict(marker=dict(color="red", size=10))


    # Cache the result
    cache_key = (
        selected_variables,
        scatter_click_data,
        matrix_click_data,
        frequency_click_data,
        selected_rows_data_table,
        data_table_data,
        selected_variables_cluster,
        method,
        n_components,
        n_neighbours,
        color_variable,
        click_data_cluster_graph,
        current_figure_cluster_graph,
        current_selected_rows_cluster_graph,
    )

    # Store the result in cache
    cache[cache_key] = (scatter_matrix_fig, cluster_fig, current_selected_rows_cluster_graph)

    return figure, scatter_matrix_fig, cluster_fig, current_selected_rows_cluster_graph
