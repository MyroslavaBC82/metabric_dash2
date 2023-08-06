import dash
from dash import html, dcc, callback, Input, Output, State
import pandas as pd
import plotly.express as px
from sklearn.preprocessing import StandardScaler
import plotly.graph_objs as go
import numpy as np

dash.register_page(__name__)
data = pd.read_csv("METABRIC_RNA_Mutation3.csv", nrows=250)
available_variables = list(data.columns)

# Boxplot for numerical variables from columns 2 to 31
numerical_variables = data.iloc[:, 1:31].select_dtypes(include='number')
numerical_variables2 = data.select_dtypes(include='number')

# Standardize the data
scaler = StandardScaler()
numerical_variables_standardized = pd.DataFrame(scaler.fit_transform(numerical_variables), columns=numerical_variables.columns)

# Melt the DataFrame to reshape it properly for the boxplot
melted_data = numerical_variables_standardized.melt(var_name='variable', value_name='value')

boxplot_fig = px.box(melted_data,
                     y='variable',
                     x='value',
                     title='Boxplot of Clinical Attributes Distribution(Standardized)',
                     orientation='h',  # Set orientation to horizontal (vertical boxplot)
                     )  

def get_filtered_data(selected_data_source):
    if selected_data_source == 'overall_survival':
        return data['overall_survival'].unique()
    elif selected_data_source == 'death_from_cancer':
        return data['death_from_cancer'].unique()
    

# Modify the update_distplot callback to take an additional input clickData and State for dropdown value
@callback(
    Output('distplot', 'figure'),
    [Input('variable-dropdown', 'value'),
     Input('data-source-radio', 'value'),
     Input('boxplot', 'clickData'),
     Input('gene-mutation-distplot', 'clickData')],
    [State('variable-dropdown', 'value')]  # Add State for dropdown value
)
def update_distplot(selected_variable, selected_data_source, box_clickData, gene_clickData, dropdown_value):
    # Check if a box in the boxplot is clicked and get the variable from the clicked box
    if box_clickData and 'points' in box_clickData:
        clicked_variable = box_clickData['points'][0]['y']
    elif gene_clickData and 'points' in gene_clickData:  # Check if a bar in the Gene Mutation Distribution barchart is clicked
        clicked_variable = gene_clickData['points'][0]['x']
    else:
        clicked_variable = None

    if selected_variable:
        dropdown_value = selected_variable
    # Use the clicked variable if available, otherwise use the selected_variable
    variable_to_plot = clicked_variable or dropdown_value or selected_variable

    fig = go.Figure()

    color_map = {0: 'red', 1: 'green', 'Living': 'green', 'Died of Disease': 'yellow', 'Died of Other Causes': 'red'}  # Define a color map for 'Death from Cancer' statuses

    for status in get_filtered_data(selected_data_source):
        filtered_data = numerical_variables2[data[selected_data_source] == status][variable_to_plot]
        # Check if the status is in the color_map, if not, use the default color 'blue'
        marker_color = color_map.get(status, 'blue')
        if selected_data_source == 'overall_survival':
            # Replace legend labels for 'Overall Survival'
            status_label = 'Survived' if status == 1 else 'Died'
            fig.add_trace(go.Histogram(x=filtered_data.reset_index(drop=True),
                                       opacity=0.7,
                                       name=f'{selected_data_source}: {status_label}',
                                       marker_color=marker_color,
                                       histnorm='density',  # Display count on the y-axis
                                       ))
        else:
            # For other data sources, use the default behavior
            fig.add_trace(go.Histogram(x=filtered_data.reset_index(drop=True),
                                       opacity=0.7,
                                       name=f'{selected_data_source}: {status}',
                                       marker_color=marker_color,
                                       histnorm='density',  # Display count on the y-axis
                                       ))

    fig.update_layout(title_text=f'Distribution of {variable_to_plot} by {selected_data_source} (Standardized)',
                      barmode='overlay',
                      xaxis=dict(title=f'{variable_to_plot} (Standardized)'),
                      yaxis=dict(title='Density'))

    fig.update_traces(opacity=0.7)

    if variable_to_plot:
        selected_variable = variable_to_plot

    return fig


def get_color_map(data_source):
    if data_source == 'overall_survival':
        return {0: 'red', 1: 'green'}
    elif data_source == 'death_from_cancer':
        return {1: 'green', 2: 'yellow', 3: 'red'}

# Modify the update_boxplot callback to take an additional input clickData and State for dropdown value
@callback(
    Output('boxplot-survival', 'figure'),
    [Input('variable-dropdown', 'value'),
     Input('data-source-radio', 'value'),
     Input('boxplot', 'clickData'),
     Input('gene-mutation-distplot', 'clickData')],
    [State('variable-dropdown', 'value')]  # Add State for dropdown value
)
def update_boxplot(selected_variable, selected_data_source, box_clickData, gene_clickData, dropdown_value):
    # Check if a box in the boxplot is clicked and get the variable from the clicked box
    if box_clickData and 'points' in box_clickData:
        clicked_variable = box_clickData['points'][0]['y']
    elif gene_clickData and 'points' in gene_clickData:  # Check if a bar in the Gene Mutation Distribution barchart is clicked
        clicked_variable = gene_clickData['points'][0]['x']
    else:
        clicked_variable = None
    if selected_variable:
        dropdown_value = selected_variable
    # Use the clicked variable if available, otherwise use the selected_variable
    variable_to_plot = clicked_variable or dropdown_value or selected_variable

    fig = px.box(data,
                 x=selected_data_source,
                 y=variable_to_plot,
                 color=selected_data_source,
                 title=f'Boxplot of {variable_to_plot} by {selected_data_source}',
                 labels={selected_data_source: selected_data_source.replace('_', ' ').title()},
                 color_discrete_map=get_color_map(selected_data_source)
                 )
    if selected_data_source == 'overall_survival':
        # Replace legend labels for 'Overall Survival'
        fig.update_traces(name='Survived', selector=dict(name='1'))
        fig.update_traces(name='Died', selector=dict(name='0'))

    if variable_to_plot:
        selected_variable = variable_to_plot

    return fig


start_column_index = 520

# Select the subset of columns from 'start_column_index' to the end
selected_columns = data.iloc[:, start_column_index:]

# Count non-zero occurrences for each gene in the mutation data
gene_mutation_count = pd.Series([np.sum(selected_columns[col].astype(str).str.strip() != '0') for col in selected_columns], index=selected_columns.columns)

gene = gene_mutation_count.index.str.replace('_mut', '')
val = gene_mutation_count.values
gene_mutation_count_df = pd.DataFrame({'Gene Name': list(gene), 'Mutation Count': list(val)})
gene_mutation_fig = px.bar(gene_mutation_count_df, x='Gene Name', y='Mutation Count',
             title='Gene Mutation Distribution',
             labels={'Gene Name': 'Gene Name', 'Mutation Count': 'Mutation Count'}
            )


# Set the default variable for the boxplot-survival and data source
default_variable = numerical_variables_standardized.columns[0]
default_data_source = 'overall_survival'

layout = html.Div(children=[

    html.Div(className='row', children=[
        html.Div(className='six columns', style={'background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'}, children=[
            dcc.Graph(id='boxplot', figure=boxplot_fig)
        ]),
        html.Div(className='six columns', style={'background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'}, children=[
            dcc.Graph(id='gene-mutation-distplot', figure=gene_mutation_fig)
        ])
    ]),

    html.Label('Select a variable from numerical_variables:'),
    dcc.Dropdown(
        id='variable-dropdown',
        options=[{'label': var, 'value': var} for var in numerical_variables2.columns],
        value=default_variable  # Set the default variable here
    ),

    html.Label('Select data source:'),
    dcc.RadioItems(
        id='data-source-radio',
        options=[
            {'label': 'Overall Survival', 'value': 'overall_survival'},
            {'label': 'Death from Cancer', 'value': 'death_from_cancer'}
        ],
        value=default_data_source,  # Set the default data source 
        labelStyle={'display': 'block'}
    ),

    html.Div(className='row', children=[
        dcc.Store(id='dropdown-value-store'),
        html.Div(className='six columns', style={'background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'}, children=[
            dcc.Graph(id='distplot', 
                      # Set the initial figure for the distplot with the default variable and data source
                      figure=update_distplot(default_variable, default_data_source, None, None, None))
        ]),
        html.Div(className='six columns', style={'background-color': 'white', 'padding': '10px', 'box-shadow': '2px 2px 10px rgba(0, 0, 0, 0.2)'}, children=[
            dcc.Graph(id='boxplot-survival', 
                      # Set the initial figure for the boxplot-survival with the default variable and data source
                      figure=update_boxplot(default_variable, default_data_source, None, None, None))
        ])
    ])

])
