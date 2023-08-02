from dash import Dash, html, dcc
import dash
import os

custom_css = '''
    body {
        font-family: "Yeseva One Regular", sans-serif;
    }
'''

# Append the custom CSS style to the external_stylesheets list
external_stylesheets = [
    'https://codepen.io/chriddyp/pen/bWLwgP.css',  # Existing stylesheet
    {'href': 'https://fonts.googleapis.com/css2?family=Yeseva+One&display=swap', 'rel': 'stylesheet'},  # Yeseva One font
    {'type': 'text/css', 'style': custom_css},  # Custom CSS style for the app
]
# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, use_pages=True, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.Div(
        [
            html.Img(
                src="/assets/logo-no-background.svg",
                style={
                    "display": "inline-block",
                    "margin": "auto",
                    "padding-top": "20px",
                    "width": "400px",
                    "background-color": "#003262",
                },
            ),
            html.H3(
                'Breast Cancer Gene Expression Profiles (METABRIC)',
                style={
                    'text-align': 'center',
                    'background-color': '#003262',
                    'color': 'white',
                    'padding': '15px 0',
                    'margin-top': 0,  # Remove default margin for H1
                    'margin-bottom': '0px',
                },
            ),
            html.P(
                'The genetics part of the dataset contains m-RNA levels z-score for 331 genes, and mutation for 175 genes.',
                style={
                    'text-align': 'center',
                    'margin-bottom': '30px',
                    'background-color': '#003262',
                    'color': 'white',
                    'padding': '0px',
                },
            ),
        ],
        style={
            'background-color': '#003262',  # Apply background color to the entire top section
        }
    ),
    html.Div(
        [
            dcc.Link(
                f"{page['name']}",
                href=page["relative_path"],
                className="tab-link",
                style={
                    "display": "inline-block",
                    "padding": "10px",
                    "margin-right": "5px",
                    "background-color": "lightgray",
                    "border-top-left-radius": "5px",
                    "border-top-right-radius": "5px",
                    "border": "1px solid gray",
                    "color": "black",
                    "text-decoration": "none",
                    "font-weight": "bold",
                    "width": "48.14%",
                    "cursor": "pointer",
                    "color": "black",
                    "background-color": "lightgray",
                    "hover": {
                        "color": "white",
                        "background-color": "#4285f4"
                    }
                }
            )
            for page in dash.page_registry.values()
        ]
    ),

    dash.page_container,
    
html.Div(
        [
            html.Img(
                src="/assets/logo-no-background.svg",
                style={
                    "display": "inline-block",
                    "margin": "auto",
                    "padding-top": "20px",
                    "width": "400px",
                    "background-color": "#003262",
                },
            ),
            html.P(
                "The project was completed as a master's dissertation by Liubchenko Myroslava, a student at the University of Glasgow.",
                style={
                    'text-align': 'left',
                    'margin-bottom': '5px',
                    'background-color': '#003262',
                    'color': 'white',
                    'padding': '5px',
                },
            ),
            html.P(
                "For inquiries, contact at: ",
                style={
                    'text-align': 'left',
                    "display": "inline",
                    'margin-bottom': '0',  # Reduce the margin to 0 for the last paragraph
                    'background-color': '#003262',
                    'color': 'white',
                    'padding': '10px',
                },
            ),
            html.A(
                "2829786L@student.gla.ac.uk",
                href="mailto:2829786L@student.gla.ac.uk",
                style={
                    'text-align': 'left',
                    'margin': '0',  # Remove the margin of the anchor element
                    'line-height': '1.5',  # Add line-height property to control spacing between lines
                    "display": "inline",
                    'background-color': '#003262',
                    'color': 'white',
                    'padding': '10px',
                    'text-decoration': 'none',
                },
            ),
        ],
        style={
            'background-color': '#003262',  # Apply background color to the entire footer section
            'margin-bottom': '0',  # Set the margin at the bottom to 0 to remove any extra gap
        }
    ),
])

if __name__ == '__main__':
    app.run_server(debug=True)
