Plotly Dash METABRIC Data Visualization App

Overview
This app allows users to interact with METABRIC Breast Cancer Dataset and creates interactive Plotly visualizations from that data. The visualizations update in real-time as the user interacts with dropdowns, sliders, and other controls.

Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

Prerequisites
You will need Python 3 and some libraries installed. It is recommended you install these into a virtual environment.

### Set up virtual environment

It is recommended to use a virtual environment for this project. This keeps the dependencies and packages for this project isolated from other projects you may have.

To create a virtual environment:

python3 -m venv dash-env

This will create a folder called `dash-env` in your current directory. 

Activate the environment:

**Linux/Mac:**
source dash-env/bin/activate

**Windows:**
dash-env\Scripts\activate

Your command prompt should now have `(dash-env)` prefixed.

### Install requirements

With the virtual environment activated, install the required packages:
pip install -r requirements.txt

This will install `dash`, `plotly` and other required packages into the virtual environment.

Now you can run and develop the app within this virtual environment.

To deactivate the environment when done:

**Linux/Mac:**
deactivate

**Windows:**
dash-env\Scripts\deactivate

Running the app
To run the app in development mode:
python app.py

Open http://localhost:8050 to view the app in your browser.

The app will reload automatically as you make changes to the code.

Built With
Dash - Main framework
Plotly Python - Used to create the interactive visualizations

Author
Myroslava Liubchenko
