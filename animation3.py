from qutrit import *
import networkx as nx
import numpy as np
import dash
import dash_cytoscape as cyto
from dash import html
from dash.dependencies import Input, Output

# Function to generate nodes and edges in a Cytoscape-compatible format
def generate_elements(central_node):
    elements = []
    neighbor_nodes = neighbors_mat(central_node)
    edge_labels_list = ['H', 'HR', 'HRHH', 'HRHHR']

    # Store full matrix representations for detailed labels
    global full_labels
    full_labels = {}  # Reset the full_labels dictionary

    # Central node position
    central_node_label = f'node_{hash(str(central_node)) % 10000}'
    full_labels[central_node_label] = str(central_node)  # Store the full label
    elements.append({'data': {'id': central_node_label, 'label': str(central_node.sde)},
                     'position': {'x': 0, 'y': 0}})  # Center position for the central node
    
    # Set up positions for neighbors
    n_neighbors = len(neighbor_nodes)
    
    for idx, (neighbor, edge_label) in enumerate(zip(neighbor_nodes, edge_labels_list)):
        neighbor_label = f'node_{hash(str(neighbor)) % 10000}'
        full_labels[neighbor_label] = str(neighbor)  # Store the full label
        
        # Position the first neighbor directly below
        if idx == 0:
            elements.append({'data': {'id': neighbor_label, 'label': str(neighbor.sde)},
                             'position': {'x': 0, 'y': -200}})
        else:
            # Position remaining neighbors above and spaced horizontally
            x_pos = (idx - 1) * 200 - (100 * (n_neighbors - 2))  # Adjust spacing horizontally
            elements.append({'data': {'id': neighbor_label, 'label': str(neighbor.sde)},
                             'position': {'x': x_pos, 'y': 200}})
        
        # Add the edge connecting to the central node
        elements.append({
            'data': {
                'source': central_node_label,
                'target': neighbor_label,
                'label': edge_label
            }
        })

    return elements

# Initialize central node
starting_node = go_stupid()
central_node = starting_node

# Initialize the app
app = dash.Dash(__name__)

# Define the layout of the app
app.layout = html.Div([
    html.H2("Interactive Network Graph with Dash Cytoscape"),
    cyto.Cytoscape(
        id='cytoscape-graph',
        layout={'name': 'preset'},  # Using preset layout for custom positions
        style={'width': '100%', 'height': '800px'},
        elements=generate_elements(central_node),
        stylesheet=[
            {
                'selector': 'node',
                'style': {
                    'content': 'data(label)',
                    'background-color': '#0074D9',
                    'color': 'white',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'width': 50,
                    'height': 50
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'content': 'data(label)',
                    'width': 2,
                    'line-color': '#888',
                    'target-arrow-color': '#888',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier'
                }
            }
        ]
    ),
    html.Div(id='hover-output', style={'padding': '10px', 'font-size': '16px', 'position': 'absolute', 'top': '0', 'right': '0'}),
    html.Div(id='click-output', style={'padding': '10px', 'font-size': '16px'})
])

# Callback to update the graph and hover information when a node is clicked or hovered
@app.callback(
    Output('cytoscape-graph', 'elements'),
    Output('hover-output', 'children'),
    Output('click-output', 'children'),
    Input('cytoscape-graph', 'tapNodeData')
)
def update_graph(node_data):
    global central_node
    
    if node_data:
        selected_node_label = node_data['id']
        # Find the corresponding node object
        if selected_node_label == f'node_{hash(str(central_node)) % 10000}':
            selected_node = central_node
        else:
            # Find the neighbor matching the selected label
            selected_node = next(
                neighbor for neighbor in neighbors_mat(central_node)
                if f'node_{hash(str(neighbor)) % 10000}' == selected_node_label
            )
        central_node = selected_node

        # Update the elements of the graph
        elements = generate_elements(central_node)
        
        # Format the full label for display
        hover_text = full_labels[selected_node_label].replace('\n', '<br>')
        click_text = f"Clicked on node: {selected_node_label}"
    else:
        elements = generate_elements(central_node)
        hover_text = "Hover over a node to view its full label."
        click_text = "Click on a node to select it."
    
    return elements, hover_text, click_text

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)





