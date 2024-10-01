from qutrit import *
import networkx as nx
import matplotlib.pyplot as plt

# Define the graph
G = nx.Graph()

# Dictionary to hold the node data (operators)
node_data = {}

# Add neighbors to the graph and return them, updating the edges with labels
def add_neighbors(graph, center_node):
    if center_node not in graph:
        graph.add_node(center_node)
    
    neighbors = neighbors_mat(center_node)
    
    # Special case: If the node's sde == 0, place all neighbors above
    if center_node.sde == 0:
        greater_than = [n[0] for n in neighbors]  # All neighbors are placed above
        less_than = []  # No node is placed below
    else:
        # Sort the neighbors as usual
        less_than = [n[0] for n in neighbors if n[0] < center_node]
        greater_than = [n[0] for n in neighbors if n[0] > center_node]
    
    # Place the smallest below and the others above (except for when sde == 0)
    edge_labels = {}  # To store the labels for the edges
    if less_than:
        graph.add_edge(center_node, less_than[0])
        edge_labels[(center_node, less_than[0])] = next(t[1] for t in neighbors if t[0] == less_than[0])
    for neighbor in greater_than[:3]:
        graph.add_edge(center_node, neighbor)
        edge_labels[(center_node, neighbor)] = next(t[1] for t in neighbors if t[0] == neighbor)
    
    # Update node_data for the custom display
    for neighbor, label in neighbors:
        if neighbor not in node_data:
            node_data[neighbor] = neighbor  # The Operator itself as the label
    
    return less_than + greater_than, edge_labels

# Manually position the nodes in a tree-like layout with closer spacing
def get_tree_layout(center_node, neighbors):
    layout = {center_node: (0, 0)}  # Center node at origin
    
    if center_node.sde == 0:
        # Special case: all neighbors are above if node.sde == 0
        for i, neighbor in enumerate(neighbors[:4]):  # Place all four neighbors above
            layout[neighbor] = (-0.2 + i * 0.2, 0.5)  # Much closer horizontally and vertically
    else:
        if len(neighbors) > 0:
            layout[neighbors[0]] = (0, -0.5)  # Closer vertically below
        for i, neighbor in enumerate(neighbors[1:4]):  # Place the larger neighbors above
            layout[neighbor] = (-0.2 + i * 0.2, 0.5)  # Much closer horizontally and vertically
    
    return layout

# Plot the graph centered on the current node with custom data and edge labels
def plot_graph(graph, center_node, neighbors, edge_labels):
    global pos  # Declare pos as global to be accessible in on_click
    plt.clf()  # Clear the current plot
    
    # Get the tree layout - set nodes closer together
    pos = get_tree_layout(center_node, neighbors)
    
    # Draw the subgraph centered on the current node
    subgraph = graph.subgraph([center_node] + neighbors)
    
    # Draw nodes and edges
    nx.draw(subgraph, pos, with_labels=False, node_color='lightblue', node_size=500)
    nx.draw_networkx_nodes(subgraph, pos, nodelist=[center_node], node_color='red', node_size=800)  # Highlight center
    nx.draw_networkx_nodes(subgraph, pos, nodelist=neighbors, node_color='green', node_size=600)  # Highlight neighbors
    
    # Add custom labels with smaller font size
    labels = {node: node_data[node] for node in [center_node] + neighbors}
    nx.draw_networkx_labels(subgraph, pos, labels, font_size=8)  # Set font size to a smaller value (e.g., 8)
    
    # Draw edge labels
    nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels, font_size=8)
    
    plt.title(f"Center Node: {center_node}")
    plt.draw()

# Handle click events
def on_click(event):
    global current_node, pos
    # Get the node closest to the click
    try:
        x, y = event.xdata, event.ydata
    except TypeError:
        return  # If click is outside axes, ignore

    clicked_node = None
    # Compare click position with current node and neighbors only (not the full graph)
    for node, (nx_pos, ny_pos) in pos.items():
        if (x - nx_pos) ** 2 + (y - ny_pos) ** 2 < 0.1:  # Check if click is close to a node
            clicked_node = node
            break

    if clicked_node is not None and clicked_node != current_node:
        current_node = clicked_node
        neighbors, edge_labels = add_neighbors(G, current_node)
        plot_graph(G, current_node, neighbors, edge_labels)

# Initialize the graph
current_node = H  # Start with the "H" operator
# Assign custom data to the initial node
node_data[current_node] = current_node  # Label for the "H" operator
neighbors, edge_labels = add_neighbors(G, current_node)

# Set up the plot
fig = plt.figure()
pos = {}  # Initialize pos globally
plot_graph(G, current_node, neighbors, edge_labels)
fig.canvas.mpl_connect('button_press_event', on_click)

plt.show()
