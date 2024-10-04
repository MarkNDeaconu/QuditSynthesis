from qu5it_R import *
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
def plot_graph(graph, center_node, neighbors, edge_labels, highlight_node=None):
    global pos  # Declare pos as global to be accessible
    plt.clf()  # Clear the current plot
    
    # Get the tree layout - set nodes closer together
    pos = get_tree_layout(center_node, neighbors)
    
    # Draw the subgraph centered on the current node
    subgraph = graph.subgraph([center_node] + neighbors)
    
    # Draw nodes and edges
    nx.draw(subgraph, pos, with_labels=False, node_color='lightblue', node_size=1000)
    nx.draw_networkx_nodes(subgraph, pos, nodelist=[center_node], node_color='black', node_size=0)  # Highlight center
    nx.draw_networkx_nodes(subgraph, pos, nodelist=neighbors, node_color='black', node_size=0)  # Highlight neighbors
    
    # Add node.sde when not hovering, full label only when hovering
    labels = {node: str(node.sde) for node in [center_node] + neighbors}
    font_sizes = {node: 24 for node in labels}  # Set font size to be triple (8*3 = 24) for non-hovering text

    if highlight_node is not None:
        labels[highlight_node] = node_data[highlight_node]
        font_sizes[highlight_node] = 8  # Keep the hover text smaller
    
    # Draw edge labels
    nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels, font_size=8)
    
    # Draw the node labels with different font sizes
    for node, label in labels.items():
        nx.draw_networkx_labels(subgraph, pos, {node: label}, font_size=font_sizes[node])
    
    # Calculate circuit length (length of edge_history)
    circuit_length = len(center_node.string)  # Length of the central node's string
    
    # Display the circuit length in a box above the edge history
    plt.text(1, 0.05, f"Circuit Length: {circuit_length}", horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='lightgrey', alpha=0.5))
    
    # Display the central node's string in the edge history box
    plt.text(1, 0, f"Edges Traveled: {center_node.string}", horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='lightgrey', alpha=0.5))
    
    plt.title(f"Center Node: {center_node}")
    plt.draw()

# Show node label when hovering over it, otherwise display node.sde with larger font size
def hover(event):
    if event.inaxes:
        x, y = event.xdata, event.ydata
        # Check which node the mouse is closest to
        for node, (nx_pos, ny_pos) in pos.items():
            if (x - nx_pos) ** 2 + (y - ny_pos) ** 2 < 0.01:  # If the distance is small, hover over this node
                # Clear the plot
                plt.clf()
                # Redraw the graph with this node's full label
                plot_graph(G, current_node, current_neighbors, current_edge_labels, highlight_node=node)
                return  # Exit after drawing the hovered node's label

    # If not hovering over any node, just show node.sde for all nodes (with large font size)
    plot_graph(G, current_node, current_neighbors, current_edge_labels)
    plt.draw()

# Handle click events and update edge history
def on_click(event):
    global current_node, current_neighbors, current_edge_labels
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
        current_neighbors, current_edge_labels = add_neighbors(G, current_node)
        plot_graph(G, current_node, current_neighbors, current_edge_labels)

# Initialize the graph
starting_node = go_stupid()
generating_string = starting_node.string
starting_node.string = ''
current_node = starting_node  # Start with the "H" operator
# Assign custom data to the initial node
node_data[current_node] = current_node  # Label for the "H" operator
current_neighbors, current_edge_labels = add_neighbors(G, current_node)

# Set up the plot
fig = plt.figure()
pos = {}  # Initialize pos globally
plot_graph(G, current_node, current_neighbors, current_edge_labels)
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('motion_notify_event', hover)

plt.show()











