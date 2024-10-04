from qu5it_R import *
import networkx as nx
import matplotlib.pyplot as plt

# Define the graph as directed
G = nx.DiGraph()

# Dictionary to hold the node data (operators)
node_data = {}

# Function to compute the string to prepend to center_node.string to get outer_node.string
def get_prepend_string(center_node, outer_node):
    center_string = center_node.string
    outer_string = outer_node.string
    if center_string == '':
        # If the center node's string is empty, the prepend is the outer node's string
        return outer_string
    elif outer_string.endswith(center_string):
        # If center_string is a suffix of outer_string, compute the prepend
        prepend_length = len(outer_string) - len(center_string)
        prepend = outer_string[:prepend_length]
        return prepend
    else:
        # If center_string is not a suffix, return the full outer_string
        return outer_string

# Add neighbors to the graph and return them, updating the edges with labels
def add_neighbors(graph, center_node):
    if center_node not in graph:
        graph.add_node(center_node)
    
    neighbors = neighbors_mat(center_node)
    
    # Special case: If the node's sde == 0, place all neighbors above
    if center_node.sde == 0:
        greater_than = [n for n in neighbors]  # All neighbors are placed above
        less_than = []  # No node is placed below
    else:
        # Sort the neighbors as usual
        less_than = [n for n in neighbors if n < center_node]
        greater_than = [n for n in neighbors if n > center_node]
    
    # Place the smallest below and the others above (except for when sde == 0)
    edge_labels = {}  # To store the labels for the edges
    if less_than:
        neighbor = less_than[0]
        graph.add_edge(center_node, neighbor)
        # Use get_prepend_string to get the edge label
        edge_label = get_prepend_string(center_node, neighbor)
        edge_labels[(center_node, neighbor)] = edge_label
    for neighbor in greater_than[:3]:
        graph.add_edge(center_node, neighbor)
        # Use get_prepend_string to get the edge label
        edge_label = get_prepend_string(center_node, neighbor)
        edge_labels[(center_node, neighbor)] = edge_label
    
    # Update node_data for the custom display
    for neighbor in neighbors:
        if neighbor not in node_data:
            node_data[neighbor] = neighbor  # The Operator itself as the label
    
    return less_than + greater_than, edge_labels

# Manually position the nodes in a tree-like layout with closer spacing
def get_tree_layout(center_node, neighbors):
    layout = {center_node: (0, 0)}  # Center node at origin

    if center_node.sde == 0:
        # Special case: all neighbors are above if node.sde == 0
        for i, neighbor in enumerate(neighbors[:4]):  # Place all four neighbors above
            layout[neighbor] = (-0.2 + i * 0.2, 0.5)  # Closer horizontally and vertically
    else:
        if len(neighbors) > 0:
            layout[neighbors[0]] = (0, -0.5)  # Below the center node
        for i, neighbor in enumerate(neighbors[1:4]):  # Place the larger neighbors above
            layout[neighbor] = (-0.2 + i * 0.2, 0.5)  # Above the center node

    return layout

# Plot the graph centered on the current node with custom data and edge labels
def plot_graph(graph, center_node, neighbors, edge_labels, highlight_node=None):
    global pos  # Declare pos as global to be accessible
    plt.clf()  # Clear the current plot

    # Get the tree layout - set nodes closer together
    pos = get_tree_layout(center_node, neighbors)

    # Create subgraph consisting only of center_node and its neighbors
    subgraph_nodes = [center_node] + neighbors
    subgraph = graph.subgraph(subgraph_nodes)

    # Draw the subgraph
    nx.draw(subgraph, pos, with_labels=False, node_color='lightblue', node_size=1000, arrows=True)

    # Draw nodes
    nx.draw_networkx_nodes(subgraph, pos, nodelist=[center_node], node_color='black', node_size=0)
    nx.draw_networkx_nodes(subgraph, pos, nodelist=neighbors, node_color='black', node_size=0)

    # Add node.sde when not hovering, full label only when hovering
    labels = {node: str(node.sde) for node in subgraph_nodes}
    font_sizes = {node: 24 for node in labels}

    if highlight_node is not None:
        labels[highlight_node] = node_data[highlight_node]
        font_sizes[highlight_node] = 8

    # Draw edge labels
    nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels, font_size=8)

    # Draw the node labels with different font sizes
    for node, label in labels.items():
        nx.draw_networkx_labels(subgraph, pos, {node: label}, font_size=font_sizes[node])

    # Display circuit length and edges traveled
    circuit_length = len(center_node.string)
    plt.text(1, 0.05, f"Circuit Length: {circuit_length}", horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='lightgrey', alpha=0.5))
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
            if (x - nx_pos) ** 2 + (y - ny_pos) ** 2 < 0.01:
                # Clear the plot
                plt.clf()
                # Redraw the graph with this node's full label
                plot_graph(G, current_node, current_neighbors, current_edge_labels, highlight_node=node)
                plt.draw()
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
        if (x - nx_pos) ** 2 + (y - ny_pos) ** 2 < 0.1:
            clicked_node = node
            break

    if clicked_node is not None and clicked_node != current_node:
        current_node = clicked_node
        current_neighbors, current_edge_labels = add_neighbors(G, current_node)
        plot_graph(G, current_node, current_neighbors, current_edge_labels)

# Initialize the graph
starting_node = go_stupid()
starting_node.string = ''
current_node = starting_node  # Start with the initial node
# Assign custom data to the initial node
node_data[current_node] = current_node  # Label for the starting operator
current_neighbors, current_edge_labels = add_neighbors(G, current_node)

# Set up the plot
fig = plt.figure()
pos = {}  # Initialize pos globally
plot_graph(G, current_node, current_neighbors, current_edge_labels)
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('motion_notify_event', hover)

plt.show()


