

#IN PROGRESS DUE TO REFACTORING
#DOES NOT RUN

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def animate(starting_node, edges):

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

    def update_visualization(central_node):
        global G, pos, central_node_label, neighbor_nodes_info, edge_label_dict, current_hover, labels, node_labels

        plt.clf()  # Clear the current figure

        # Get neighbors_mat
        # neighbor_nodes = central_node.neighbors_mat(edges, edges)
        neighbor_nodes = central_node.neighbors_mat(edges, edges)
        edge_labels_list = ['C', 'HRC', 'HRHHRC']  # Ensure this list matches the neighbors_mat

        # Create graph
        G = nx.Graph()

        # Use node labels as identifiers
        central_node_label = str(central_node)
        
        # Initialize lists for neighbor nodes info
        neighbor_nodes_info = []

        # Collect neighbor nodes along with their edge labels
        for neighbor, edge_label in zip(neighbor_nodes, edge_labels_list):
            neighbor_nodes_info.append({
                'node': neighbor,
                'label': str(neighbor),
                'sde': neighbor.sde,
                'edge_label': edge_label
            })

        # Classify neighbor nodes based on comparison
        greater_nodes_info = []
        less_nodes_info = []

        for info in neighbor_nodes_info:
            neighbor = info['node']
            if neighbor == central_node:
                pass
            elif not neighbor.__lt__(central_node):
                greater_nodes_info.append(info)
            else:
                less_nodes_info.append(info)

        # Initialize the edge label dictionary
        edge_label_dict = {}

        # Add the central node
        G.add_node(central_node_label, sde=central_node.sde)

        # Add greater nodes above and less nodes below the central node
        for info in greater_nodes_info + less_nodes_info:
            G.add_node(info['label'], sde=info['sde'])
            edge_key = (central_node_label, info['label'])
            G.add_edge(central_node_label, info['label'])
            edge_label_dict[edge_key] = info['edge_label']

        # Set the positions
        pos = {central_node_label: (0, 0)}  # Central node at the center

        # Position greater nodes equally spaced above the central node
        n_greater = len(greater_nodes_info)
        if n_greater > 0:
            x_coords_greater = np.linspace(-1, 1, n_greater)
            y_coords_greater = np.full(n_greater, 1)
            for i, info in enumerate(greater_nodes_info):
                pos[info['label']] = (x_coords_greater[i], y_coords_greater[i])

        # Position less nodes centered below the central node
        n_less = len(less_nodes_info)
        if n_less > 0:
            x_coords_less = np.linspace(-1, 1, n_less) if n_less > 1 else [0]
            y_coords_less = np.full(n_less, -1)
            for i, info in enumerate(less_nodes_info):
                pos[info['label']] = (x_coords_less[i], y_coords_less[i])

        # Draw nodes and edges only once
        nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=500)
        nx.draw_networkx_edges(G, pos)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_label_dict)

        # Draw initial node labels with sde
        node_labels = {node: G.nodes[node]['sde'] for node in G.nodes()}
        labels = nx.draw_networkx_labels(G, pos, labels=node_labels)

        # Remove axis
        plt.axis('off')

        # Set plot limits to expand to the full screen
        plt.xlim(-2, 2)  # Expands left and right limits
        plt.ylim(-2, 2)  # Expands top and bottom limits

        # Display text at the bottom
        text = get_prepend_string(starting_node, central_node)
        plt.text(0, -1.8, str(len(text)) + ' - ' + text, ha='center', va='top')

        # Use tight layout to further maximize screen space usage
        plt.tight_layout(pad=0)

        # Connect hover event with updated labels and node_labels
        current_hover = None  # Track the currently hovered node to prevent unnecessary redraws
        fig.canvas.mpl_connect('motion_notify_event', lambda event: smooth_hover(event))
        
        # Redraw the canvas
        plt.draw()

    def on_click(event):
        global central_node
        ax = plt.gca()
        if event.inaxes == ax:
            for node_label in G.nodes():
                x, y = pos[node_label]
                dx = event.xdata - x
                dy = event.ydata - y
                distance = np.sqrt(dx**2 + dy**2)
                if distance < 0.1:  # Adjust threshold as needed
                    if node_label == central_node_label:
                        selected_node = central_node
                    else:
                        selected_node = next(info['node'] for info in neighbor_nodes_info if info['label'] == node_label)
                    central_node = selected_node
                    update_visualization(selected_node)
                    break

    def smooth_hover(event):
        global current_hover, labels, node_labels
        ax = plt.gca()
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
            closest_node = None

            for node_label, (nx_pos, ny_pos) in pos.items():
                if (x - nx_pos) ** 2 + (y - ny_pos) ** 2 < 0.01:
                    closest_node = node_label
                    break

            if closest_node != current_hover:
                # Update the label for the newly hovered node
                current_hover = closest_node
                for node, text_obj in labels.items():
                    text_obj.set_text(node_labels[node] if node != current_hover else node)  # Show full label on hover

                plt.draw()

    # Initial central node
    # starting_node = z3.from_orbit([H,S,R])

    # Create a global figure with fullscreen mode
    fig, ax = plt.subplots(figsize=(16, 9))  # Full HD aspect ratio
    fig.canvas.manager.full_screen_toggle()   # Toggle fullscreen mode

    # Initialize the visualization with the starting node
    update_visualization(starting_node)

    # Connect the click event
    fig.canvas.mpl_connect('button_press_event', on_click)

    # Show the plot
    plt.show()









