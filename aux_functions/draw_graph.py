import matplotlib.pyplot as plt
import networkx as nx

def draw_graph(graph: dict[str, dict[str, int]]) -> None:
    """
    Draws the graph
    :param graph: adjacency list representation of graph
        - key: node
        - value: Dictionary of neighbors and the weight of the edge
            - key: neighbor
            - value: weight
    :return: None
    """
    G = nx.DiGraph()

    # Add pairs/NDDs as nodes
    for pair in graph['pairs_NDDs']:
        G.add_node(pair['pair_id'], is_ndd=pair['is_NDD'], donor_blood_type=pair['donor_blood_type'], 
        patient_blood_type=pair['patient_blood_type'], patient_vpra=pair['patient_vpra'])

    # Add arcs as edges with weights
    for arc in data['arcs']:
        G.add_edge(arc['donor_pair_id'], arc['patient_pair_id'], weight=arc['weight'])

    # Draw ndoes (pairs/NDDs)
    nd_nodes = [n for n, d in G.nodes(data=True) if d['is_ndd'] == 1]
    pr_nodes = [n for n, d in G.nodes(data=True) if d['is_ndd'] == 0]

    nx.draw_networkx_nodes(G, pos, nodelist=nd_nodes, node_color='red', label='NDD', node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=pr_nodes, node_color='blue', label='Pair', node_size=100)

    # Draw edges with weights
    edges = G.edges(data=True)
    nx.draw_networkx_edges(G, pos, edgelist=edges, arrowstyle='->', arrowsize=10, alpha=0.25)
    edge_labels = {(u, v): f"{d['weight']}" for u, v, d in edges}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    # Draw labels for the nodes
    nx.draw_networkx_labels(G, pos, font_size=10, font_color="white")

    # Display the graph
    plt.legend(scatterpoints=1)
    plt.title('Kidney Exchange Graph')
    plt.show()

    # Create a PyVis network
    net = Network(height='750px', width='100%', directed=True)

    # Convert NetworkX graph to PyVis network
    net.from_nx(G)