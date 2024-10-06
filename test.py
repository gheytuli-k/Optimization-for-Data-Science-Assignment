from aux_functions import Loader
from aux_functions import find_disjoint_cycles


graph = {
    'A': {'B'},
    'B': {'C'},
    'C': {'D'},
    'D': {'A'}
}

# drawing the grpah

import networkx as nx
import matplotlib.pyplot as plt

# drawing the graph
G = nx.DiGraph()
G.add_edges_from([(u, v) for u, neighbors in graph.items() for v in neighbors])

pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_size=700, node_color="skyblue", node_shape="s", alpha=0.5, linewidths=40)
plt.show()


print(graph)
cycles = find_disjoint_cycles(graph, 1)
print(cycles)

# data = Loader("Instance Files\Delorme_50_NDD_Unit_2.txt")


# def build_graph(data: dict[int, int, int, list[dict], list[dict], list[dict]]) -> dict[int, list[dict]]:
#     """
#     returns adjacency list representation of graph
#     with also weight5s of the edges
#     """


        
#     print(graph)

# build_graph(data)

