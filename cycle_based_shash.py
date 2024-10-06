from aux_functions import Loader
from aux_functions import data_to_adjacency
import gurobipy as gp
from gurobipy import GRB
import networkx as nx

data = Loader("Instance Files\Delorme_50_NDD_Unit_2.txt")
# making adjacency out of the data and connecting all the pairs to god donor
adjacency = data_to_adjacency(data, True, 0)
NDDs = [pair['pair_id'] for pair in data['pairs_NDDs']]

# def find_cycles(G, max_length, NDDs):
#     cycles = []
    
#     def dfs(current_node, path):
#         if len(path) > max_length:
#             return
        
#         for neighbor in G.successors(current_node):
#             if neighbor in NDDs:
#                 continue  # skip NDDs
#             if neighbor == path[0] and len(path) > 2:  # valid cycle
#                 cycles.append(path + [neighbor])
#             elif neighbor not in path:
#                 dfs(neighbor, path + [neighbor])

#     for node in G.nodes():
#         if node in NDDs:
#             continue  # skip NDDs
#         dfs(node, [node])
    
#     return cycles

def find_cycles(G, max_length, NDDs):
    cycles = []
    for cycle in nx.simple_cycles(G):
        if len(cycle) <= max_length and not any(node in NDDs for node in cycle):
            cycles.append(cycle)
    return cycles


def find_chains(G, NDDs, max_chain_length):
    chains = []

    def dfs_chain(current_node, path):
        if len(path) > max_chain_length:
            return
        for neighbor in G.successors(current_node):
            if neighbor not in path:
                chains.append(path + [neighbor])
                dfs_chain(neighbor, path + [neighbor])

    for ndd in NDDs:
        dfs_chain(ndd, [ndd])

    return chains

G = nx.DiGraph(adjacency)

# Assign weight of 1 to all the edges
for u, v in G.edges():
    G[u][v]['weight'] = 1

# Remove edges pointing to NDDs
for ndd in NDDs:
    incoming_edges = list(G.in_edges(ndd))
    G.remove_edges_from(incoming_edges)

max_cycle_length = 15 
max_chain_length = 20 

cycles = find_cycles(G, max_cycle_length, NDDs)
chains = find_chains(G, NDDs, max_chain_length)

sequences = cycles + chains

def optimize_sequences(sequences, G, NDDs):
    model = gp.Model("Kidney_Exchange_Optimization")

    # variables and weights
    sequence_vars = {}
    sequence_weights = {}
    for i, seq in enumerate(sequences):
        weight = 0
        for j in range(len(seq) - 1):
            u = seq[j]
            v = seq[j + 1]
            if G.has_edge(u, v):
                weight += G[u][v]['weight']
            else:
                print(f"Edge from {u} to {v} not found in graph.")
        sequence_weights[i] = weight
        sequence_vars[i] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}")

    # objective
    model.setObjective(
        gp.quicksum(sequence_weights[i] * sequence_vars[i] for i in sequence_vars),
        GRB.MAXIMIZE
    )

    # constraints for nodes (excluding NDDs)
    vertex_sequences = {}
    for i, seq in enumerate(sequences):
        # exclude NDDs from vertex constraints
        seq_nodes = set(seq) - set(NDDs)
        for vertex in seq_nodes:
            if vertex not in vertex_sequences:
                vertex_sequences[vertex] = []
            vertex_sequences[vertex].append(i)

    for vertex, seq_indices in vertex_sequences.items():
        model.addConstr(
            gp.quicksum(sequence_vars[i] for i in seq_indices) <= 1,
            name=f"vertex_{vertex}_constraint"
        )

    # Constraints for NDDs
    for ndd in NDDs:
        ndd_sequences = [i for i, seq in enumerate(sequences) if seq[0] == ndd]
        model.addConstr(
            gp.quicksum(sequence_vars[i] for i in ndd_sequences) <= 1,
            name=f"NDD_{ndd}_constraint"
        )

    # Optimize
    model.optimize()

    # Retrieve selected sequences
    selected_sequences = [sequences[i] for i in sequence_vars if sequence_vars[i].X > 0.5]
    total_weight = sum(sequence_weights[i] for i in sequence_vars if sequence_vars[i].X > 0.5)
    print(f"Total weight of selected sequences: {total_weight}")
    return selected_sequences

selected_sequences = optimize_sequences(sequences, G, NDDs)

# Print results
print("Selected Sequences (Cycles and Chains):")
for seq in selected_sequences:
    print(seq)
