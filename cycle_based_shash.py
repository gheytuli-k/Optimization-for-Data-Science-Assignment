from aux_functions import Loader
from aux_functions import data_to_adjacency
import gurobipy as gp
from gurobipy import GRB
import networkx as nx

data = Loader("Instance Files\Delorme_50_NDD_Unit_2.txt")
# making adjacency out of the data and connecting all the pairs to god donor
adjacency = data_to_adjacency(data, True, 1)
NDDs = [pair['pair_id'] for pair in data['pairs_NDDs']]

# precomputing the cycles upto length 5 in the network
def find_cycles(G, max_length):
    cycles = []
    
    def dfs(current_node, path):
        if len(path) > max_length:
            return
        
        for neighbor in G.successors(current_node):
            if neighbor == path[0] and len(path) > 2:  # Valid cycle found
                cycles.append(path + [neighbor])
            elif neighbor not in path:
                dfs(neighbor, path + [neighbor])

    for node in G.nodes():
        dfs(node, [node])
    
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
cycles = find_cycles(G, 6)

max_cycle_length = 6  
max_chain_length = 10 
cycles = find_cycles(G, max_cycle_length)
chains = find_chains(G, NDDs, max_chain_length)

sequences = cycles + chains

def optimize_sequences(sequences, G, NDDs):
    model = gp.Model("Kidney_Exchange_Optimization")

    # Variables and weights
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

    # Objective function
    model.setObjective(
        gp.quicksum(sequence_weights[i] * sequence_vars[i] for i in sequence_vars),
        GRB.MAXIMIZE
    )

    # Constraints for nodes
    vertex_sequences = {}
    for i, seq in enumerate(sequences):
        for vertex in set(seq):
            if vertex not in vertex_sequences:
                vertex_sequences[vertex] = []
            vertex_sequences[vertex].append(i)

    for vertex, seq_indices in vertex_sequences.items():
        model.addConstr(
            gp.quicksum(sequence_vars[i] for i in seq_indices) <= 1,
            name=f"vertex_{vertex}_constraint"
        )

    # Constraints for NDDs (if necessary)
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