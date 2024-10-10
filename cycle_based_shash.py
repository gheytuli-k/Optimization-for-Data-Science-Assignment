from aux_functions import Loader
from aux_functions import data_to_adjacency
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import time as time

max_cycle_length = 4

data = Loader("Instance Files\Delorme_50_NoNDD_Weight_2.txt")
start_time = time.time()
# making adjacency out of the data and connecting all the pairs to god donor
adjacency = data_to_adjacency(data, True, 0)
NDDs = [pair['pair_id'] for pair in data['pairs_NDDs']]

G = nx.DiGraph(adjacency)
# Assign weight of 1 to all the edges
# for u, v in G.edges():
#     G[u][v]['weight'] = 1

all_nodes = set(G.nodes())
donors = all_nodes - set(NDDs)

for donor in donors:
    for ndd in NDDs:
        if not G.has_edge(donor, ndd):
            G.add_edge(donor, ndd, weight=1) 

def find_cycles(G, max_length, NDDs): #use built in func
    cts = time.time()
    cycles = []
    for cycle in nx.simple_cycles(G, max_length):
        if len(cycle) <= max_length:
            cycles.append(cycle)
    cyc_time = time.time() - cts    
    return cycles, cyc_time

cycles, cycle_time = find_cycles(G, max_cycle_length, NDDs)

sequences = cycles

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
end_time = time.time()

# Print results
print("Max cycle Length:", max_cycle_length)
if selected_sequences:
    print("Selected Cycles and their weights:")
    for seq in selected_sequences:
        weight = 0
        # Sum weights of consecutive edges in the cycle
        for j in range(len(seq) - 1):
            u = seq[j]
            v = seq[j + 1]
            if G.has_edge(u, v):
                weight += G[u][v]['weight']
        u = seq[-1]
        v = seq[0]
        if G.has_edge(u, v):
            weight += G[u][v]['weight']
        print(f"Cycle {seq} has weight: {weight}")

else:
    print("No sequences were selected.")
print(f"Cycle generation time: {cycle_time:.2f} seconds")
print(f"Time taken: {end_time - start_time:.2f} seconds")