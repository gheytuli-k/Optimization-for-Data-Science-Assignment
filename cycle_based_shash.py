#imports and functions from aux_functions
from aux_functions import Loader
from aux_functions import data_to_adjacency
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import time as time

def load_data(data_path):
    data = Loader(data_path)

    #dict with pair_id as key and patient_vpra as value
    pair_vpra = {}
    for pair in data['pairs_NDDs']:
        pair_vpra[pair['pair_id']] = pair['patient_vpra']
    for pair in data['pairs_NoNDDs']:
        pair_vpra[pair['pair_id']] = pair['patient_vpra']
    data['pair_vpra'] = pair_vpra

    return data

def create_graph(data, god_donor=True, god_donor_weight=0):
    #making adjacency out of the data and connecting all the pairs to god donor
    adjacency = data_to_adjacency(data, god_donor, god_donor_weight)
    NDDs = [pair['pair_id'] for pair in data['pairs_NDDs']]
    G = nx.DiGraph(adjacency)

    #connect all the pairs to god donor with weight god_donor_weight
    all_nodes = set(G.nodes())
    donors = all_nodes - set(NDDs)
    for donor in donors:
        for ndd in NDDs:
            if not G.has_edge(donor, ndd):
                G.add_edge(donor, ndd, weight=god_donor_weight)
    return G

def find_cycles(G, max_length): #use built in func
    cts = time.time()
    cycles = []

    #only generate cycle of length <= max_length
    for cycle in nx.simple_cycles(G, max_length):
        if len(cycle) <= max_length:
            cycles.append(cycle)
    cyc_time = time.time() - cts    
    return cycles, cyc_time

def optimize_sequences(sequences, G, beta, pair_vpra):
    #initialise model
    model = gp.Model("Kidney_Exchange_Optimization")

    # variables and weights
    sequence_vars = {} #dict with cycle as key and variable as value
    sequence_weights = {} #dict with cycle as key and weight as value
    for i, seq in enumerate(sequences):
        weight = 0
        for j in range(len(seq) - 1):
            u = seq[j]
            v = seq[j + 1]
            if G.has_edge(u, v): #for each edge in the cycle
                edge_weight = G[u][v]['weight']
                #adjust weight based on beta (and pair_vpra)
                adjusted_weight = edge_weight * (1 + beta * pair_vpra.get(v, 0)) #if no vpra, use 0
                weight += adjusted_weight
            else:
                print(f"Edge from {u} to {v} not found in graph.")

        #add the weight of the last edge to complete the cycle
        if seq[0] == seq[-1] and len(seq) > 1:
            u = seq[-1]
            v = seq[0]
            if G.has_edge(u, v):
                edge_weight = G[u][v]['weight']
                adjusted_weight = edge_weight * (1 + beta * pair_vpra.get(v, 0))
                weight += adjusted_weight
        sequence_weights[i] = weight
        #add binary variable for each cycle
        sequence_vars[i] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}")

    # objective function -> maximise total weight of selected cycles (inlcudes beta and pair_vpra factor) 
    model.setObjective(
        gp.quicksum(sequence_weights[i] * sequence_vars[i] for i in sequence_vars),
        GRB.MAXIMIZE
    )

    # constraints for nodes - ensure each node is in at most one cycle
    vertex_sequences = {}
    for i, seq in enumerate(sequences):
        seq_nodes = set(seq)
        for vertex in seq_nodes: #for each vertex in the cycle add the cycle index to the vertex_sequences dict
            if vertex not in vertex_sequences: #if not there add to dict
                vertex_sequences[vertex] = []
            vertex_sequences[vertex].append(i)

    for vertex, seq_indices in vertex_sequences.items():
        model.addConstr(
            gp.quicksum(sequence_vars[i] for i in seq_indices) <= 1,
            name=f"vertex_{vertex}_constraint"
        )

    #Optimise
    model.optimize()
    
    #retrieve selected sequences
    selected_sequences = [sequences[i] for i in sequence_vars if sequence_vars[i].X > 0.5]
    total_weight = sum(sequence_weights[i] for i in sequence_vars if sequence_vars[i].X > 0.5)
    print(f"Total weight of selected sequences (model): {total_weight}")
    return selected_sequences

def output_results(selected_sequences, G, max_cycle_length, cycle_time, start_time):
    print("Set max cycle Length:", max_cycle_length)
    if selected_sequences:
        print("Selected Cycles and their weights:")
        total_weight = 0
        for seq in selected_sequences:
            weight = 0
            #sum weights of consecutive edges in the cycle
            for j in range(len(seq) - 1):
                u = seq[j]
                v = seq[j + 1]
                if G.has_edge(u, v):
                    weight += G[u][v]['weight']
            u = seq[-1] #add last edge to complete the cycle
            v = seq[0]
            if G.has_edge(u, v):
                weight += G[u][v]['weight']
                total_weight += weight
            print(f"Cycle {seq} has weight: {weight}") #per cycle weight
        print(f"Total weight of selected cycles: {total_weight}") #total weight (no beta or vpra factor)
    else:
        print("No sequences were selected.")
    print(f"Cycle generation time: {cycle_time:.2f} seconds")
    print(f"Total time taken: {time.time() - start_time:.2f} seconds")

def cycle_based_based_optimisation(data_path, max_cycle_length: int = 3, god_donor: bool = True, god_donor_weight: int = 0, beta: float = 0.0):
    """
    Optimizes the cycles in the graph
    :param data_path: Path to the data file
    :param max_cycle_length: Maximum length of cycles to consider
    :return: None (but prints results)
    """
    start_time = time.time()
    
    #load data
    data = load_data(data_path)
    
    #create graph
    G = create_graph(data, god_donor, god_donor_weight)

    #find cycles
    cycles, cycle_time = find_cycles(G, max_cycle_length)

    #optimize
    selected_sequences = optimize_sequences(cycles, G, beta, data['pair_vpra'])
    
    #output results
    output_results(selected_sequences, G, max_cycle_length, cycle_time, start_time)


cycle_based_based_optimisation("Instance Files\Saidman_50_NDD_Weight_2.txt", 3, True, 0, 0)

    
    

