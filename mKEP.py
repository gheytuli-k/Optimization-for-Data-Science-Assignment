from aux_functions import Loader
from aux_functions import data_to_adjacency
import random
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from gurobipy import Model
import time

# data_delorme_1000_NDD_Unit_0 = Loader("Instance Files/Delorme_50_NDD_Unit_0.txt")
# adjacency_delorme_1000_NDD_Unit_0 = data_to_adjacency(data_delorme_1000_NDD_Unit_0, True, 1)

def split_nodes(nodes, proportions):
    subset_nodes = []
    total_nodes = len(nodes)
    
    # Calculate the exact number of nodes to assign based on proportions
    allocated_count = [int(total_nodes * prop) for prop in proportions]
    
    # Adjust any rounding differences to ensure all nodes are assigned
    difference = total_nodes - sum(allocated_count)
    
    # Distribute the remaining nodes due to rounding
    for i in range(difference):
        allocated_count[i] += 1
    
    # Split the nodes based on the adjusted counts
    start_idx = 0
    for count in allocated_count:
        end_idx = start_idx + count
        subset_nodes.append(nodes[start_idx:end_idx])
        start_idx = end_idx

    return subset_nodes

def split_adjacency_list(adj_list: dict, n_subsets: int, proportions: list[float], ndd_count: int):

    random.seed(0)
    all_nodes = list(adj_list.keys())
    ndd_nodes = all_nodes[-ndd_count:]
    non_ndd_nodes = all_nodes[:-ndd_count]

    random.shuffle(ndd_nodes)
    random.shuffle(non_ndd_nodes)  

    subset_non_ndd = split_nodes(non_ndd_nodes, proportions)
    subset_ndd = split_nodes(ndd_nodes, proportions)

    subset_adj_list = [{} for _ in range(n_subsets)]

    for i in range(n_subsets):
        for node in subset_non_ndd[i]:
            subset_adj_list[i][node] = {}
            for k, v in adj_list[node].items():
                if k in subset_non_ndd[i] or k in subset_ndd[i]:
                    if node not in subset_adj_list[i]:
                        subset_adj_list[i][node] = {}
                    subset_adj_list[i][node][k] = v
        for node in subset_ndd[i]:
            subset_adj_list[i][node] = {}
            for k, v in adj_list[node].items():
                if k in subset_ndd[i] or k in subset_non_ndd[i]:
                    if node not in subset_adj_list[i]:
                        subset_adj_list[i][node] = {}
                    subset_adj_list[i][node][k] = v
    return subset_adj_list







########################

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
    return G, NDDs

def find_cycles(G, max_length): #use built in func
    cts = time.time()
    cycles = []

    #only generate cycle of length <= max_length
    for cycle in nx.simple_cycles(G, max_length):
        if len(cycle) <= max_length:
            cycles.append(cycle)
    cyc_time = time.time() - cts    
    return cycles, cyc_time

def optimize_sequences(sequences, G, NDDs, beta, pair_vpra):
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

def cycle_based_mkep(data_path, country_proportions, max_cycle_length=3, god_donor=True, god_donor_weight=0, beta=0.0):

    start_time = time.time()

    data = load_data(data_path)

    # Create the full graph with NDD nodes
    G, NDDs = create_graph(data, god_donor, god_donor_weight)
    print(G.edges(data=True))
    # Split nodes into countries based on provided proportions (excluding NDD nodes)
    all_nodes = list(set(G.nodes()) - set(NDDs))
    random.shuffle(all_nodes)
    node_partitions = split_nodes(all_nodes, country_proportions)

    random.shuffle(NDDs)
    NDD_partitions = split_nodes(NDDs, country_proportions)

    print(node_partitions[0])

    # Add the NDD nodes to their respective country partitions
    for i, NDD_subset in enumerate(NDD_partitions):
        node_partitions[i].extend(NDD_subset)
    
    # Create mapping of nodes to countries
    node_to_country = {}
    for i, country_nodes in enumerate(node_partitions):
        for node in country_nodes:
            node_to_country[node] = i
    
    # Generate all cycles once
    full_cycles, full_cycle_time = find_cycles(G, max_cycle_length)

    print(f"Total number of cycles found in the full graph: {len(full_cycles)}")
    
    # Filter cycles for each country
    country_cycles = [[] for _ in country_proportions]
    for seq in full_cycles:
        countries_in_cycle = set(node_to_country.get(node) for node in seq)
        if len(countries_in_cycle) == 1:
            country_index = countries_in_cycle.pop()
            country_cycles[country_index].append(seq)
    
    country_weights = []
    # Proceed with country-level optimizations
    for i, cycles in enumerate(country_cycles):
        print(f"\nOptimizing for country {i+1} with {len(node_partitions[i])} nodes...")
        country_nodes = node_partitions[i]
        country_G = G.subgraph(country_nodes).copy()
        
        NDD_subset = NDD_partitions[i]
        
        # Optimize sequences for this country
        selected_sequences = optimize_sequences(cycles, country_G, NDD_subset, beta, data['pair_vpra'])
    
        # Output results for this country 
        output_results(selected_sequences, country_G, max_cycle_length, full_cycle_time, start_time)
    
        # Calculate total weight achieved by the country
        total_weight = sum(
            sum(
                G[u][v]['weight'] * (1 + beta * data['pair_vpra'].get(v, 0))
                for u, v in zip(seq, seq[1:] + [seq[0]])
            )
            for seq in selected_sequences
        )
        print(f"Total weight achieved by country {i+1}: {total_weight}")
        print("#######################################################################")
    
        country_weights.append(total_weight)

    # Global optimization: use the same cycles
    print("\nProceeding with global optimization using the same cycles...\n")
    
    # Variables and weights for the global graph
    model = gp.Model("mKEP_Global_Optimization")

    sequence_vars = {}
    sequence_weights = {}

    for i, seq in enumerate(full_cycles):
        weight = sum(
            G[u][v]['weight'] * (1 + beta * data['pair_vpra'].get(v, 0)) for u, v in zip(seq, seq[1:] + [seq[0]])
        )
        sequence_vars[i] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}")
        sequence_weights[i] = weight

    # Objective function
    model.setObjective(
        gp.quicksum(sequence_weights[i] * sequence_vars[i] for i in sequence_vars),
        GRB.MAXIMIZE
    )

    # Global constraint: each node should appear in at most one cycle
    vertex_sequences = {}
    for i, seq in enumerate(full_cycles):
        for vertex in set(seq):
            if vertex not in vertex_sequences:
                vertex_sequences[vertex] = []
            vertex_sequences[vertex].append(i)
    
    for vertex, seq_indices in vertex_sequences.items():
        model.addConstr(
            gp.quicksum(sequence_vars[i] for i in seq_indices) <= 1,
            name=f"vertex_{vertex}_constraint"
        ) 

    # Country-specific constraints: ensure each country gets at least as much weight as their individual solution
    for i, country_weight in enumerate(country_weights):
        country_sequences = []
        for j, seq in enumerate(full_cycles):
            countries_in_seq = set(node_to_country.get(node) for node in seq)
            if i in countries_in_seq:
                country_sequences.append(j)
        model.addConstr(
            gp.quicksum(sequence_weights[j] * sequence_vars[j] for j in country_sequences) >= country_weight,
            name=f"country_{i+1}_constraint"
        )

    # Optimize globally
    model.optimize()

    # Retrieve selected sequences
    if model.Status == GRB.OPTIMAL or model.Status == GRB.TIME_LIMIT:
        selected_sequences = [full_cycles[i] for i in sequence_vars if sequence_vars[i].X > 0.5]
        total_global_weight = sum(sequence_weights[i] for i in sequence_vars if sequence_vars[i].X > 0.5)

        print(f"\nTotal weight of selected sequences (global model): {total_global_weight}")
        print(f"Global cycle generation time: {full_cycle_time:.2f} seconds")
        print(f"Total time taken for global optimization: {time.time() - start_time:.2f} seconds")

        # Count the number of nodes selected per country
        nodes_selected = set()
        for seq in selected_sequences:
            nodes_selected.update(seq)

        country_node_counts = {i: 0 for i in range(len(country_proportions))}
        for node in nodes_selected:
            country_index = node_to_country.get(node)
            if country_index is not None:
                country_node_counts[country_index] += 1

        print("\nNumber of nodes selected per country in the final global model:")
        for i in range(len(country_proportions)):
            print(f"Country {i+1}: {country_node_counts[i]} nodes selected")

        return selected_sequences, total_global_weight
    else:
        print("No optimal solution found.")
        return None, None


cycle_based_mkep("Instance Files/Delorme_200_NDD_Unit_0.txt", [0.3, 0.3, 0.4], 3, True, 0, 0.0)