#imports and functions from aux_functions
from aux_functions import Loader
from aux_functions import data_to_adjacency
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import time as time
import os
from pathlib import Path
import csv

def load_data(data_path):
    """
    Loads data from the specified path and processes it to create a dictionary and also maps pair IDs to patient VPRA values.
    Args:
        data_path (str): The path to the data file.
    Returns:
        dict: A dictionary containing the loaded data with an additional key 
              'pair_vpra' which maps pair IDs to patient VPRA values.
    """
    
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
    """
    Creates a directed graph from the given data, optionally connecting all pairs to a 'god donor'.
    Args:
        data (dict): The input data containing pairs and non-directed donors (NDDs).
        god_donor (bool, optional): If True, connects all pairs to a 'god donor'. Defaults to True.
        god_donor_weight (int, optional): The weight of the edges connecting pairs to the 'god donor'. Defaults to 0.
    Returns:
        networkx.DiGraph: A directed graph representing the connections between pairs and NDDs.
    """

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
    """
    Finds all simple cycles in a directed graph G with a length less than or equal to max_length.
    Parameters:
    G (networkx.DiGraph): A directed graph in which to find cycles.
    max_length (int): The maximum length of cycles to be considered.
    Returns:
    tuple: A tuple containing:
        - cycles (list of lists): A list of cycles, where each cycle is represented as a list of nodes.
        - cyc_time (float): The time taken to find the cycles in seconds.
    """

    cts = time.time()
    cycles = []

    #only generate cycle of length <= max_length
    for cycle in nx.simple_cycles(G, max_length):
        if len(cycle) <= max_length:
            cycles.append(cycle)
    cyc_time = time.time() - cts    
    return cycles, cyc_time

def optimize_sequences(sequences, G, beta, pair_vpra):
    """
    Optimize sequences in a kidney exchange problem using Gurobi.
    Parameters:
    sequences (list of lists): A list of sequences (cycles), where each sequence is a list of nodes.
    G (networkx.Graph): A directed graph where nodes represent patients/donors and edges represent possible exchanges.
    beta (float): A factor to adjust the weight of edges based on the pair_vpra values.
    pair_vpra (dict): A dictionary where keys are node identifiers and values are VPRA (Virtual PRA) values.
    Returns:
    selected_sequences (list of lists): A list of selected sequences (cycles) after optimization.
    total_model_weight (float): The total weight of the selected sequences.
    """

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
    total_model_weight = sum(sequence_weights[i] for i in sequence_vars if sequence_vars[i].X > 0.5)
    print(f"Total weight of selected sequences (model): {total_model_weight}")
    return selected_sequences, total_model_weight

def output_results(selected_sequences, G, max_cycle_length, cycle_time, start_time, experiment, total_model_weight):
    """
    Outputs the results of the cycle selection process and optionally returns performance metrics.
    Parameters:
    selected_sequences (list of lists): A list of sequences, where each sequence is a list of nodes representing a cycle.
    G (networkx.Graph): The graph containing nodes and weighted edges.
    max_cycle_length (int): The maximum allowed length of a cycle.
    cycle_time (float): The time taken to generate cycles.
    start_time (float): The start time of the entire process.
    experiment (bool): Flag indicating whether to return performance metrics.
    total_model_weight (float): The total weight of the model used in the experiment.
    Returns:
    tuple: If experiment is True, returns a tuple containing total_weight, cycle_time, total_time, and total_model_weight.
    """

    print("Set max cycle Length:", max_cycle_length)
    if selected_sequences:
        print("Selected Cycles and their weights:")
        total_weight = 0
        for seq in selected_sequences:
            weight = 0
            # Sum weights of consecutive edges in the sequence
            for j in range(len(seq) - 1):
                u = seq[j]
                v = seq[j + 1]
                if G.has_edge(u, v):
                    weight += G[u][v]['weight']
            # Add the weight of the last edge if it's a cycle
            if seq[0] == seq[-1] and len(seq) > 1:
                u = seq[-1]
                v = seq[0]
                if G.has_edge(u, v):
                    weight += G[u][v]['weight']
            total_weight += weight
            print(f"Cycle {seq} has weight: {weight}")  # Per cycle weight
        print(f"Total weight of selected cycles: {total_weight}")  # Total weight (no beta or vpra factor)
    else:
        print("No sequences were selected.")
    print(f"Cycle generation time: {cycle_time:.2f} seconds")
    total_time = time.time() - start_time
    print(f"Total time taken: {total_time:.2f} seconds")
    if experiment:
        return total_weight, cycle_time, total_time, total_model_weight


def cycle_based_based_optimisation(data_path, max_cycle_length: int = 3, god_donor: bool = True, god_donor_weight: int = 10, beta: float = 0.0, experiment = False):
    """
    Perform cycle-based optimization on a given dataset.
    Parameters:
    data_path (str): Path to the input data file.
    max_cycle_length (int, optional): Maximum length of cycles to consider. Default is 3.
    god_donor (bool, optional): Whether to include a god donor in the graph. Default is True.
    god_donor_weight (int, optional): Weight assigned to the god donor. Default is 0.
    beta (float, optional): Regularization parameter for optimization. Default is 0.0.
    experiment (bool, optional): Flag to indicate if the function is run in experimental mode. Default is False.
    Returns:
    tuple: A tuple containing:
        - total_weight (float): The total weight of the selected sequences.
        - cycle_time (float): The time taken to find cycles.
        - total_time (float): The total time taken for the optimization process.
        - total_model_weight (float): The total model weight after optimization.
    """
    
    start_time = time.time()
    
    #load data
    data = load_data(data_path)
    
    #create graph
    G = create_graph(data, god_donor, god_donor_weight)

    #find cycles
    cycles, cycle_time = find_cycles(G, max_cycle_length)

    #optimize
    selected_sequences, total_weight = optimize_sequences(cycles, G, beta, data['pair_vpra'])
    
    #output results
    total_weight, cycle_time, total_time, total_model_weight = output_results(selected_sequences, G, max_cycle_length, cycle_time, start_time, experiment, total_weight)

    return  total_weight, cycle_time, total_time, total_model_weight


#cycle_based_based_optimisation("Instance Files\Saidman_50_NDD_Unit_0.txt", 3, True, 0, 100, False)


def run_experiments(dataset_folder, max_cycle_length: int = 3, god_donor: bool = True, god_donor_weight: int = 0, beta: float = 0.0):
    """
    Run optimization experiments on all dataset files in the specified folder.
    Parameters:
    dataset_folder (str): Path to the folder containing dataset files.
    max_cycle_length (int, optional): Maximum cycle length for the optimization. Default is 3.
    god_donor (bool, optional): Flag to indicate if a god donor should be used. Default is True.
    god_donor_weight (int, optional): Weight of the god donor. Default is 0.
    beta (float, optional): Beta parameter for the optimization. Default is 0.0.
    Returns:
    dict: A dictionary containing the results of the optimization for each dataset file.
          The keys are the filenames, and the values are dictionaries with the following keys:
          - "Total Weight": The total weight of the optimized solution.
          - "Cycle Time (s)": The time taken to find the optimal cycle.
          - "Total Time (s)": The total time taken for the optimization.
          - "Total Model Weight": The total weight of the model used in the optimization.
    """

    
    results_dict = {}

    dataset_path = Path(dataset_folder)

    for file_path in dataset_path.glob("*.txt"): #loop through all files in the folder

        try:
            # run optimization on the current file
            metrics = cycle_based_based_optimisation(
                data_path=str(file_path),
                max_cycle_length=max_cycle_length,
                god_donor=god_donor,
                god_donor_weight=god_donor_weight,
                beta=beta,
                experiment=True
            )
        
            if metrics: #if metrics are returned, add them to the results dictionary
                total_weight, cycle_time, total_time, total_model_weight = metrics
                
                results_dict[file_path.name[:-4]] = {
                    "Total Weight": total_weight,
                    "Cycle Time (s)": cycle_time,
                    "Total Time (s)": total_time,
                    "Total Model Weight": total_model_weight,
                    "Max Cycle Length": max_cycle_length

                }

        except Exception as e:
            pass

    
    return results_dict, 
        
def dict_to_csv(results_dict, output_file):
    """
    Write the results of the optimization experiments to a CSV file.
    Parameters:
    results_dict (dict): A dictionary containing the results of the optimization experiments.
    output_file (str): The path to the output CSV file.
    """

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Dataset", "Total Weight", "Cycle Time (s)", "Total Time (s)", "Total Model Weight", "Max Cycle Length"])
        for dataset, metrics in results_dict.items():
            writer.writerow([dataset, metrics["Total Weight"], metrics["Cycle Time (s)"], metrics["Total Time (s)"], metrics["Total Model Weight"], metrics["Max Cycle Length"]])

e = run_experiments("Instance Files Test", 3, True, 0, 0)
dict_to_csv(e[0], "results.csv")
