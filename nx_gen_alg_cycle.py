import networkx as nx
import matplotlib.pyplot as plt
import os
import random
from collections import defaultdict
import numpy as np
from collections import Counter
import time
import pandas as pd
from tqdm import tqdm

from google.colab import drive
drive.mount('/content/drive')

def Loader(file_path: str) -> dict[int, int, int, list[dict], list[dict], list[dict]]:
    """
    Load the data from the input file and return it as a dictionary.
    :param file_path: The path to the input file.
    :return: A dictionary containing the following keys:
        - 'nr_pairs': The number of pairs.
        - 'nr_NDDs': The number of NDDs.
        - 'nr_arcs': The number of arcs.
        - 'pairs_NDDs': A list of dictionaries, each containing the
            following keys:
                - 'pair_id': The ID of the pair.
                - 'is_NDD': 0 if the pair is a pair, 1 if it is an NDD.
                - 'donor_blood_type': The blood type of the donor.
                - 'patient_blood_type': The blood type of the patient.
                - 'patient_vpra': The VPRA of the patient.
        - 'pairs_NoNDDs': A list of dictionaries, with the same keys as above.
        - 'arcs': A list of dictionaries, each containing the following keys:
                - 'donor_pair_id': The ID of the donor pair.
                - 'patient_pair_id': The ID of the patient pair.
                - 'weight': The weight of the arc.
    """

    with open(file_path, "r") as file:
        # Step 1: Read the header
        pair_header = file.readline().strip().split(" ")[-1]
        nr_pairs = int(pair_header)
        print(f"Number of pairs: {nr_pairs}")

        NDD_header = file.readline().strip().split(" ")[-1]
        nr_NDDs = int(NDD_header)
        print(f"Number of NDDs: {nr_NDDs}")

        arc_header = file.readline().strip().split(" ")[-1]
        nr_arcs = int(arc_header)
        print(f"Number of arcs: {nr_arcs}")

        # Step 2: Read the pairs and NDDs
        pairs_NDDs = []
        pairs_NoNDDs = []
        for _ in range(nr_pairs):
            line = file.readline().strip().split(",")
            pair_id = int(line[0])
            is_NDD = int(line[1])  # 0 for pair, 1 for NDD
            # 0 for O, 1 for A, 2 for B, 3 for AB
            donor_blood_type = int(line[2])
            # 0 for O, 1 for A, 2 for B, 3 for AB
            patient_blood_type = int(line[3])
            # 0 for <0.5, 1 for 0.5-0.85, 2 for >0.85
            patient_vpra = int(line[4])

            # store each pair/NDD as a dictionary
            if is_NDD == 0:
                pairs_NoNDDs.append({
                    "pair_id": pair_id,
                    "is_NDD": is_NDD,
                    "donor_blood_type": donor_blood_type,
                    "patient_blood_type": patient_blood_type,
                    "patient_vpra": patient_vpra})
            else:
                pairs_NDDs.append({
                    "pair_id": pair_id,
                    "is_NDD": is_NDD,
                    "donor_blood_type": donor_blood_type,
                    "patient_blood_type": patient_blood_type,
                    "patient_vpra": patient_vpra})

        # Step 3: Read the arcs
        arcs = []
        for _ in range(nr_arcs):
            arc_line = file.readline().strip().split(',')
            # Strip the opening parenthesis
            donor_pair_id = int(arc_line[0][1:])
            patient_pair_id = int(arc_line[1].split(')')[
                                  0])  # Strip the closing parenthesis
            weight = int(arc_line[2])

            arcs.append({
                'donor_pair_id': donor_pair_id,
                'patient_pair_id': patient_pair_id,
                'weight': weight
            })

        # close the file
        file.close()

    # Step 4: Return the parsed data
    return {
        'nr_pairs': nr_pairs,
        'nr_NDDs': nr_NDDs,
        'nr_arcs': nr_arcs,
        'pairs_NDDs': pairs_NDDs,
        'pairs_NoNDDs': pairs_NoNDDs,
        'arcs': arcs
    }

def data_to_adjacency(data: dict[int, int, int, list[dict], list[dict], list[dict]], god_donor: bool=False, god_donor_weight: int=0) -> dict[int, list[dict]]:
    """
    Turns the input data into a adjacency representation of the graph
    :param data: dictionary with the following keys
        - 'nr_pairs': The number of pairs.
        - 'nr_NDDs': The number of NDDs.
        - 'nr_arcs': The number of arcs.
        - 'pairs_NDDs': A list of dictionaries, each containing the
            following keys:
                - 'pair_id': The ID of the pair.
                - 'is_NDD': 0 if the pair is a pair, 1 if it is an NDD.
                - 'donor_blood_type': The blood type of the donor.
                - 'patient_blood_type': The blood type of the patient.
                - 'patient_vpra': The VPRA of the patient.
        - 'pairs_NoNDDs': A list of dictionaries, with the same keys as above.
        - 'arcs': A list of dictionaries, each containing the following keys:
                - 'donor_pair_id': The ID of the donor pair.
                - 'patient_pair_id': The ID of the patient pair.
                - 'weight': The weight of the arc.

    :param god_donor: boolean to indicate if all the noNDD pairs should be connected to god donors
    :param god_donor_weight: weight of the edge between god donor and noNDD pairs when an edge is not present
    :return: adjacency representation of the graph
    """

    graph = {i: {} for i in range(data["nr_pairs"]+data["nr_NDDs"])}


    for arc in data["arcs"]:
        donor_pair_id = arc["donor_pair_id"]
        patient_pair_id = arc["patient_pair_id"]
        weight = arc["weight"]
        graph[donor_pair_id][patient_pair_id] = {'weight': weight}

    NDD_pairs = pair_ids = [pair['pair_id'] for pair in data['pairs_NDDs']]

    if god_donor:
        for i in range(data["nr_pairs"]):
            for j in NDD_pairs:
                if j not in graph[i]:
                    graph[i][j] = {'weight': god_donor_weight}

    return graph

def load_data(data_path):
    data = Loader(data_path)

    # Dict with pair_id as key and patient_vpra as value
    pair_vpra = {}
    for pair in data['pairs_NDDs']:
        pair_vpra[pair['pair_id']] = pair['patient_vpra']
    for pair in data['pairs_NoNDDs']:
        pair_vpra[pair['pair_id']] = pair['patient_vpra']
    data['pair_vpra'] = pair_vpra

    return data

def create_graph(data, god_donor=True, god_donor_weight=0):
    # Making adjacency out of the data and connecting all the pairs to god donor
    adjacency = data_to_adjacency(data, god_donor, god_donor_weight)
    NDDs = [pair['pair_id'] for pair in data['pairs_NDDs']]
    G = nx.DiGraph(adjacency)

    # Connect all the pairs to god donor with weight god_donor_weight
    all_nodes = set(G.nodes())
    donors = all_nodes - set(NDDs)
    for donor in donors:
        for ndd in NDDs:
            if not G.has_edge(donor, ndd):
                G.add_edge(donor, ndd, weight=god_donor_weight)
    return G, NDDs

# DFS-based function to find all cycles of length <= k
def find_cycles_of_length_leq_k(G, k):
    cycles = nx.simple_cycles(G, k)

    # Using a set for efficient cycle storage and uniqueness
    cycle_set = set()
    for cycle in cycles:
        min_node = min(cycle)
        min_index = cycle.index(min_node)
        canonical_cycle = tuple(cycle[min_index: ] + cycle[: min_index])  # Convert to canonical form
        cycle_set.add(canonical_cycle)  # Add the canonical cycle to the set
    return list(cycle_set)

# Function to randomly initialize the population matrix
def initialize_population(num_paths, population_size, path_dct, path_keys):
    # Initialize the population matrix with zeros
    population = np.zeros((population_size, num_paths), dtype=int)
    # Initialize the fitness result with zeros
    fitness_result = np.zeros(population_size)

    # Assign random edges for each individual
    for i in tqdm(range(population_size)):
        # Only add valid solution
        while True:
            checking = True
            option = np.zeros(num_paths)
            size = np.random.randint(low = 1, high = 6, size = 1)[0]
            for j in np.random.randint(low = 0, high = num_paths, size = size):
                option[j] = 1
            selected_paths = find_valid_path(option, path_keys)
            if selected_paths:
                break
        total_weight = sum(path_dct[path] for path in selected_paths)
        population[i] = option
        fitness_result[i] = total_weight

    return population, fitness_result

def find_valid_path(option, path_keys):
    # The array should not have just 0 or just 1 (If there are many possible arrays)
    if len(np.unique(option)) == 1 and len(option) > 5:
        return []
    # If any overlapping exists, return empty list immediately
    seen_node = set()
    selected_paths = []
    for i in range(len(option)):
        if option[i] == 1:
            cycle = path_keys[i]
            for node in cycle:
                if node in seen_node:
                    return []
                else:
                    seen_node.add(node)
            selected_paths.append(cycle)
    return selected_paths

def selection(population, fitness_result, path_dct, path_keys):
    # Roulette wheel
    # Generate a random number x between 0 and 1
    # Check it in cumulative normalized fitness array
    # Find the index i such that array[i - 1] < x < array[i]
    selection_probs = fitness_result / sum(fitness_result)
    non_zero_prob_count = len(np.where(selection_probs != 0)[0])
    sorted_probs = -np.sort(-selection_probs)[: non_zero_prob_count]
    cum_probs = np.cumsum(sorted_probs)
    corresponding_solution = np.argsort(-selection_probs)[: non_zero_prob_count]

    # Create 2 parents for crossover and mutation
    parent_1, parent_2 = [population[corresponding_solution[np.searchsorted(cum_probs, random.uniform(0, 1))]] for _ in range(2)]
    child_lst = []

    crossover_lst = []
    seen_idx = set()
    patience = 0
    while True:
        # Crossover
        # Create a random index to split
        split_idx = np.random.randint(low = 1, high = len(parent_1) - 1, size = 1)[0]
        # Skip seen index that does not work
        if split_idx in seen_idx:
            continue
        seen_idx.add(split_idx)

        # Combine different parts from 2 parents to create new children
        child_1 = np.concatenate((parent_1[: split_idx], parent_2[split_idx: ]))
        child_2 = np.concatenate((parent_2[: split_idx], parent_1[split_idx: ]))
        temp_child = [child_1, child_2]

        # Only use valid child that is different from its parents
        for child in temp_child:
            if not(all(child == parent_1) or all(child == parent_2)):
                selected_paths = find_valid_path(child, path_keys)
                if selected_paths:
                    crossover_lst.append(selected_paths)
                    child_lst.append(child)
        if len(crossover_lst) > 0 or len(seen_idx) == len(parent_1) - 2 or patience == 500:
            break
        patience += 1

    mutation_lst = []
    attempt = 0
    while True:
        # Mutation
        # Create a random array with the same shape as original_array
        random_flips = np.random.rand(*parent_1.shape)
        child_3 = parent_1.copy()
        child_4 = parent_2.copy()
        # Flip probability is 1 / number of cycles
        prob = 1 / len(parent_1)
        # Flip the elements where random_flips is less than flip_probability
        child_3[random_flips < prob] = 1 - child_3[random_flips < prob]
        child_4[random_flips < prob] = 1 - child_4[random_flips < prob]
        temp_child = [child_3, child_4]

        # Only use valid child that is different from its parents
        for child in temp_child:
            if not(all(child == parent_1) or all(child == parent_2)):
                selected_paths = find_valid_path(child, path_keys)
                if selected_paths:
                    child_lst.append(child)
                    mutation_lst.append(selected_paths)
        if len(mutation_lst) > 0 or attempt == 20:
            break
        attempt += 1

    # Evaluate the child's fitness
    child_combo = crossover_lst + mutation_lst
    child_fitness = [sum(path_dct[path] for path in child) for child in child_combo]

    # Append the new children's result
    if child_combo:
        for fit, combo in zip(child_fitness, child_lst):
            if fit > np.min(fitness_result):
                replace_idx = np.argmin(fitness_result)
                population[replace_idx] = combo
                fitness_result[replace_idx] = fit
    return population, fitness_result

def genetic_algorithm(file_path, max_length = 5, population_size = 10000, constant_solution_nr = 20, beta = 0):
    fit_start = time.time()
    # Load the data
    data = load_data(file_path)
    G, NDDs = create_graph(data)

    # Find vpra of every pair and store it into a dictionary
    pair_vpra = data['pair_vpra']

    # Find cycles of length <= max_cycle_length
    cycles = find_cycles_of_length_leq_k(G, max_length)
    cycle_nr = len(cycles)

    # Build cycle/chain-weight dictionary, where the weight is adjusted to account for vpra pair
    path_dct = {}
    real_path_dct = {}
    for cycle in cycles:
        total_weight = 0
        total_real_weight = 0
        cycle_length = len(cycle)
        for i in range(cycle_length):
            u = cycle[i]
            v = cycle[(i + 1) % cycle_length]  # Next node in the cycle, wrapping around
            total_weight += G[u][v]["weight"] * (1 + beta * pair_vpra.get(v, 0))
            total_real_weight += G[u][v]["weight"]
        path_dct[tuple(cycle)] = total_weight
        real_path_dct[tuple(cycle)] = total_real_weight

    path_keys = list(path_dct.keys())

    if cycle_nr == 1:
        final_time = time.time() - fit_start
        return file_path[: -4], real_path_dct[path_keys[0]], final_time, final_time, path_dct[path_keys[0]], max_length, len(path_keys[0]), [path_keys[0]], beta, population_size

    if cycle_nr == 2:
        seen_nodes = set()
        for key in path_keys:
            for node in key:
                if node in seen_nodes:
                    best_cycle = max(path_dct, key = path_dct.get)
                    final_time = time.time() - fit_start
                    return file_path[: -4], real_path_dct[best_cycle], final_time, final_time, path_dct[best_cycle], max_length, len(best_cycle), [best_cycle], beta, population_size
                seen_nodes.add(node)
        final_time = time.time() - fit_start
        return file_path[: -4], sum(real_path_dct[i] for i in path_keys), final_time, final_time, sum(path_dct[i] for i in path_keys), max_length, sum(len(i) for i in path_keys), path_keys, beta, population_size

    # Initialize the population and fitness result given a population size
    # All solutions must be valid
    population, fitness_result = initialize_population(len(path_keys), population_size, path_dct, path_keys)

    fit_time = time.time() - fit_start

    gen_start = time.time()

    # Initialize variables
    max_fitness = float('-inf')
    unchange_nr = 0
    iter_nr = 0

    # Continue generating new solution until the best solution remains the same
    while unchange_nr < constant_solution_nr:
        population, fitness_result = selection(population, fitness_result, path_dct, path_keys)
        new_fitness = np.max(fitness_result)
        if max_fitness == new_fitness:
            unchange_nr += 1
        else:
            unchange_nr = 0
            max_fitness = new_fitness
        iter_nr += 1
        print(f"Iteration {iter_nr} done")

    # Return the paths
    all_idxs = np.where(fitness_result == max_fitness)[0]
    result = []
    real_weight_lst = []
    for idx in all_idxs:
        cycle_sol = []
        temp_result = population[idx]
        real_weight = 0
        for i in range(len(temp_result)):
            if temp_result[i] == 1:
                real_weight += real_path_dct[path_keys[i]]
                cycle_sol.append(path_keys[i])
        if cycle_sol not in result:
            result.append(cycle_sol)
            real_weight_lst.append(real_weight)

    # If there exists solutions with the same adjusted weights, select one with higher actual weight
    max_weight = np.max(real_weight_lst)
    eqv_idx = np.argmax(real_weight_lst)
    final_result = result[eqv_idx]
    exchange_nr = sum(len(i) for i in final_result)
    gen_time = time.time() - gen_start

    # Dataset name, original weight, time to generate cycles, total time, adjusted weight, max cycle length, exchange number, cycle, beta, population
    return file_path[: -4], max_weight, fit_time, fit_time + gen_time, max_fitness, max_length, exchange_nr, final_result, beta, population_size

running_files = [
    'Delorme_50_NDD_Unit_0.txt',
    'Delorme_50_NDD_Unit_1.txt',
    'Delorme_50_NDD_Unit_2.txt',
    'Delorme_50_NDD_Unit_3.txt',
    'Delorme_50_NDD_Unit_4.txt',
    'Delorme_50_NDD_Weight_0.txt',
    'Delorme_50_NDD_Weight_1.txt',
    'Delorme_50_NDD_Weight_2.txt',
    'Delorme_50_NDD_Weight_3.txt',
    'Delorme_50_NDD_Weight_4.txt',
    'Delorme_50_NoNDD_Unit_0.txt',
    'Delorme_50_NoNDD_Unit_1.txt',
    'Delorme_50_NoNDD_Unit_2.txt',
    'Delorme_50_NoNDD_Unit_3.txt',
    'Delorme_50_NoNDD_Unit_4.txt',
    'Delorme_50_NoNDD_Weight_0.txt',
    'Delorme_50_NoNDD_Weight_1.txt',
    'Delorme_50_NoNDD_Weight_2.txt',
    'Delorme_50_NoNDD_Weight_3.txt',
    'Delorme_50_NoNDD_Weight_4.txt',
    'Saidman_50_NDD_Unit_0.txt',
    'Saidman_50_NDD_Unit_1.txt',
    'Saidman_50_NDD_Unit_2.txt',
    'Saidman_50_NDD_Unit_3.txt',
    'Saidman_50_NDD_Unit_4.txt',
    'Saidman_50_NDD_Weight_0.txt',
    'Saidman_50_NDD_Weight_1.txt',
    'Saidman_50_NDD_Weight_2.txt',
    'Saidman_50_NDD_Weight_3.txt',
    'Saidman_50_NDD_Weight_4.txt',
    'Saidman_50_NoNDD_Unit_0.txt',
    'Saidman_50_NoNDD_Unit_1.txt',
    'Saidman_50_NoNDD_Unit_2.txt',
    'Saidman_50_NoNDD_Unit_3.txt',
    'Saidman_50_NoNDD_Unit_4.txt',
    'Saidman_50_NoNDD_Weight_0.txt',
    'Saidman_50_NoNDD_Weight_1.txt',
    'Saidman_50_NoNDD_Weight_2.txt',
    'Saidman_50_NoNDD_Weight_3.txt',
    'Saidman_50_NoNDD_Weight_4.txt',
    'Delorme_200_NDD_Unit_0.txt',
    'Delorme_200_NDD_Unit_1.txt',
    'Delorme_200_NDD_Unit_2.txt',
    'Delorme_200_NDD_Unit_3.txt',
    'Delorme_200_NDD_Unit_4.txt',
    'Delorme_200_NDD_Weight_0.txt',
    'Delorme_200_NDD_Weight_1.txt',
    'Delorme_200_NDD_Weight_2.txt',
    'Delorme_200_NDD_Weight_3.txt',
    'Delorme_200_NDD_Weight_4.txt',
    'Delorme_200_NoNDD_Unit_0.txt',
    'Delorme_200_NoNDD_Unit_1.txt',
    'Delorme_200_NoNDD_Unit_2.txt',
    'Delorme_200_NoNDD_Unit_3.txt',
    'Delorme_200_NoNDD_Unit_4.txt',
    'Delorme_200_NoNDD_Weight_0.txt',
    'Delorme_200_NoNDD_Weight_1.txt',
    'Delorme_200_NoNDD_Weight_2.txt',
    'Delorme_200_NoNDD_Weight_3.txt',
    'Delorme_200_NoNDD_Weight_4.txt',
    'Saidman_200_NDD_Unit_0.txt',
    'Saidman_200_NDD_Unit_1.txt',
    'Saidman_200_NDD_Unit_2.txt',
    'Saidman_200_NDD_Unit_3.txt',
    'Saidman_200_NDD_Unit_4.txt',
    'Saidman_200_NDD_Weight_0.txt',
    'Saidman_200_NDD_Weight_1.txt',
    'Saidman_200_NDD_Weight_2.txt',
    'Saidman_200_NDD_Weight_3.txt',
    'Saidman_200_NDD_Weight_4.txt',
    'Saidman_200_NoNDD_Unit_0.txt',
    'Saidman_200_NoNDD_Unit_1.txt',
    'Saidman_200_NoNDD_Unit_2.txt',
    'Saidman_200_NoNDD_Unit_3.txt',
    'Saidman_200_NoNDD_Unit_4.txt',
    'Saidman_200_NoNDD_Weight_0.txt',
    'Saidman_200_NoNDD_Weight_1.txt',
    'Saidman_200_NoNDD_Weight_2.txt',
    'Saidman_200_NoNDD_Weight_3.txt',
    'Saidman_200_NoNDD_Weight_4.txt',
]

sparse_files = [
    'RandomSparse_50_NoNDD_Unit_0.txt',
    'RandomSparse_50_NoNDD_Unit_1.txt',
    'RandomSparse_50_NoNDD_Unit_2.txt',
    'RandomSparse_50_NoNDD_Unit_3.txt',
    'RandomSparse_50_NoNDD_Unit_4.txt',
    'RandomSparse_50_NoNDD_Weight_0.txt',
    'RandomSparse_50_NoNDD_Weight_1.txt',
    'RandomSparse_50_NoNDD_Weight_2.txt',
    'RandomSparse_50_NoNDD_Weight_3.txt',
    'RandomSparse_50_NoNDD_Weight_4.txt',
    'RandomSparse_200_NDD_Unit_0.txt',
    'RandomSparse_200_NDD_Unit_1.txt',
    'RandomSparse_200_NDD_Unit_2.txt',
    'RandomSparse_200_NDD_Unit_3.txt',
    'RandomSparse_200_NDD_Unit_4.txt',
    'RandomSparse_200_NDD_Weight_0.txt',
    'RandomSparse_200_NDD_Weight_1.txt',
    'RandomSparse_200_NDD_Weight_2.txt',
    'RandomSparse_200_NDD_Weight_3.txt',
    'RandomSparse_200_NDD_Weight_4.txt',
    'RandomSparse_200_NoNDD_Unit_0.txt',
    'RandomSparse_200_NoNDD_Unit_1.txt',
    'RandomSparse_200_NoNDD_Unit_2.txt',
    'RandomSparse_200_NoNDD_Unit_3.txt',
    'RandomSparse_200_NoNDD_Unit_4.txt',
    'RandomSparse_200_NoNDD_Weight_0.txt',
    'RandomSparse_200_NoNDD_Weight_1.txt',
    'RandomSparse_200_NoNDD_Weight_2.txt',
    'RandomSparse_200_NoNDD_Weight_3.txt',
    'RandomSparse_200_NoNDD_Weight_4.txt',
    'RandomSparse_500_NDD_Unit_0.txt',
    'RandomSparse_500_NDD_Unit_1.txt',
    'RandomSparse_500_NDD_Unit_2.txt',
    'RandomSparse_500_NDD_Unit_3.txt',
    'RandomSparse_500_NDD_Unit_4.txt',
    'RandomSparse_500_NDD_Weight_0.txt',
    'RandomSparse_500_NDD_Weight_1.txt',
    'RandomSparse_500_NDD_Weight_2.txt',
    'RandomSparse_500_NDD_Weight_3.txt',
    'RandomSparse_500_NDD_Weight_4.txt',
    'RandomSparse_500_NoNDD_Unit_0.txt',
    'RandomSparse_500_NoNDD_Unit_1.txt',
    'RandomSparse_500_NoNDD_Unit_2.txt',
    'RandomSparse_500_NoNDD_Unit_3.txt',
    'RandomSparse_500_NoNDD_Unit_4.txt',
    'RandomSparse_500_NoNDD_Weight_0.txt',
    'RandomSparse_500_NoNDD_Weight_1.txt',
    'RandomSparse_500_NoNDD_Weight_2.txt',
    'RandomSparse_500_NoNDD_Weight_3.txt',
    'RandomSparse_500_NoNDD_Weight_4.txt'
]

beta_lst = [0, 0.5, 1]
population_size = 1000
max_length_lst = [3, 4, 5]

config_lst = []
for file_path in running_files:
    for beta in beta_lst:
        for max_length in max_length_lst:
            config_lst.append((file_path, beta, population_size, max_length))

# Experiment Loop
# 640

starting_idx = 719
col_names = ["Dataset", "Total Weight", "Cycle Time (s)", "Total Time (s)", "Total Model Weight", "Max Cycle Length", "Exchange Number", "Selected Cycles", "Beta", "Population"]
os.chdir('/content/drive/MyDrive/Optimization/data/Instance Files')

for idx, (file_path, beta, population_size, max_length) in enumerate(config_lst[starting_idx: ]):
    print(f"Experiment {idx + starting_idx} starts!")
    if max_length == 5 or max_length == 4:
        continue
    res = genetic_algorithm(file_path = file_path, max_length = max_length, population_size = population_size, beta = beta)
    try:
        final_df = pd.read_excel("/content/drive/MyDrive/Optimization/Genetic Result/genetic_result.xlsx")
        df = pd.DataFrame([res], columns = col_names)
        final_df = pd.concat([final_df, df])
    except:
        final_df = pd.DataFrame([res], columns = col_names)

    final_df.to_excel("/content/drive/MyDrive/Optimization/Genetic Result/genetic_result.xlsx", index = False)

