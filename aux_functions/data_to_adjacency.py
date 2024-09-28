def data_to_adjacency(data: dict[int, int, int, list[dict], list[dict], list[dict]], god_donor: bool=False, god_donor_weight: int=0) -> dict[int, list[dict]]:
    """
    Turns the input data into a adjacency representation of the graph
    :paramm data: dictionary with the following keys
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
        graph[donor_pair_id][patient_pair_id] = weight
    
    NDD_pairs = pair_ids = [pair['pair_id'] for pair in data['pairs_NDDs']]

    if god_donor:
        for i in range(data["nr_pairs"]):
            for j in NDD_pairs:
                if j not in graph[i]:
                    graph[i][j] = god_donor_weight
            
    return graph
    