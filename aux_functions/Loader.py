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
