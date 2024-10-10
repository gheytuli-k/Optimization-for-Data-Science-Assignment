import networkx as nx
import gurobipy as gp
from gurobipy import Model, GRB, quicksum


def cycle_based_based_optimization(G: nx.DiGraph, length_limit: int = 5) -> list[list[object]]:
    """
    Solves the kidney exchange problem using a cycle-based approach based on the given input graph.

    :param G: The network of the kidney exchange problem.
    :param length_limit: The maximum length of the cycles to consider.
    :return: A list of cycles that form a solution to the kidney exchange problem based on the given setup.
    """

    # Precomputing all the cycles in the graph upto a certain length
    all_cycles = list(nx.simple_cycles(G))
    all_cycles = [cycle for cycle in all_cycles if 2 <=
                  len(cycle) <= length_limit]
    print(f"Total cycles found (length <= {length_limit}): {len(all_cycles)}")

    model = Model("cycle_based_KEP")  # Initialize the model

    # Creating variables for each cycle and compute their weights using edge weights
    cycle_vars = {}
    cycle_weights = {}
    for i, cycle in enumerate(all_cycles):
        # The weight of the cycles is the sum of the arcs in the cycle
        weight = 0
        cycle_valid = True  # Flag to check if the cycle is valid

        for j in range(len(cycle) - 1):
            u = cycle[j]
            # Ensures wrap-around to form a cycle
            v = cycle[(j + 1) % len(cycle)]
            if G.has_edge(u, v):
                weight += G[u][v]['weight']
            else:
                print(
                    f"Edge from {u} to {v} not found in graph. Skipping cycle.")
                cycle_valid = False
                break  # Skip this cycle if any edge is missing

        if cycle_valid:
            cycle_weights[i] = weight
            cycle_vars[i] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}")

    if not cycle_vars:
        print("No valid cycles found. Exiting optimization.")
        return []

    # Setting the objective so it maximize the total weight of selected cycles
    objective = gp.quicksum(
        cycle_weights[i] * cycle_vars[i] for i in cycle_vars)
    model.setObjective(objective, GRB.MAXIMIZE)

    # Adding constraints to ensure vertex-disjoint cycles
    vertex_cycles = {}
    for i, cycle in enumerate(all_cycles):
        if i not in cycle_vars:
            continue  # Skipping invalid cycles
        for vertex in cycle:
            if vertex not in vertex_cycles:
                vertex_cycles[vertex] = []
            vertex_cycles[vertex].append(i)

    for vertex, cycle_indices in vertex_cycles.items():
        model.addConstr(
            gp.quicksum(cycle_vars[i] for i in cycle_indices) <= 1,
            name=f"vertex_{vertex}_constraint"
        )

    # Optimize the 'GOD' model that I have created and try to save lives :)
    model.optimize()

    # Retrieve and return the optimal solution
    if model.status == GRB.OPTIMAL:
        selected_cycles = [all_cycles[i]
                           for i in cycle_vars if cycle_vars[i].X > 0.5]
        total_weight = sum(cycle_weights[i]
                           for i in cycle_vars if cycle_vars[i].X > 0.5)
        print(f"Total weight of selected cycles: {total_weight}")
        print(f"Selected cycles: {selected_cycles}")
        return selected_cycles

    else:
        print("No optimal solution found.")
        return []
