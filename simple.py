from aux_functions import Loader
from aux_functions import draw_graph
from aux_functions import data_to_adjacency
from gurobipy import Model, GRB
import networkx as nx
import matplotlib.pyplot as plt

def draw_kidney_exchange(matches):
    """
    Draws a directed graph for kidney exchange using networkx.
    
    :param matches: A dictionary representing donor-recipient matches
                    (e.g., {donor: recipient})
    """
    # Create a directed graph
    G = nx.DiGraph()

    # Add edges to the graph based on the matches
    for donor, recipient in matches.items():
        G.add_edge(donor, recipient)

    # Generate a layout for our graph
    pos = nx.spring_layout(G)

    # Draw the graph
    plt.figure(figsize=(10, 8))
    
    # Draw the nodes
    nx.draw_networkx_nodes(G, pos, node_size=700, node_color='skyblue', alpha=0.9)
    
    # Draw the edges
    nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=20, edge_color='black', width=2)

    # Draw labels for nodes
    nx.draw_networkx_labels(G, pos, font_size=12, font_color='black')

    # Add edge labels (optional)
    edge_labels = {(donor, recipient): f"{donor} -> {recipient}" for donor, recipient in matches.items()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red', font_size=10)

    # Display the plot
    plt.title("Kidney Exchange Matchings")
    plt.axis('off')  # Turn off the axis
    plt.show()

# data = Loader("Instance Files\Delorme_200_NDD_Unit_4.txt")

# model = Model("Kidney Exchange")

# # x = {}
# # for donor in adjacency:
# #     for recipient in adjacency[donor]:
# #         x[donor, recipient] = model.addVar(vtype=GRB.BINARY, name=f"x_{donor}_{recipient}")

# # # maximize total weight of donation
# # model.setObjective(sum(adjacency[donor][recipient] * x[donor, recipient] for donor in adjacency for recipient in adjacency[donor]), GRB.MAXIMIZE)

# # # Add constraints 
# # # Each donor can donate to at most one recipient
# # for donor in adjacency:
# #     model.addConstr(sum(x[donor, recipient] for recipient in adjacency[donor]) <= 1, name=f"donor_{donor}")

# # # Each recipient can receive from at most one donor
# # recipients = set()
# # for donor in adjacency:
# #     recipients.update(adjacency[donor].keys())

# # for recipient in recipients:
# #     model.addConstr(sum(x[donor, recipient] for donor in adjacency if recipient in adjacency[donor]) <= 1, name=f"recipient_{recipient}")

# # model.optimize()


# # transplant = {}
# # if model.status == GRB.OPTIMAL:
# #     print("Optimal Solution Found")
# #     for donor in adjacency:
# #         for recipient in adjacency[donor]:
# #             if x[donor, recipient].x > 0.5:
# #                 transplant[donor] = recipient
# #                 print(f"Match: {donor} -> {recipient}")

# # draw_kidney_exchange(transplant)

# def optimize_kidney_exchange(adjacency, cycle_limit=4):
#     # Create a directed graph from the matches (adjacency list)
#     G = nx.DiGraph(adjacency)

#     # step 1: find all simples cycles of length <= cycle_limit
#     all_cycles = [cycle for cycle in nx.simple_cycles(G) if len(cycle) < cycle_limit]

#     # Initilize the model
#     model = Model("Kidney Exchange")

#     # step 2: create a binary variable for each cycle
#     cycle_vars = {}
#     for idx, cycle in enumerate(all_cycles):
#         cycle_vars[idx] = model.addVar(vtype=GRB.BINARY, name=f"cycle_{idx}")

#     # Step 3: objective function - maximize the number of cycles (or weighted c ycles if needede)
#     model.setObjective(sum(cycle_vars[idx] for idx in cycle_vars), GRB.MAXIMIZE)

#     # Step 4: Add constraints - ensure that each donor or recipient appears in at most one cycle
#     # This ensures each donor-recipient pair is in at most one cycles(simple cycle)
#     for node in G.nodes:
#         model.addConstr(sum(cycle_vars[idx] for idx, cycle in enumerate(all_cycles) if node in cycle) <= 1,
#         name=f"node_{node}_constraint")

#     # Step 5: Optimize the model
#     model.optimize()

#     # Step 6: Extract the cycles from the model
#     if model.status == GRB.OPTIMAL:
#         print("\nOptimal solution found:")
#         for idx, cycle in enumerate(all_cycles):
#             if cycle_vars[idx].x > 0.5:  # If the cycle is selected
#                 print(f"Cycle: {' -> '.join(map(str, cycle))}")
#         else:
#             print("No optimal solution found.")

# optimize_kidney_exchange(adjacency, cycle_limit=4)

data = Loader("Instance Files\Delorme_50_NDD_Unit_2.txt")
adjacency = data_to_adjacency(data, False, 0)

G = nx.DiGraph(adjacency)

# print("Nodes:", G.nodes)
# print("Edges:", G.edges)

def find_chains(G, max_length=4):
    chains = []

    for node in G.nodes:
        stack = [(node, [node])]

        while stack:
            current_node, path = stack.pop()

            if len(path) < max_length:
                for neighbor in G.successors(current_node):
                    if neighbor not in path:
                        stack.append((neighbor, path + [neighbor]))
                    else:
                        chains.append(path)
    
    return chains

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

max_chain_length = 4
max_cycle_length = 4

chains = find_chains(G, max_chain_length)
cycles = find_cycles(G, max_cycle_length)

# print("Chains:", chains)
# print("Cycles:", cycles)

def setup_gurobi_model(cycles):
    # Create a Gurobi model
    model = Model("Maximize_Kidney_Donations")

    # Convert cycles to a set of tuples to remove duplicates
    unique_cycles = list(set(map(tuple, cycles)))  # Remove duplicate cycles

    # Create binary variables for each unique cycle
    cycle_vars = model.addVars(unique_cycles, vtype=GRB.BINARY, name="cycle")

    # Define the objective function
    model.setObjective(
        cycle_vars.sum('*'),  # Maximize total donations from cycles
        GRB.MAXIMIZE
    )

    # Add constraints to ensure no node is included in more than one cycle
    all_nodes = set(node for cycle in unique_cycles for node in cycle)
    for node in all_nodes:
        model.addConstr(
            sum(cycle_vars[c] for c in unique_cycles if node in c) <= 1,
            f"Node_Constraint_{node}"
        )

    return model


# Example usage
gurobi_model = setup_gurobi_model(cycles)

# Print the model to verify
gurobi_model.write("kidney_donations.lp")

gurobi_model.optimize()

# Extract the optimal solution
if gurobi_model.status == GRB.OPTIMAL:
    print("\nOptimal solution found:")
    for cycle in gurobi_model.getVars():
        if cycle.x > 0.5:
            print(f"Cycle: {cycle.varName}")
else:
    print("No optimal solution found.")