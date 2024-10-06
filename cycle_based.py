from aux_functions import Loader
from aux_functions import data_to_adjacency
import gurobipy as gp
from gurobipy import GRB
import networkx as nx

data = Loader("Instance Files\Delorme_50_NDD_Unit_2.txt")
# making adjacency out of the data and connecting all the pairs to god donor
adjacency = data_to_adjacency(data, True, 1)

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

G = nx.DiGraph(adjacency)

# Assign weight of 1 to all the edges
for u, v in G.edges():
    G[u][v]['weight'] = 1
cycles = find_cycles(G, 6)


# def optimize_cycles(cycles, vertex_weights):
#     # Create a new Gurobi model
#     model = gp.Model("Kidney_Exchange_Optimization")

#     # Step 1: Create variables for each cycle
#     # Binary decision variables x_c for each cycle in the cycles list
#     cycle_vars = {}
#     for i, cycle in enumerate(cycles):
#         cycle_vars[i] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}")
#     # Step 2: Set the objective to maximize the total weight of selected cycles
#     # Each cycle's weight is the sum of the weights of the arcs in the cycle
#     objective = gp.quicksum(vertex_weights[i] * cycle_vars[i] for i in range(len(cycles)))
#     model.setObjective(objective, GRB.MAXIMIZE)

#     # Step 3: Add constraints to ensure vertex-disjoint cycles
#     # For each vertex, ensure that it's in at most one cycle
#     vertex_cycles = {}
    
#     for i, cycle in enumerate(cycles):
#         for vertex in cycle:
#             if vertex not in vertex_cycles:
#                 vertex_cycles[vertex] = []
#             vertex_cycles[vertex].append(i)
    
#     for vertex, cycle_indices in vertex_cycles.items():
#         model.addConstr(gp.quicksum(cycle_vars[i] for i in cycle_indices) <= 1, name=f"vertex_{vertex}_constraint")

#     # Step 4: Optimize the model
#     model.optimize()

#     # Step 5: Retrieve and return the optimal solution
#     selected_cycles = [cycles[i] for i in range(len(cycles)) if cycle_vars[i].x > 0.5]
#     return selected_cycles

# vertex_weights = {i: len(cycle) for i, cycle in enumerate(cycles)}
# print(optimize_cycles(cycles, vertex_weights))

def optimize_cycles(cycles, G):
    # Create a new Gurobi model
    model = gp.Model("Kidney_Exchange_Optimization")

    # Step 1: Create variables for each cycle and compute their weights using edge weights
    cycle_vars = {}
    cycle_weights = {}
    for i, cycle in enumerate(cycles):
        # Compute the weight of the cycle as the sum of the weights of the arcs
        weight = 0
        for j in range(len(cycle) - 1):
            u = cycle[j]
            v = cycle[j + 1]
            if G.has_edge(u, v):
                weight += G[u][v]['weight']
            else:
                print(f"Edge from {u} to {v} not found in graph.")
        # If it's a cycle (returns to the starting node), add the weight of the closing edge
        if G.has_edge(cycle[-1], cycle[0]):
            weight += G[cycle[-1]][cycle[0]]['weight']
        # else:
        #     print(f"Edge from {cycle[-1]} to {cycle[0]} not found in graph.")
        cycle_weights[i] = weight
        cycle_vars[i] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}")

    # Step 2: Set the objective to maximize the total weight of selected cycles
    objective = gp.quicksum(cycle_weights[i] * cycle_vars[i] for i in range(len(cycles)))
    model.setObjective(objective, GRB.MAXIMIZE)

    # Step 3: Add constraints to ensure vertex-disjoint cycles
    vertex_cycles = {}
    for i, cycle in enumerate(cycles):
        for vertex in set(cycle):  # Use set to avoid double-counting the starting node
            if vertex not in vertex_cycles:
                vertex_cycles[vertex] = []
            vertex_cycles[vertex].append(i)

    for vertex, cycle_indices in vertex_cycles.items():
        model.addConstr(
            gp.quicksum(cycle_vars[i] for i in cycle_indices) <= 1,
            name=f"vertex_{vertex}_constraint"
        )

    # Step 4: Optimize the model
    model.optimize()

    # Step 5: Retrieve and return the optimal solution
    selected_cycles = [cycles[i] for i in range(len(cycles)) if cycle_vars[i].X > 0.5]
    total_weight = sum(cycle_weights[i] for i in range(len(cycles)) if cycle_vars[i].X > 0.5)
    print(f"Total weight of selected cycles: {total_weight}")
    return selected_cycles

print(optimize_cycles(cycles, G))