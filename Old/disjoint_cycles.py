def find_disjoint_cycles(graph: dict[str, dict[str, int]], n: int) -> list[list[str]]:
    """
    Find all the disjoint cycles of length n in the graph
    :param graph: adjacency list representation of graph
        - key: node
        - value: Dictionary of neighbors and the weight of the edge
            - key: neighbor
            - value: weight
    :param n: length of cycles
    :return: list of disjoint cycles
    """
    # Set to track nodes used in disjoint cycles
    used_nodes = set()
    disjoint_cycles = []

    def dfs(node, start, depth, visited, path):
        if depth == n:
            if start in graph[node]:
                # Check if the cycles is disjoint
                if all(p not in used_nodes for p in path):
                    disjoint_cycles.append(path + [start])
                    used_nodes.update(path)  # Mark these nodes as used

            return

        visited.add(node)

        for neighbor in graph[node]:
            if neighbor not in visited and neighbor not in used_nodes:
                dfs(neighbor, start, depth + 1, visited, path + [neighbor])

        visited.remove(node)

    for node in graph:
        dfs(node, node, 1, set(), [node])

    return disjoint_cycles
