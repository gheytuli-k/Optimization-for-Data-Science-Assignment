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

