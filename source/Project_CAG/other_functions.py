import networkx as nx


# Unused functions created while developing the code

def detect_cycle(start_n, user_func):
    # Function used to detect cycles in a graph
    # Uses rabbit-hare algorithm
    def get_next_node(node):
        next_nodes = user_func(node)
        return next(iter(next_nodes.keys())) if next_nodes else None

    slow_pointer = fast_pointer = start_n

    while fast_pointer is not None and get_next_node(fast_pointer) is not None:
        slow_pointer = get_next_node(slow_pointer)  # Move one step
        fast_pointer = get_next_node(get_next_node(fast_pointer))  # Move two steps

        if fast_pointer is not None and slow_pointer == fast_pointer:
            cycle_start = start_n
            while cycle_start != slow_pointer:
                cycle_start = get_next_node(cycle_start)
                slow_pointer = get_next_node(slow_pointer)

            cycle_nodes = [slow_pointer]
            slow_pointer = get_next_node(slow_pointer)

            while slow_pointer != cycle_nodes[0]:
                cycle_nodes.append(slow_pointer)
                slow_pointer = get_next_node(slow_pointer)

            return cycle_nodes

    return None  # No cycle detected


def detect_cycle_2(graph):
    # Function used to detect cycles in the input graph
    # Uses networkx library
    G = nx.DiGraph(graph)
    return list(nx.simple_cycles(G))


def pagerank_algo(graph):
    # Uses a pagerank algorithm to output probabilities of ending up at a node in a given graph
    G = nx.DiGraph()
    for node, neighbors in graph.items():
        for neighbor, weight in neighbors.items():
            G.add_edge(node, neighbor, weight=weight)

    pr = nx.pagerank(G, weight='weight')

    return pr


def calc_cycle_r(cycle, edge_list):
    # Function used to calculate changes in probability caused as a result of a cycle
    # Based on the principle that a cycle introduces an infinite geometric progression series
    cycle.append(cycle[0])
    r = 1
    for idx, cycle_node in enumerate(cycle):
        if idx == (len(cycle) - 1):
            break

        r *= edge_list[cycle_node][cycle[idx + 1]]

    return r


# code used to update and handle cycles
#
# detected_cycles = set()
#
# cycles = detect_cycle_2(edge_list)
# cycles = set([frozenset(c) for c in cycles if not frozenset(c) in detected_cycles])
# if cycles:
#     flag = False
#     for c in cycles:
#         if len(c) > 1:
#             flag = True
#             break
#     detected_cycles = detected_cycles | cycles
#
#     if flag:
#         p_queue = update_priorities(edge_list, node_passing_p, p_queue, visited)

