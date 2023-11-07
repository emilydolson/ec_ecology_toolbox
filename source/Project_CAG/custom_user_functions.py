import copy


def coexistence(n):
    # Example user function used as input
    n_species = 10
    edges = {}
    for i in range(n_species):
        if i not in n.members:
            new_node = copy.deepcopy(n)
            new_node.members.add(i)
            edges[new_node] = 1/(n_species - len(n.members))
    return edges


def custom_user_func_1(n):
    # Example user function used as input
    graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
    return graph[n]


def custom_user_func_2(n):
    # Example user function used as input
    graph = {1: {2: 0.3, 3: 0.7}, 2: {2: 0.8, 3: 0.2}, 3: {5: 0.7, 6: 0.3}, 4: {2: 0.5, 4: 0.5},
             5: {4: 0.8, 5: 0.2}, 6: {5: 0.5, 6: 0.5}}
    return graph[n]


def custom_user_func_3(n):
    # Example user function used as input
    graph = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6},
             4: {7: 0.4}, 5: {4: 0.5, 10: 0.2}, 6: {5: 0.6}, 7: {8: 0.4}, 8: {2: 0.5},
             9: {6: 0.5}, 10: {9: 0.7}}
    return graph[n]


def custom_user_func_4(n):
    # Example user function used as input
    graph = {1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
             5: {4: 0.8, 8: 0.2}, 6: {}, 7: {}, 8: {}}
    return graph[n]
