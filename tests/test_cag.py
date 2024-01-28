import pytest
from ec_ecology_toolbox import community_assembly_graph, example_nodes


def create_chr_nodes(graph):
    new_graph = {}
    for key1, val1 in graph.items():
        temp_dict = {}
        for key2, val2 in val1.items():
            temp_dict[example_nodes.Chr_Node(key2)] = val2

        new_graph[example_nodes.Chr_Node(key1)] = temp_dict

    return new_graph


def create_set_nodes(graph):
    new_graph = {}
    for key1, val1 in graph.items():
        temp_dict = {}
        for key2, val2 in val1.items():
            temp_dict[example_nodes.Set_Node(key2)] = val2

        new_graph[example_nodes.Set_Node(key1)] = temp_dict

    return new_graph


@pytest.mark.exclude
def test_graph_output():
    def custom_user_func(n):
        # Example user function used as input
        graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
        edges = graph[n.val]
        new_edges = {}
        for key, val in edges.items():
            new_edges[example_nodes.Chr_Node(key)] = val
        return new_edges

    expected1 = create_chr_nodes({1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}})
    expected2 = create_chr_nodes({1: {2: 0.4, 3: 0.6}, 3: {5: 0.7}, 5: {}})
    expected3 = create_chr_nodes({1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 5: {}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 5, False)
    result2 = community_assembly_graph.CAG(start_n, custom_user_func, 3, False)
    result3 = community_assembly_graph.CAG(start_n, custom_user_func, 4, False)

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3


def test_graph_output_cycle():
    def custom_user_func(n):
        # Example user function used as input
        graph = {1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
                 5: {4: 0.8, 8: 0.2}, 6: {}, 7: {}, 8: {}}
        edges = graph[n.val]
        new_edges = {}
        for node, edge in edges.items():
            new_edges[example_nodes.Chr_Node(node)] = edge
        return new_edges

    expected1 = create_chr_nodes({1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
                                  5: {4: 0.8, 8: 0.2}, 6: {}, 7: {}, 8: {}})
    expected2 = create_chr_nodes({1: {2: 0.3, 3: 0.7}, 3: {5: 0.7, 6: 0.2}, 5: {4: 0.8, 8: 0.2}})
    expected3 = create_chr_nodes({1: {2: 0.3, 3: 0.7}, 3: {5: 0.7, 6: 0.2}, 5: {4: 0.8, 8: 0.2}, 4: {2: 0.5, 7: 0.2}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 8, False)
    result2 = community_assembly_graph.CAG(start_n, custom_user_func, 3, False)
    result3 = community_assembly_graph.CAG(start_n, custom_user_func, 4, False)

    temp_graph = create_chr_nodes({1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
                                   5: {4: 0.8, 8: 0.2}, 6: {6: 1.0}, 7: {7: 1.0}, 8: {8: 1.0}})
    prob_list, node_dir = community_assembly_graph.random_walk_matrix_multiplication(temp_graph, start_n, 1500)
    result4 = {}
    for i in node_dir:
        result4[i] = round(prob_list[node_dir[i]], 8)

    exp_4 = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0.16101695, 7: 0.09016949, 8: 0.11271186}
    expected4 = {}
    for key, val in exp_4.items():
        expected4[example_nodes.Chr_Node(key)] = val

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3
    assert result4 == expected4


def test_graph_output_normalized():
    def custom_user_func(n):
        graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
        edges = graph[n.val]
        new_edges = {}
        for node, edge in edges.items():
            new_edges[example_nodes.Chr_Node(node)] = edge
        return new_edges

    expected1 = create_chr_nodes(
        {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7, 3: 0.3}, 4: {4: 1.0}, 5: {5: 1.0}})
    expected2 = create_chr_nodes({1: {2: 0.4, 3: 0.6}, 3: {5: 0.7, 3: 0.3}, 5: {5: 1.0}})
    expected3 = create_chr_nodes({1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7, 3: 0.3}, 5: {5: 1.0}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 5, True)
    result2 = community_assembly_graph.CAG(start_n, custom_user_func, 3, True)
    result3 = community_assembly_graph.CAG(start_n, custom_user_func, 4, True)

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3


def test_graph_output_empty_queue():
    def custom_user_func(n):
        graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
        edges = graph[n.val]
        new_edges = {}
        for node, edge in edges.items():
            new_edges[example_nodes.Chr_Node(node)] = edge
        return new_edges

    expected1 = create_chr_nodes({1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 50, False)

    assert result1 == expected1


def test_graph_output_inactive_nodes_in_queue():
    def custom_user_func(n):
        graph = {1: {2: 0.25, 3: 0.75}, 2: {}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.1}, 5: {}}
        edges = graph[n.val]
        new_edges = {}
        for node, edge in edges.items():
            new_edges[example_nodes.Chr_Node(node)] = edge
        return new_edges

    expected1 = create_chr_nodes({1: {2: 0.25, 3: 0.75}, 2: {}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.1}, 5: {}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 10, False)

    assert result1 == expected1


def test_graph_output_full_graph_explored():
    def custom_user_func(n):
        graph = {1: {2: 0.1, 3: 0.9}, 2: {5: 0.8}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.7}, 5: {}}
        edges = graph[n.val]
        new_edges = {}
        for node, edge in edges.items():
            new_edges[example_nodes.Chr_Node(node)] = edge
        return new_edges

    expected1 = create_chr_nodes({1: {2: 0.1, 3: 0.9}, 2: {5: 0.8}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.7}, 5: {}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 10, False)

    assert result1 == expected1


def test_graph_output_general():
    def custom_user_func(n):
        # Example user function used as input
        graph = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6},
                 4: {7: 0.4}, 5: {4: 0.5, 10: 0.2}, 6: {5: 0.6}, 7: {8: 0.4}, 8: {2: 0.5},
                 9: {6: 0.5}, 10: {9: 0.7}}
        edges = graph[n.val]
        new_edges = {}
        for node, edge in edges.items():
            new_edges[example_nodes.Chr_Node(node)] = edge
        return new_edges

    expected1 = create_chr_nodes({1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6},
                                  4: {7: 0.4}, 5: {4: 0.5, 10: 0.2}, 6: {5: 0.6}, 7: {8: 0.4}, 8: {2: 0.5},
                                  9: {6: 0.5}, 10: {9: 0.7}})
    expected2 = create_chr_nodes({1: {2: 0.3, 3: 0.6, 8: 0.1}, 3: {5: 0.2, 6: 0.6}, 5: {4: 0.5, 10: 0.2},
                                  6: {5: 0.6}})
    expected3 = create_chr_nodes({1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6, 3: 0.2},
                                  4: {7: 0.4, 4: 0.6}, 5: {4: 0.5, 10: 0.2, 5: 0.3}, 6: {5: 0.6, 6: 0.4},
                                  7: {8: 0.4, 7: 0.6}})

    start_n = example_nodes.Chr_Node(1)

    result1 = community_assembly_graph.CAG(start_n, custom_user_func, 10, False)
    result2 = community_assembly_graph.CAG(start_n, custom_user_func, 4, False)
    result3 = community_assembly_graph.CAG(start_n, custom_user_func, 7, True)

    temp_graph = create_chr_nodes({1: {2: 0.3, 3: 0.6, 8: 0.1}, 3: {5: 0.2, 6: 0.6}, 6: {5: 0.6}, 5: {4: 0.5, 10: 0.2},
                                   2: {3: 0.2, 4: 0.4, 7: 0.4}, 8: {8: 1.0}, 10: {10: 1.0}, 4: {4: 1.0}, 7: {7: 1.0}})
    prob_list, node_dir = community_assembly_graph.random_walk_matrix_multiplication(temp_graph, start_n, 1500)
    result4 = {}
    for i in node_dir:
        result4[i] = round(prob_list[node_dir[i]], 8)

    exp_4 = {1: 0, 2: 0, 3: 0, 4: 0.3048, 5: 0, 6: 0, 7: 0.12, 8: 0.1, 10: 0.07392}
    expected4 = {}
    for key, val in exp_4.items():
        expected4[example_nodes.Chr_Node(key)] = val

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3
    assert result4 == expected4


def test_graph_output_set_node():
    import copy
    n_species = 10

    def coexistence(n):
        edges = {}
        for i in range(n_species):
            if i not in n.members:
                new_node = copy.deepcopy(n)
                new_node.members.add(i)
                edges[new_node] = 1 / (n_species - len(n.members))
        return edges

    expected = create_set_nodes({(3, 4, 7): {(0, 3, 4, 7): 0.14285714285714285, (1, 3, 4, 7): 0.14285714285714285,
                                             (2, 3, 4, 7): 0.14285714285714285, (3, 4, 5, 7): 0.14285714285714285,
                                             (3, 4, 6, 7): 0.14285714285714285, (8, 3, 4, 7): 0.14285714285714285,
                                             (9, 3, 4, 7): 0.14285714285714285},
                                 (0, 3, 4, 7): {(0, 1, 3, 4, 7): 0.16666666666666666,
                                                (0, 2, 3, 4, 7): 0.16666666666666666,
                                                (0, 3, 4, 5, 7): 0.16666666666666666,
                                                (0, 3, 4, 6, 7): 0.16666666666666666,
                                                (0, 3, 4, 7, 8): 0.16666666666666666,
                                                (0, 3, 4, 7, 9): 0.16666666666666666},
                                 (2, 3, 4, 7): {(0, 2, 3, 4, 7): 0.16666666666666666,
                                                (1, 2, 3, 4, 7): 0.16666666666666666,
                                                (2, 3, 4, 5, 7): 0.16666666666666666,
                                                (2, 3, 4, 6, 7): 0.16666666666666666,
                                                (2, 3, 4, 7, 8): 0.16666666666666666,
                                                (2, 3, 4, 7, 9): 0.16666666666666666}})

    start_n = example_nodes.Set_Node({3, 4, 7})

    result = community_assembly_graph.CAG(start_n, coexistence, 3, False)

    assert result == expected
