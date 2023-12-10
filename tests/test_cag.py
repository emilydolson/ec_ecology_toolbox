import pytest
from ec_ecology_toolbox import community_assembly_graph
from queue import PriorityQueue


@pytest.mark.exclude
def test_graph_output():
    def custom_user_func(n):
        # Example user function used as input
        graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
        return graph[n]

    expected1 = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
    expected2 = {1: {2: 0.4, 3: 0.6}, 3: {5: 0.7}, 5: {}}
    expected3 = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 5: {}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 5, False)
    result2 = community_assembly_graph.CAG(1, custom_user_func, 3, False)
    result3 = community_assembly_graph.CAG(1, custom_user_func, 4, False)

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3


def test_graph_output_cycle():
    def custom_user_func(n):
        # Example user function used as input
        graph = {1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
                 5: {4: 0.8, 8: 0.2}, 6: {}, 7: {}, 8: {}}
        return graph[n]

    expected1 = {1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
                 5: {4: 0.8, 8: 0.2}, 6: {}, 7: {}, 8: {}}
    expected2 = {1: {2: 0.3, 3: 0.7}, 3: {5: 0.7, 6: 0.2}, 5: {4: 0.8, 8: 0.2}}
    expected3 = {1: {2: 0.3, 3: 0.7}, 3: {5: 0.7, 6: 0.2}, 5: {4: 0.8, 8: 0.2}, 4: {2: 0.5, 7: 0.2}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 8, False)
    result2 = community_assembly_graph.CAG(1, custom_user_func, 3, False)
    result3 = community_assembly_graph.CAG(1, custom_user_func, 4, False)

    temp_graph = {1: {2: 0.3, 3: 0.7}, 2: {3: 0.2}, 3: {5: 0.7, 6: 0.2}, 4: {2: 0.5, 7: 0.2},
                 5: {4: 0.8, 8: 0.2}, 6: {6: 1.0}, 7: {7: 1.0}, 8: {8: 1.0}}
    prob_list, node_dir = community_assembly_graph.random_walk_matrix_multiplication(temp_graph, 1, 1500)
    result4 = {}
    for i in node_dir:
        result4[i] = round(prob_list[node_dir[i]], 8)

    expected4 = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0.16101695, 7: 0.09016949, 8: 0.11271186}

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3
    assert result4 == expected4


def test_graph_output_normalized():
    def custom_user_func(n):
        graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
        return graph[n]

    expected1 = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7, 3: 0.3}, 4: {4: 1.0}, 5: {5: 1.0}}
    expected2 = {1: {2: 0.4, 3: 0.6}, 3: {5: 0.7, 3: 0.3}, 5: {5: 1.0}}
    expected3 = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7, 3: 0.3}, 5: {5: 1.0}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 5, True)
    result2 = community_assembly_graph.CAG(1, custom_user_func, 3, True)
    result3 = community_assembly_graph.CAG(1, custom_user_func, 4, True)

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3


def test_graph_output_empty_queue():
    def custom_user_func(n):
        graph = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}
        return graph[n]

    expected1 = {1: {2: 0.4, 3: 0.6}, 2: {3: 0.2, 4: 0.5, 5: 0.3}, 3: {5: 0.7}, 4: {}, 5: {}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 50, False)

    assert result1 == expected1


def test_graph_output_inactive_nodes_in_queue():
    def custom_user_func(n):
        graph = {1: {2: 0.25, 3: 0.75}, 2: {}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.1}, 5: {}}
        return graph[n]

    expected1 = {1: {2: 0.25, 3: 0.75}, 2: {}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.1}, 5: {}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 10, False)

    assert result1 == expected1


def test_graph_output_full_graph_explored():
    def custom_user_func(n):
        graph = {1: {2: 0.1, 3: 0.9}, 2: {5: 0.8}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.7}, 5: {}}
        return graph[n]

    expected1 = {1: {2: 0.1, 3: 0.9}, 2: {5: 0.8}, 3: {2: 0.5, 4: 0.4}, 4: {5: 0.7}, 5: {}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 10, False)

    assert result1 == expected1


def test_graph_output_general():
    def custom_user_func(n):
        # Example user function used as input
        graph = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6},
                 4: {7: 0.4}, 5: {4: 0.5, 10: 0.2}, 6: {5: 0.6}, 7: {8: 0.4}, 8: {2: 0.5},
                 9: {6: 0.5}, 10: {9: 0.7}}
        return graph[n]

    expected1 = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6},
                 4: {7: 0.4}, 5: {4: 0.5, 10: 0.2}, 6: {5: 0.6}, 7: {8: 0.4}, 8: {2: 0.5},
                 9: {6: 0.5}, 10: {9: 0.7}}
    expected2 = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 3: {5: 0.2, 6: 0.6}, 5: {4: 0.5, 10: 0.2},
                 6: {5: 0.6}}
    expected3 = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 2: {3: 0.2, 4: 0.4, 7: 0.4}, 3: {5: 0.2, 6: 0.6, 3: 0.2},
                 4: {7: 0.4, 4: 0.6}, 5: {4: 0.5, 10: 0.2, 5: 0.3}, 6: {5: 0.6, 6: 0.4}, 7: {8: 0.4, 7: 0.6}}

    result1 = community_assembly_graph.CAG(1, custom_user_func, 10, False)
    result2 = community_assembly_graph.CAG(1, custom_user_func, 4, False)
    result3 = community_assembly_graph.CAG(1, custom_user_func, 7, True)

    temp_graph = {1: {2: 0.3, 3: 0.6, 8: 0.1}, 3: {5: 0.2, 6: 0.6}, 6: {5: 0.6}, 5: {4: 0.5, 10: 0.2},
                  2: {3: 0.2, 4: 0.4, 7: 0.4}, 8: {8: 1.0}, 10: {10: 1.0}, 4: {4: 1.0}, 7: {7: 1.0}}
    prob_list, node_dir = community_assembly_graph.random_walk_matrix_multiplication(temp_graph, 1, 1500)
    result4 = {}
    for i in node_dir:
        result4[i] = round(prob_list[node_dir[i]], 8)

    expected4 = {1: 0, 2: 0, 3: 0, 4: 0.3048, 5: 0, 6: 0, 7: 0.12, 8: 0.1, 10: 0.07392}

    assert result1 == expected1
    assert result2 == expected2
    assert result3 == expected3
    assert result4 == expected4

