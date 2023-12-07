import pytest
from ec_ecology_toolbox import community_assembly_graph


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

