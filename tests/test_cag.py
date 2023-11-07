import pytest
from ec_ecology_toolbox.source.Project_CAG import community_assembly_graph


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
