import pytest
from ec_ecology_toolbox import community_assembly_graph, example_nodes


@pytest.mark.exclude
def test_chr_nodes():
    node1 = example_nodes.Chr_Node(1)
    node2 = example_nodes.Chr_Node(2)
    node1_dummy = example_nodes.Chr_Node(1)

    assert str(node1) == str(1)
    assert str(node2) == str(2)

    assert node1 == node1_dummy
    assert node1 != node2

    assert node1 < node2
    assert node2 > node1

    assert node1 <= node2
    assert node1 <= node1_dummy

    assert node2 >= node1

    temp_dict = {node1: node1.val, node2: node2.val}
    assert temp_dict[node1] == node1.val
    assert temp_dict[node2] == node2.val


def test_set_nodes():
    node1 = example_nodes.Set_Node({7, 9})
    node2 = example_nodes.Set_Node({1, 4, 6})
    node1_dummy = example_nodes.Set_Node({7, 9})

    assert str(node1) == str({7, 9})
    assert str(node2) == str({1, 4, 6})

    assert node1 == node1_dummy
    assert node1 != node2

    assert node1 < node2
    assert node2 > node1

    assert node1 <= node2
    assert node1 <= node1_dummy

    assert node2 >= node1

    temp_dict = {node1: node1.members, node2: node2.members}
    assert temp_dict[node1] == node1.members
    assert temp_dict[node2] == node2.members
