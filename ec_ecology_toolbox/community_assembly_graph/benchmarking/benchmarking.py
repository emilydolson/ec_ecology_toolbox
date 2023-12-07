import random
from ec_ecology_toolbox.community_assembly_graph.example_nodes import Chr_Node


def benchmarking_func1(n):
    if n.val > 50000000:
        return {}
    elif n.val == 50000000:
        return {Chr_Node(100000000): 0.5}
    return {Chr_Node(n*2): 0.5, Chr_Node(n*2+1): 0.5}


def benchmarking_func2(n):
    random_int = random.randint(0, 15)
    edges = {}
    for i in range(random_int):
        node = random.randint(1, 2500)
        weight = round(random.random(), 2)
        edges[Chr_Node(node)] = weight
    return edges
