import copy

class Node:
    def __init__(self):
        self.members = set()

def coexistance(n):
    n_species = 10
    edges = {}
    for i in range(n_species):
        if i not in n.members:
            new_node = copy.deepcopy(n)
            new_node.members.add(i)
            edges[new_node] = 1/(n_species - len(n.members))
    return edges
