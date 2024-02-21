# Example classes for each node within the CAG

class Set_Node:
    """
    Sample class for user defined Node in graphs
    Nodes are characterized based on a set of elements
    Equality of nodes are based on having the same set of elements
    Inequality is compared using length of the set i.e. number of elements in the set
    Node class is hashable and uses frozen sets
    """
    def __init__(self, *args):
        self.members = set()
        if len(args) == 1:
            self.members = set(args[0])

    def __repr__(self):
        return str(self.members)

    def __eq__(self, other):
        return self.members == other.members

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return len(self.members) < len(other.members)

    def __gt__(self, other):
        return len(self.members) > len(other.members)

    def __le__(self, other):
        return len(self.members) < len(other.members) or self.members == other.members

    def __ge__(self, other):
        return len(self.members) > len(other.members) or self.members == other.members

    def __hash__(self):
        return hash(frozenset(self.members))


class Chr_Node:
    """
    Sample class for user defined Node in graphs
    Nodes are characterized based on a single character e.g. integers
    Equality of nodes are based on having the same character
    Inequality is compared using number values if integers, ASCII values if string character
    Node class is hashable
    """

    def __init__(self, *args):
        self.val = None
        if len(args) == 1:
            self.val = args[0]

    def __repr__(self):
        return str(self.val)

    def __eq__(self, other):
        return self.val == other.val

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.val < other.val

    def __gt__(self, other):
        return self.val > other.val

    def __le__(self, other):
        return self.val < other.val or self.val == other.val

    def __ge__(self, other):
        return self.val > other.val or self.val == other.val

    def __hash__(self):
        return hash(self.val)




