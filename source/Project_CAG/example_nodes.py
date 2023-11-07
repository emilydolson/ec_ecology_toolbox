# Example classes for each node within the CAG

class Node:
    # Represented by a set of members
    def __init__(self):
        self.members = {}

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
