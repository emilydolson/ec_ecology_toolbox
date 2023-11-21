import sys
import csv
from collections import Counter
import json
import ec_ecology_toolbox as eco
import numpy as np
import itertools
import heapq
import glob

gene_phen_map = {}
phen_str_to_phen = {}


cut_off = int(sys.argv[2])
file_root = sys.argv[1].replace("/", "_").replace("*","") + "_cutoff_" + str(cut_off)

def print_phenotype_adjacency_list(map):
    outfile = open(file_root + "phen_adjaceny.csv", "w")
    node_file = open(file_root + "phen_info.csv", "w")
    outfile.write("from, to, weight\n")
    node_file.write("id, phenotype, neutrality\n")
    for el in map:
        denom = sum(map[el].full_adj.values())
        if denom != 0:
            neut = map[el].full_adj[el]/denom
        else:
            neut = map[el].full_adj[el]
        id = map[el].id
        node_file.write(",".join([str(id), str(el), str(neut)])+"\n")

        for a in map[el].full_adj:
            outfile.write(",".join([str(i) for i in [map[el].id, a.id, map[el].full_adj[a]/denom]]) + "\n")
    outfile.close()
    node_file.close()

def print_community_adjacency_list(map):
    edgefile = open(file_root + "comm_edges.csv", "w")
    nodefile = open(file_root + "comm_nodes.csv", "w")
    edgefile.write("From, To, Weight\n")
    nodefile.write("id, members, n_members, evaluated, actual_sink\n")
    for el in map:
        if not map[el].stable:
            continue
        members = str([i.id for i in map[el].members]).replace(",", " ")
        mem1 = str(map[el].id)
        actual_sink = map[el].evaluated and sum(map[el].edges.values()) == 0
        nodefile.write(",".join([mem1, members, str(len(map[el].members)), str(int(map[el].evaluated)), str(actual_sink)]) + "\n")
        denom = sum(map[el].edges.values())
        # print()
        # print(mem1, map[el].edges)
        for a in map[el].edges:
            edge_count = map[el].edges[a]
            while a.transitions_to != a:
                a = a.transitions_to
            assert a.stable
            assert a.members in map
            # mem2 = str([i.id for i in a.members]).replace(",", " ")             
            mem2 = str(a.id)
            # print(mem2)
            edgefile.write(",".join([str(i) for i in [mem1, mem2, edge_count/denom]]) + "\n")
        # if mem1 == "0":
        #     exit()
    edgefile.close()
    nodefile.close()

class Phenotype:
    next_id = 0

    def __init__(self, phen):
        self.phen_str = phen
        self.phen = json.loads(phen)
        self.id = Phenotype.next_id
        self.dominates = set()
        Phenotype.next_id += 1
        # print(self.id, Phenotype.next_id)
        self.adj = Counter()
        self.full_adj = Counter()

    def add_adj(self, other):
        # print(self.phen, other.phen)
        self.full_adj[other] += 1
        if not all([self.phen[i] >= other.phen[i] for i in range(len(self.phen))]):
            # print("Adding!")
            self.adj[other] += 1
        elif other.phen != self.phen:
            self.dominates.add(other)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self) -> int:
        return hash(self.id)

    def __repr__(self) -> str:
        return str(self.id) + ": " + self.phen_str


class Community:
    next_id = 0
    S = 1000
    G = 5
    thresh = .5
    seen = {}
    queue = []

    def __init__(self, members, priority, source, transitions_to=None):
        self.id = Community.next_id
        Community.next_id += 1
        self.added_to_queue = False
        self.members = frozenset(members)
        self.edges = {}
        self.in_priorities = {source: priority}
        self.priority = priority
        self.evaluated = False
        assert self.members not in Community.seen
        # if self.members in Community.seen:
        #     self.transitions_to = Community.seen[self.members].members

        if transitions_to is None:

            dominated = set(itertools.chain(*[i.dominates for i in members]))
            transitions_to = list(self.members - dominated)
            # print("transitions_to", transitions_to, dominated, self.members)
            if frozenset(transitions_to) in Community.seen:
                self.transitions_to = Community.seen[frozenset(transitions_to)]
            else:
                count = 0
                transient_states = [frozenset(transitions_to)]
                P_survival = [0]

                while any([i < Community.thresh for i in P_survival]) and count < Community.G:
                    # print("looping")
                    count += 1
                    lex_input = [m.phen for m in transitions_to]
                    plex = eco.LexicaseFitness(lex_input)
                    P_survival = (np.ones(len(plex)) - (np.ones(len(plex)) - plex)**Community.S)**Community.G
                    survivors = frozenset([transitions_to[i] for i in range(len(transitions_to)) if P_survival[i] > Community.thresh])
                    # print("survivors", survivors)
                    if survivors in Community.seen:
                        self.transitions_to = Community.seen[survivors].transitions_to
                        break

                    transitions_to = list(survivors)
                    transient_states.append(survivors)

                if frozenset(transitions_to) not in Community.seen:
                    if frozenset(transitions_to) == self.members:
                        Community.seen[frozenset(transitions_to)] = self
                    else:
                        Community.seen[frozenset(transitions_to)] = Community(transitions_to, priority, self, transitions_to)
                    self.transitions_to = Community.seen[frozenset(transitions_to)]
                for state in transient_states:
                    if state not in Community.seen and state != self.members:
                        Community.seen[state] = Community(state, priority, self, transitions_to)
        else:
            if frozenset(transitions_to) == self.members:
                self.transitions_to = self
            else:
                self.transitions_to = Community.seen[frozenset(transitions_to)]
        self.stable = (self.members == self.transitions_to.members)
        if (not self.stable):
            self.edges[self.transitions_to] = 1
        self.tested = True
        Community.seen[self.members] = self
        # print("Add?", to_add, to_add.added_to_queue)
        if self.transitions_to.stable and not self.transitions_to.added_to_queue:
            # print("added", self.transitions_to, self.transitions_to.priority)
            heapq.heappush(Community.queue, (-1 * priority * self.transitions_to.priority, self.transitions_to))
            self.transitions_to.added_to_queue = True

    def __eq__(self, __value: object) -> bool:
        return self.members == __value.members
    
    def __hash__(self) -> int:
        return hash(self.members)
    
    def __repr__(self) -> str:
        return "Stable: " + str(self.stable) + ", Members: " + str([i.id for i in self.members])

    def __lt__(self, other):
        return len(self.members) < len(other.members)

    def update_prob(self, prob, source, already_done=None):
        if already_done is None:
            already_done = set()
        already_done.add(self)
        self.in_priorities[source] = prob
        self.priority = sum(self.in_priorities.values())
        if not self.evaluated:
            heapq.heappush(Community.queue, (-1 * self.priority, self))            
        denom = sum(self.edges.values())
        for adj in self.edges:
            if adj not in already_done:
                adj.update_prob(prob * self.edges[adj]/denom, self, already_done)

for filename in glob.glob(sys.argv[1]):
    with open(filename) as infile:
        count = 0
        reader = csv.reader(infile)
        header = next(reader)
        genotype_id_col = header.index("genotype_id")
        original_col = header.index("is_original")
        scores_col = header.index("training_case_scores")

        for line in reader:
            # if count > 100:
            #     break
            # count += 1
            phen_str = line[scores_col]
            if phen_str not in phen_str_to_phen:
                phen = Phenotype(phen_str)
                phen_str_to_phen[phen_str] = phen
                # id_to_phen[next_id] = phen
            else:
                phen = phen_str_to_phen[phen_str]

            genotype_id = line[genotype_id_col]
            if line[original_col] == "1":
                gene_phen_map[genotype_id] = phen
            else:
                # phen.add_adj(gene_phen_map[genotype_id])
                gene_phen_map[genotype_id].add_adj(phen)

print_phenotype_adjacency_list(phen_str_to_phen)

#curr_str = "[1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,1,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,0,0,1,0,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,1,0,1,1]"
curr_str = str([0]*100).replace(" ", "")
print(phen_str_to_phen[curr_str])
curr_phen = phen_str_to_phen[curr_str]
comm = Community([curr_phen], 1, None)

count = 0
# print(count, Community.queue)

while Community.queue and count < cut_off:
    count += 1
    # print(Community.queue)

    priority, curr = heapq.heappop(Community.queue)
    # print(curr)
    while curr.evaluated:
        if Community.queue:
            curr = heapq.heappop(Community.queue)[1]
        else:
            print_community_adjacency_list(Community.seen)
            exit(0)

    curr.evaluated = True
    # print(len(queue), curr)
    adj = Counter()
    # print(curr.members)
    dominated = set(itertools.chain(*[i.dominates for i in curr.members]))
    for phen in curr.members:
        # print("Finding adj", phen, phen.adj.keys())
        adj.update({i: phen.adj[i] + adj[i] for i in phen.adj.keys() if i not in dominated})

    print(len(Community.queue), priority, len(adj), curr)
    # print(adj)
    for phen in adj:
        # print("adding", phen)
        members = set(list(curr.members) + [phen])
        new_members = frozenset(members - set(phen.dominates))

        if (new_members not in Community.seen):
            # print("not seen")

            c = Community(new_members, curr.priority * adj[phen]/sum(adj.values()), curr)
        else:
            # print("seen")
            c = Community.seen[new_members]
            c.update_prob(curr.priority * adj[phen]/sum(adj.values()), curr)

        if c in curr.edges:
            curr.edges[c] += adj[phen]/sum(adj.values())
        else:
            # print("Adding", c, "to", curr)
            # exit()
            curr.edges[c] = adj[phen]/sum(adj.values())

print_community_adjacency_list(Community.seen)