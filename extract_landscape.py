import sys
import csv
from collections import Counter
import json

gene_phen_map = {}
phen_str_to_id = {}
id_to_phen = {}


class Phenotype:
    def __init__(self, phen, id):
        self.phen_str = phen
        self.phen = json.loads(phen)
        self.id = id
        self.adj = Counter()

    def add_adj(self, other_id):
        self.adj[other_id] += 1

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self) -> int:
        return hash(self.id)

    def __repr__(self) -> str:
        return str(self.id) + ": " + self.phen_str


with open(sys.argv[1]) as infile:
    next_id = 0
    reader = csv.reader(infile)
    header = next(reader)
    genotype_id_col = header.index("genotype_id")
    original_col = header.index("is_original")
    scores_col = header.index("training_case_scores")

    for line in reader:
        phen_str = line[scores_col]
        if phen_str not in phen_str_to_id:
            phen_id = next_id
            phen_str_to_id[phen_str] = phen_id
            id_to_phen[next_id] = Phenotype(phen_str, next_id)
        else:
            phen_id = phen_str_to_id[phen_str]

        genotype_id = line[genotype_id_col]
        if line[original_col] == "1":
            gene_phen_map[genotype_id] = phen_id
        elif gene_phen_map[genotype_id] != phen_id:
            id_to_phen[phen_id].add_adj(gene_phen_map[genotype_id])
            id_to_phen[gene_phen_map[genotype_id]].add_adj(phen_id)

        next_id += 1

curr_str = "[1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,1,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,0,0,1,0,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,1,0,1,1]"
print(phen_str_to_id[curr_str])
index = phen_str_to_id[curr_str]
curr_phen = id_to_phen[index]
queue = [set([index])]

while queue and len(queue) < 10000000:
    print(len(queue))
    curr = queue.pop()
    curr_phens = [id_to_phen[i] for i in curr]
    adj = set()
    
    for phen in curr_phens:
        # print(phen, phen.adj)
        adj.update(phen.adj.keys())

    for phen_id in adj:
        # print(phen_id)
        queue.append(curr.union([phen_id]))