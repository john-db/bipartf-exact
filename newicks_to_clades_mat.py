from utils import newick_to_subtrees
import numpy as np

path = "/data/CDSLSahinalp/bridgersjd/bpf-exact/10_leaf_clts.newick"
with open(path) as f:
    for tree in f:
        print(''.join(map(str, np.array(newick_to_subtrees(tree)).astype(int).flatten())))
