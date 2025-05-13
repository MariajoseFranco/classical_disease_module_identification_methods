import os
import sys

# Absolute path to the folder containing the module
sys.path.append(os.path.abspath('./state_of_art_repos'))

from DIAMOnD.DIAMOnD import run_diamond_from_args


class DIAMOND:
    def __init__(self, alpha=1):
        self.alpha = alpha

    def run_diamond(self, ppi, seed_nodes, n):
        seed_genes, added_nodes = run_diamond_from_args([
            ppi,
            seed_nodes,
            str(n)
        ])
        added_genes = set([gene[0] for gene in added_nodes])
        result = {
            'seed_nodes': list(seed_genes)
        }
        result['seed_nodes_module_1'] = list(added_genes)
        return result
