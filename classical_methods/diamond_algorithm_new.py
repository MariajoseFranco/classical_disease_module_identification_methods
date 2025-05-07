import os
import sys

# Absolute path to the folder containing the module
sys.path.append(os.path.abspath('/Users/mariajosefranco/Desktop/Data Science - UPM/TFM/project/state_of_art_repos'))

from DIAMOnD.DIAMOnD import run_diamond_from_args


class DIAMOND:
    def __init__(self, alpha=1):
        self.alpha = alpha

    def run_diamond(self, ppi, seed_nodes, n):
        added_nodes = run_diamond_from_args([
            ppi,
            seed_nodes,
            str(n)
        ])
        return added_nodes
