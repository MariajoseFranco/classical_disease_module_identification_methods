import os
import sys
import tempfile
from pathlib import Path

import networkx as nx
import pandas as pd

# Absolute path to the folder containing the module
sys.path.append(os.path.abspath('/Users/mariajosefranco/Desktop/Data Science - UPM/TFM/project/state_of_art_repos/DOMINO'))
from src.core.domino import main as domino_main


class DOMINO:
    def __init__(self):
        pass

    def _run_domino_pipeline(self, G, seed_nodes):
        # Save PPI to temp file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as ppi_file:
            for u, v in G.edges():
                ppi_file.write(f"{u}\t{v}\n")
            ppi_path = ppi_file.name

        # Save seeds to temp file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as seeds_file:
            for node in seed_nodes:
                seeds_file.write(f"{node}\n")
            seeds_path = seeds_file.name

        # Call domino_main with paths
        domino_main(seeds_path, ppi_path)

        # Clean up
        os.remove(ppi_path)
        os.remove(seeds_path)

    def run_domino(self, ppi_path: str, seed_nodes_path: str):
        # Step 1: Load inputs
        if not Path(ppi_path).is_file():
            raise FileNotFoundError(f"PPI file not found: {ppi_path}")
        if not Path(seed_nodes_path).is_file():
            raise FileNotFoundError(f"Seed node file not found: {seed_nodes_path}")

        # Load PPI as dataframe
        ppi_df = pd.read_csv(ppi_path, sep="\t", header=None)
        ppi_df.columns = ["protein1", "protein2"]

        # Load seed nodes
        with open(seed_nodes_path, "r") as f:
            seed_nodes = [line.strip() for line in f if line.strip()]

        # Build graph
        G = nx.from_pandas_edgelist(ppi_df, source="protein1", target="protein2")

        # Run the DOMINO pipeline
        self._run_domino_pipeline(G, seed_nodes)
