import os
import sys
import tempfile
from pathlib import Path

import networkx as nx
import pandas as pd

# Absolute path to the folder containing the module
sys.path.append(os.path.abspath('/Users/mariajosefranco/Desktop/Data Science - UPM/TFM/project/state_of_art_repos/DOMINO'))
from src.runner import main_domino, main_slicer

# from src.core.domino import main as main_domino
# from src.core.preprocess_slices import main as main_slicer


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

        # Create slices file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as slices_file:
            slices_path = slices_file.name

        # Create output folder
        output_folder = tempfile.mkdtemp()

        # Step 1: run slicer
        sys.argv = [
            "preprocess_slices.py",
            "-n", ppi_path,
            "-o", slices_path
        ]
        main_slicer()

        # Step 2: run DOMINO with full args
        sys.argv = [
            "runner.py",
            "-a", seeds_path,
            "-n", ppi_path,
            "-s", slices_path,
            "-o", output_folder,
            "-c", "false",              # avoid pickle cache
            "-p", "1",                  # single-thread
            "-v", "false",              # skip visualization
            "-sth", "0.3",              # slice threshold
            "-mth", "0.05"              # module threshold
        ]
        G_final_modules = main_domino()

        # Clean up if needed:
        os.remove(ppi_path)
        os.remove(seeds_path)
        os.remove(slices_path)
        return G_final_modules

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
        G_final_modules = self._run_domino_pipeline(G, seed_nodes)
        return G_final_modules
