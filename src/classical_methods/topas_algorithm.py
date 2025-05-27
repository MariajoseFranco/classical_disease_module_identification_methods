from multiprocessing import Pool, cpu_count

import networkx as nx
import numpy as np
import pandas as pd

from classical_methods.lcc_algorithm import LCC

# Globals for worker processes
global G_shared, seeds_shared, exp_shared
G_shared = None
seeds_shared = None
exp_shared = None


def init_worker(graph, seeds, expansion_steps):
    """
    Initializer for worker subprocesses: sets up shared globals.
    """
    global G_shared, seeds_shared, exp_shared
    G_shared = graph
    seeds_shared = seeds
    exp_shared = expansion_steps


def sp_compute_global(source_v):
    """
    Computes intermediate nodes on shortest paths from a seed to other seeds
    within expansion_steps + 1, without building the full all-pairs matrix.
    Uses globals G_shared, seeds_shared, exp_shared.
    """
    # Single-source distances up to cutoff
    lengths = nx.single_source_shortest_path_length(
        G_shared,
        source=source_v,
        cutoff=exp_shared + 1
    )
    # Destinations among seeds at appropriate distance
    dest_v = [v for v, d in lengths.items()
              if 1 < d <= exp_shared + 1 and v in seeds_shared]

    paths = []
    # Gather interior nodes of each shortest path
    for v in dest_v:
        try:
            path = nx.shortest_path(G_shared, source=source_v, target=v)
            paths.extend(path[1:-1])
        except nx.NetworkXNoPath:
            continue
    return paths


class TOPAS:
    def __init__(self, expansion_steps=2, cores=None, tol=1e-6, restart_prob=0.75, max_iter=1000):
        self.expansion_steps = expansion_steps
        self.cores = cores if cores is not None else cpu_count()
        if not isinstance(self.cores, int) or self.cores < 1 or self.cores > cpu_count():
            raise ValueError(f"'cores' must be an integer between 1 and {cpu_count()}")
        self.restart_prob = restart_prob
        self.max_iter = max_iter
        self.tol = tol
        self.lcc = LCC()

    def _compute_connectors(self, G_lcc, seeds):
        with Pool(
            processes=self.cores,
            initializer=init_worker,
            initargs=(G_lcc, seeds, self.expansion_steps)
        ) as pool:
            connector_lists = pool.map(sp_compute_global, list(seeds))

        connectors = set()
        for conn in connector_lists:
            connectors.update(conn)
        return connectors

    def _random_walk_prune(self, subgraph, seeds):
        # transition matrix row-normalized
        A = nx.to_numpy_array(subgraph)
        row_sums = A.sum(axis=1, keepdims=True)
        M = np.divide(A, row_sums, where=row_sums != 0)

        # restart vector (un-normalized)
        r = np.array([1.0 if node in seeds else 0.0 for node in subgraph.nodes])
        # initial prob vector
        p = r.copy()

        # iterate until convergence
        for _ in range(self.max_iter):
            p_new = (1 - self.restart_prob) * (M.T @ p) + self.restart_prob * r
            if np.linalg.norm(p_new - p, 1) < self.tol:
                p = p_new
                break
            p = p_new

        df = pd.DataFrame({
            "node": list(subgraph.nodes),
            "is_seed": [1 if node in seeds else 0 for node in subgraph.nodes],
            "prob": p
        }).sort_values("prob")

        # iterative pruning
        for _, row in df[df["is_seed"] == 0].iterrows():
            v = row["node"]
            if v not in subgraph:
                continue
            temp = subgraph.copy()
            temp.remove_node(v)
            if nx.number_connected_components(temp) > 1:
                # extract LCC based on seed count
                comp = self.lcc.run_lcc(temp)
                remaining = seeds & set(comp.nodes)
                if len(remaining) == len(seeds):
                    subgraph = comp
                    seeds = remaining
            else:
                subgraph = temp
        return subgraph

    def read(self, network_file, seeds_file):
        network_df = pd.read_csv(network_file, sep="\t", header=None)
        seeds = pd.read_csv(seeds_file, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
        return network_df, set(seeds)

    def run(self, network_file, seeds_file):
        network_df, seeds = self.read(network_file, seeds_file)
        if network_df is None or network_df.shape[1] < 2:
            raise ValueError("Network file is missing.")
        if not seeds:
            raise ValueError("Seeds are missing.")

        # build graph and simplify
        G = nx.Graph()
        G.add_edges_from(network_df.iloc[:, :2].astype(str).values)
        G.remove_edges_from(nx.selfloop_edges(G))

        # STEP 1: full seed network
        # STEP 2: largest connected module (LCC)
        G_lcc = self.lcc.run_lcc_topas_style(G, seeds)
        seeds &= set(G_lcc.nodes)

        # STEP 2b: compute connectors
        connectors = self._compute_connectors(G_lcc, seeds)

        # extract subgraph of seeds + connectors
        sub_nodes = seeds.union(connectors)
        subgraph = G_lcc.subgraph(sub_nodes).copy()
        if subgraph.number_of_edges() == 0:
            print("No module found!")
            return None

        # ensure largest seed-based CC
        subgraph = self.lcc.run_lcc(subgraph)
        seeds &= set(subgraph.nodes)

        # STEP 3: random walk pruning
        pruned = self._random_walk_prune(subgraph, seeds)

        # format results
        edgelist = nx.to_pandas_edgelist(pruned)
        nodes = pd.unique(edgelist[['source', 'target']].values.ravel())
        return {
            'edgelist': edgelist,
            'module_nodes': list(nodes)
        }
