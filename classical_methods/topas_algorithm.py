import networkx as nx
import numpy as np
from networkx import Graph

from classical_methods.lcc_algorithm import LCC


class TOPAS:
    def __init__(self) -> None:
        self.lcc = LCC()

    def random_walk_with_restart(self, G, seeds, restart_prob=0.75, max_iter=100, tol=1e-6):
        nodes = list(G.nodes())
        idx_map = {node: i for i, node in enumerate(nodes)}
        n = len(nodes)

        A = nx.to_numpy_array(G, nodelist=nodes)
        A = A / A.sum(axis=1, keepdims=True)

        p0 = np.zeros(n)
        for seed in seeds:
            if seed in idx_map:
                p0[idx_map[seed]] = 1.0
        p0 /= p0.sum()

        pt = p0.copy()
        for _ in range(max_iter):
            pt_next = (1 - restart_prob) * A @ pt + restart_prob * p0
            if np.linalg.norm(pt_next - pt, 1) < tol:
                break
            pt = pt_next

        return dict(zip(nodes, pt))

    def find_potential_connectors(self, G: Graph, seed_nodes: list, max_dist: int = 3) -> list:
        """_summary_

        Args:
            G: the graph with the protein-protein interaction
            seed_nodes: the seed nodes of the disease of interest
            max_dist: how far to search from the disease-related seed nodes. Defaults to 3.

        Returns:
            list: list of the connectors that are not present in the seed nodes
        """
        connectors = set()
        for i, s1 in enumerate(seed_nodes):
            for s2 in seed_nodes[i+1:]:
                try:
                    path = nx.shortest_path(G, s1, s2)
                    if 2 <= len(path) - 1 <= max_dist:
                        connectors.update(set(path[1:-1]))
                except nx.NetworkXNoPath:
                    continue
        return list(connectors - set(seed_nodes))

    def prune_module(self, G, seeds, rwr_scores):
        non_seeds = [n for n in G.nodes() if n not in seeds]
        sorted_by_score = sorted(non_seeds, key=lambda x: rwr_scores.get(x, 0))

        for node in sorted_by_score:
            if node not in G:
                continue
            G_tmp = G.copy()
            G_tmp.remove_node(node)
            components = list(nx.connected_components(G_tmp))
            best_cc = max(components, key=lambda comp: len(set(comp) & set(seeds)))
            if len(set(seeds) & best_cc) == len(seeds):
                G = G_tmp.subgraph(best_cc).copy()
        return G

    def run_topas(
            self, G: Graph, seed_nodes: list, max_dist: int = 3, top_percent: float = 0.3
    ) -> Graph:
        """
        Run the TOPAS method

        Args:
            G: the graph with the protein-protein interaction
            seed_nodes: the seed nodes of the disease of interest
            max_dist: how far to search from the disease-related seed nodes
            top_percent: how many of the nodes found should be retained for the final module

        Returns:
            _type_: _description_
        """
        # Step 1: LCC extraction
        G_lcc = self.lcc.run_lcc_per_disease(G, seed_nodes)
        seed_nodes_lcc = list(G_lcc.nodes)

        # Step 2: Connector discovery
        connectors = self.find_potential_connectors(G, seed_nodes_lcc, max_dist)  # G o G_lcc?
        if not connectors:
            return G_lcc.subgraph(seed_nodes_lcc).copy()

        # Step 3: RWR scoring
        rwr_scores = self.random_walk_with_restart(G_lcc, seed_nodes_lcc)
        ranked_connectors = sorted(connectors, key=lambda x: -rwr_scores.get(x, 0))
        n_top = max(1, int(top_percent * len(ranked_connectors)))
        top_connectors = ranked_connectors[:n_top]

        # Step 4: Build and prune module
        module_nodes = set(seed_nodes_lcc) | set(top_connectors)
        G_module = G_lcc.subgraph(module_nodes).copy()
        G_pruned = self.prune_module(G_module, seed_nodes_lcc, rwr_scores)

        return G_pruned
