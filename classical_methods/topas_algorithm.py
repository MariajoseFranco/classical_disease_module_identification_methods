import networkx as nx
import numpy as np


class TOPAS():
    def __init__(self) -> None:
        pass

    def random_walk_with_restart(self, G, seeds, restart_prob=0.75, max_iter=100, tol=1e-6):
        nodes = list(G.nodes())
        idx_map = {node: i for i, node in enumerate(nodes)}
        n = len(nodes)

        A = nx.to_numpy_array(G, nodelist=nodes)
        A = A / A.sum(axis=1, keepdims=True)  # row-normalized

        # Initial probabilities (1 for seed nodes)
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

    def find_potential_connectors(self, G, seed_nodes, max_dist=3):
        connectors = set()
        for i, s1 in enumerate(seed_nodes):
            for s2 in seed_nodes[i+1:]:
                try:
                    path = nx.shortest_path(G, s1, s2)
                    if 2 <= len(path) - 1 <= max_dist:
                        connectors.update(set(path[1:-1]))  # exclude the seeds themselves
                except nx.NetworkXNoPath:
                    continue
        return list(connectors - set(seed_nodes))

    def topas_with_top_percent_rwr(self, G, seed_nodes, max_dist=3, top_percent=0.1):
        """
        TOPAS variant: keep top X% RWR-ranked connector nodes.
        """
        G = G.copy()
        seed_nodes = [s for s in seed_nodes if s in G.nodes]  # filter invalid seeds

        # Step 1: Find potential connectors
        connectors = self.find_potential_connectors(G, seed_nodes, max_dist)
        print(f"Found {len(connectors)} connector candidates")

        if not connectors:
            print("No connector candidates found.")
            return G.subgraph(seed_nodes).copy()

        # Step 2: Rank by RWR scores
        rwr_scores = self.random_walk_with_restart(G, seed_nodes)
        ranked_connectors = sorted(connectors, key=lambda x: -rwr_scores.get(x, 0))

        # Step 3: Keep top X% of connectors
        n_top = max(1, int(top_percent * len(ranked_connectors)))
        top_connectors = ranked_connectors[:n_top]
        print(f"Retaining top {n_top} connectors ({top_percent*100:.0f}%) by RWR score")

        # Step 4: Build module
        module_nodes = set(seed_nodes) | set(top_connectors)
        Gs = G.subgraph(module_nodes).copy()

        return Gs
