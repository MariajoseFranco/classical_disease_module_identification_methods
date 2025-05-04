import networkx as nx
import numpy as np


class ROBUST():
    def __init__(self) -> None:
        pass

    def run_robust_fallback(self, G, seeds, n_trees=30, alpha=0.25, beta=0.9, s_threshold=0.1):
        """
        Approximate ROBUST using shortest paths + consensus across trees.

        Parameters:
            G (nx.Graph): interaction graph
            seeds (list): seed nodes
            n_trees (int): how many consensus trees to build
            alpha (float): background score
            beta (float): decay factor for repeated nodes
            s_threshold (float): support threshold (e.g., 0.1 for 10%)
        """
        import random
        from collections import defaultdict

        # Step 1: Largest connected component
        G_lcc = max(nx.connected_components(G), key=len)
        G = G.subgraph(G_lcc).copy()

        # Step 2: Filter seeds
        seeds = [s for s in seeds if s in G]
        if not seeds:
            raise ValueError("No valid seed genes found in graph.")

        # Step 3: Assign artificial "prizes"
        prizes = defaultdict(lambda: alpha)
        for s in seeds:
            prizes[s] = 1.0

        # Step 4: Generate diverse subgraphs
        node_support = defaultdict(int)
        for _ in range(n_trees):
            sampled_seeds = random.sample(seeds, min(len(seeds), 5))

            # Build pairwise shortest paths
            selected_nodes = set(sampled_seeds)
            for i in range(len(sampled_seeds)):
                for j in range(i + 1, len(sampled_seeds)):
                    try:
                        path = nx.shortest_path(G, sampled_seeds[i], sampled_seeds[j])
                        selected_nodes.update(path)
                    except nx.NetworkXNoPath:
                        continue

            # Connect with a minimum spanning tree
            subG = G.subgraph(selected_nodes)
            if nx.is_connected(subG):
                mst = nx.minimum_spanning_tree(subG)
                for node in mst.nodes:
                    node_support[node] += 1

        # Step 5: Consensus threshold
        min_support = int(np.ceil(s_threshold * n_trees))
        final_nodes = [n for n, count in node_support.items() if count >= min_support]

        return G.subgraph(final_nodes).copy()
