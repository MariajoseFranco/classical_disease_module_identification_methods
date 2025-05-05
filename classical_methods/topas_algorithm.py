import networkx as nx
import numpy as np


class TOPAS:
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

    def extract_lcc(self, G, seeds):
        components = list(nx.connected_components(G))
        best_cc = max(components, key=lambda comp: len(set(comp) & set(seeds)))
        return G.subgraph(best_cc).copy()

    def find_potential_connectors(self, G, seed_nodes, max_dist=3):
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

    def topas_full(self, G, seed_nodes, max_dist=3, top_percent=0.3):
        G = G.copy()
        seed_nodes = [s for s in seed_nodes if s in G.nodes]
        if len(seed_nodes) == 0:
            return nx.Graph()

        # Step 1: LCC extraction
        G_lcc = self.extract_lcc(G.subgraph(seed_nodes + list(G.nodes())), seed_nodes)
        seed_nodes = [s for s in seed_nodes if s in G_lcc.nodes]

        # Step 2: Connector discovery
        connectors = self.find_potential_connectors(G_lcc, seed_nodes, max_dist)
        if not connectors:
            return G_lcc.subgraph(seed_nodes).copy()

        # Step 3: RWR scoring
        rwr_scores = self.random_walk_with_restart(G_lcc, seed_nodes)
        ranked_connectors = sorted(connectors, key=lambda x: -rwr_scores.get(x, 0))
        n_top = max(1, int(top_percent * len(ranked_connectors)))
        top_connectors = ranked_connectors[:n_top]

        # Step 4: Build and prune module
        module_nodes = set(seed_nodes) | set(top_connectors)
        G_module = G_lcc.subgraph(module_nodes).copy()
        G_pruned = self.prune_module(G_module, seed_nodes, rwr_scores)

        return G_pruned
