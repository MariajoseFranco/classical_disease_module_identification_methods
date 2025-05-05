from collections import defaultdict

import numpy as np
from scipy.special import gammaln

# from scipy.stats import hypergeom


class DIAMOND:
    def __init__(self, alpha=1):
        self.alpha = alpha

    def compute_all_gamma_ln(self, N):
        return {i: gammaln(i) for i in range(1, N + 2)}

    def logchoose(self, n, k, gamma_ln):
        if n - k + 1 <= 0:
            return np.inf
        return gamma_ln[n + 1] - (gamma_ln[k + 1] + gamma_ln[n - k + 1])

    def gauss_hypergeom(self, x, r, b, n, gamma_ln):
        return np.exp(self.logchoose(r, x, gamma_ln) +
                      self.logchoose(b, n - x, gamma_ln) -
                      self.logchoose(r + b, n, gamma_ln))

    def pvalue(self, kb, k, N, s, gamma_ln):
        p = 0.0
        for n in range(kb, k + 1):
            if n > s:
                break
            prob = self.gauss_hypergeom(n, s, N - s, k, gamma_ln)
            p += prob
        return min(p, 1.0)

    def get_neighbors_and_degrees(self, G):
        neighbors = {node: set(G.neighbors(node)) for node in G.nodes()}
        degrees = dict(G.degree())
        return neighbors, degrees

    def reduce_candidates(self, degrees, neighbors, not_in_cluster, cluster_nodes, alpha):
        reduced = {}
        kb2k = defaultdict(dict)
        for node in not_in_cluster:
            k = degrees[node]
            kb = sum(1 for nbr in neighbors[node] if nbr in cluster_nodes)
            k += (alpha - 1) * kb
            kb += (alpha - 1) * kb
            kb2k[kb][k] = node

        k2kb = defaultdict(dict)
        for kb, k_map in kb2k.items():
            min_k = min(k_map.keys())
            k2kb[min_k][kb] = k_map[min_k]

        for k, kb_map in k2kb.items():
            max_kb = max(kb_map.keys())
            node = kb_map[max_kb]
            reduced[node] = (max_kb, k)

        return reduced

    def run_diamond(self, G, seed_nodes, num_iterations=100, save_path=None):
        G = G.copy()
        all_nodes = set(G.nodes())
        seed_nodes = set(seed_nodes) & all_nodes
        N = G.number_of_nodes()
        neighbors, degrees = self.get_neighbors_and_degrees(G)
        gamma_ln = self.compute_all_gamma_ln(N + 1)

        cluster = set(seed_nodes)
        s0 = len(cluster)
        s0 += (self.alpha - 1) * s0
        N += (self.alpha - 1) * s0

        not_in_cluster = set()
        for node in cluster:
            not_in_cluster |= neighbors[node]
        not_in_cluster -= cluster

        results = []
        all_p = {}

        while len(results) < num_iterations:
            candidates = self.reduce_candidates(
                degrees, neighbors, not_in_cluster, cluster, self.alpha
            )

            pmin, next_node, info = 10, None, {}
            for node, (kb, k) in candidates.items():
                key = (k, kb, s0)
                p = all_p.get(key)
                if p is None:
                    p = self.pvalue(kb, k, N, s0, gamma_ln)
                    all_p[key] = p
                if p < pmin:
                    pmin = p
                    next_node = node
                info[node] = (k, kb, p)

            if next_node is None:
                break

            results.append((next_node, *info[next_node]))
            cluster.add(next_node)
            s0 = len(cluster)
            not_in_cluster |= neighbors[next_node] - cluster
            not_in_cluster.remove(next_node)

        if save_path:
            with open(save_path, 'w') as fout:
                fout.write('#rank\tDIAMOnD_node\n')
                for i, (node, _, _, p) in enumerate(results, 1):
                    fout.write(f'{i}\t{node}\n')

        return results

# class DIAMOND():
#     def __init__(self) -> None:
#         pass

#     def compute_p_value(self, k, ks, N, s):
#         """
#         Compute the hypergeometric p-value of a node having ks links to s seed nodes
#         in a network of N total nodes and degree k.
#         """
#         return hypergeom.sf(ks - 1, N, s, k)

#     def run_diamond(self, G, seed_nodes, num_iterations=100):
#         """
#         Run the DIAMOnD algorithm on graph G starting from seed_nodes.
#         Returns a list of nodes ranked by disease relevance.
#         """
#         N = G.number_of_nodes()
#         all_nodes = set(G.nodes())
#         seed_set = set(seed_nodes)
#         diamond_nodes = []

#         for _ in range(num_iterations):
#             candidate_nodes = list(all_nodes - seed_set)
#             p_values = {}

#             for node in candidate_nodes:
#                 neighbors = set(G.neighbors(node))
#                 ks = len(neighbors & seed_set)
#                 k = len(neighbors)

#                 if ks == 0:
#                     continue  # skip if no connection to seed

#                 p_val = self.compute_p_value(k, ks, N, len(seed_set))
#                 p_values[node] = p_val

#             if not p_values:
#                 print("No more candidates with ks > 0.")
#                 break

#             # Add the node with lowest p-value to the module
#             next_node = min(p_values, key=p_values.get)
#             seed_set.add(next_node)
#             diamond_nodes.append((next_node, p_values[next_node]))
#         return diamond_nodes
