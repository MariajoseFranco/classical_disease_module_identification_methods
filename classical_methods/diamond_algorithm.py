from collections import defaultdict
from typing import Union

import numpy as np
from networkx import Graph
from scipy.special import gammaln


class DIAMOND:
    def __init__(self, alpha=1):
        self.alpha = alpha

    def compute_all_gamma_ln(self, N: int) -> dict:
        """
        Computes the logarithm of the gamma function for all integers from 1 to N + 1

        Args:
            N: total number of nodes in G

        Returns:
            dict: dictionary with each node and its logarithm of the gamma function
        """
        return {i: gammaln(i) for i in range(1, N + 2)}

    def logchoose(self, n: int, k: int, gamma_ln: dict) -> float:
        """
        Computes the natural logarithm of the binomial coefficient "n choose k"
        using precomputed logarithmic gamma values

        Args:
            n: Total number of elements
            k: Number of selected elements
            gamma_ln: Dictionary mapping integers to their log-gamma values

        Returns:
            float: The logarithm of the binomial coefficient C(n, k). Returns infinity if k > n.
        """
        if n - k + 1 <= 0:
            return np.inf
        return gamma_ln[n + 1] - (gamma_ln[k + 1] + gamma_ln[n - k + 1])

    def gauss_hypergeom(self, x: int, r: int, b: int, n: int, gamma_ln: dict) -> float:
        """
        Computes the value of the hypergeometric probability mass function (PMF)
        using logarithmic gamma functions for numerical stability.

        This represents the probability of drawing exactly 'x' successful outcomes
        (seed genes) in 'n' draws, without replacement, from a population with:
            - 'r' total successes (number of seed genes),
            - 'b' total failures (non-seed nodes),
            - 'r + b' total population size.

        Args:
            x: Number of observed successes (edges to seed nodes).
            r: Total number of possible successes in the population (seed nodes).
            b: Total number of failures in the population (non-seed nodes).
            n: Number of draws (node degree).
            gamma_ln: Precomputed log-gamma values for efficient binomial coefficient calculation.

        Returns:
            float: Hypergeometric probability value for observing 'x' successes.
        """
        return np.exp(self.logchoose(r, x, gamma_ln) +
                      self.logchoose(b, n - x, gamma_ln) -
                      self.logchoose(r + b, n, gamma_ln))

    def pvalue(self, kb: int, k: int, N: int, s: int, gamma_ln: dict) -> float:
        """
        Computes the p-value of observing at least 'kb' connections between a candidate node
        and the current seed set (cluster), given its total degree 'k'. A lower p-value means
        the candidate node is more likely to be significantly connected to the current
        disease module

        Args:
            kb: Number of edges the candidate node has to the current seed set
            k: Total degree of the candidate node (possibly adjusted with alpha)
            N: Total number of nodes in the graph (adjusted for alpha-weighted seeds)
            s: Number of seed nodes currently in the cluster (adjusted).
            gamma_ln: Precomputed dictionary of log-gamma values to speed up calculations.

        Returns:
            float: The p-value of observing 'kb' or more links to the cluster, under a
            hypergeometric null model. The value is capped at 1.0.
        """
        p = 0.0
        for n in range(kb, k + 1):
            if n > s:
                break
            prob = self.gauss_hypergeom(n, s, N - s, k, gamma_ln)
            p += prob
        return min(p, 1.0)

    def get_neighbors_and_degrees(self, G: Graph) -> Union[dict, dict]:
        """

        Args:
            G: the graph with the protein-protein interaction

        Returns:
            dict: dictionary with adjacent nodes of eahc node
            dict: dictionary containing the mapping of each node to its degree
        """
        neighbors = {node: set(G.neighbors(node)) for node in G.nodes()}
        degrees = dict(G.degree())
        return neighbors, degrees

    def reduce_candidates(
            self,
            degrees: dict,
            neighbors: dict,
            not_in_cluster: set,
            cluster_nodes: set,
            alpha: int
    ) -> dict:
        """
        Selects a reduced subset of candidate nodes (not in the cluster) by:
            - Calculating how well each candidate connects to the current cluster
            - Ranking them based on their connectivity to the cluster (`kb`)
              and their total degree (`k`)
            - Adjusting both values using a seed-weight factor (`alpha`)
            - Selecting the best node per (kb, k) combination to reduce redundant
              computations in DIAMOnD

        Args:
            degrees: A mapping from each node to its degree in the graph
            neighbors: A mapping from each node to the set of its neighboring nodes
            not_in_cluster: Nodes that are candidates to be added to the cluster (not yet in it)
            cluster_nodes: The current set of nodes in the disease module cluster
            alpha: Weight multiplier applied to edges connecting to seed nodes,
            emphasizing their importance

        Returns:
            dict: A dictionary mapping each selected candidate node to a tuple (kb, k),
                where 'kb' is the adjusted count of neighbors in the cluster and
                'k' is the adjusted total degree, used for downstream ranking
        """
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

    def run_diamond(self, G: Graph, seed_nodes: list, num_iterations: int = 100) -> list:
        """
        Run the DIAMOND method

        Args:
            G: the graph with the protein-protein interaction
            seed_nodes: the seed nodes of the disease of interest
            num_iterations: desired number of nodes to be added. Defaults to 100.

        Returns:
            list: results obtained
        """
        seed_nodes = set(seed_nodes)
        N = G.number_of_nodes()
        neighbors, degrees = self.get_neighbors_and_degrees(G)
        gamma_ln = self.compute_all_gamma_ln(N + 1)

        cluster = seed_nodes.copy()
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

        return results
