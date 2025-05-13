from typing import List

import networkx as nx

from classical_methods.lcc_algorithm import LCC


class UnitEdgeWeight:
    def __call__(self, u, v, d):
        return 1  # Assign unit weight to every edge


class PPIInstance:
    def __init__(self, ppi_graph: nx.Graph, terminals: List[str], edge_weights):
        self.graph = ppi_graph
        self.terminals = set(terminals)
        self.edge_weights = edge_weights
        self.assign_weights()

    def assign_weights(self):
        for u, v in self.graph.edges:
            self.graph[u][v]['weight'] = self.edge_weights(u, v, self.graph[u][v])


class ExpMinMaxDiverseSteinerTreeComputer:
    def __init__(self, initial_fraction=0.25, reduction_factor=0.9, initial_terminal_multiple=2):
        self.initial_fraction = initial_fraction
        self.reduction_factor = reduction_factor
        self.initial_terminal_multiple = initial_terminal_multiple

    def __call__(self, ppi_instance: PPIInstance, n: int = 30) -> List[nx.Graph]:
        # Dummy diverse steiner tree generator using shortest paths
        G = ppi_instance.graph
        terminals = list(ppi_instance.terminals)
        if not terminals:
            return []

        # Use a simple greedy path-growing method for baseline trees
        solutions = []
        for i in range(n):
            subgraph = nx.Graph()
            used = set()
            center = terminals[i % len(terminals)]
            used.add(center)
            subgraph.add_node(center)

            while len(used) < min(len(terminals), int(len(terminals) * self.initial_fraction)):
                dists = nx.single_source_dijkstra_path_length(G, center, weight='weight')
                candidates = [t for t in terminals if t not in used and t in dists]
                if not candidates:
                    break
                closest = min(candidates, key=lambda x: dists[x])
                path = nx.shortest_path(G, source=center, target=closest, weight='weight')
                nx.add_path(subgraph, path)
                used.update(path)

            solutions.append(subgraph)

        return solutions


class ROBUST:
    def __init__(self) -> None:
        self.lcc = LCC()

    def run_robust(self, G: nx.Graph, seed_nodes: List[str]) -> List[nx.Graph]:
        """
        Run the ROBUST method

        Args:
            G (nx.Graph): _description_
            seed_nodes (List[str]): _description_

        Returns:
            List[nx.Graph]: _description_
        """
        G_lcc = self.lcc.run_lcc(G)
        edge_weights = UnitEdgeWeight()
        ppi_instance = PPIInstance(ppi_graph=G_lcc, terminals=seed_nodes, edge_weights=edge_weights)

        tree_computer = ExpMinMaxDiverseSteinerTreeComputer(
            initial_fraction=0.25,
            reduction_factor=0.9,
            initial_terminal_multiple=2
        )

        solution_set = tree_computer(ppi_instance, n=30)
        return solution_set
