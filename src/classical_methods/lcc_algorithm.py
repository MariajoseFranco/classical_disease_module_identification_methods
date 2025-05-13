import networkx as nx
from networkx import Graph


class LCC():
    def __init__(self) -> None:
        pass

    def run_lcc_per_disease(self, G: Graph, seed_nodes: list) -> Graph:
        """
        Run the LCC method

        Args:
            G: the graph with the protein-protein interaction
            seed_nodes: the seed nodes of the disease of interest

        Returns:
            Graph: subgraph containing only the nodes that are
            present in the largest_cc
        """
        subgraph = G.subgraph([node for node in seed_nodes if node in G.nodes()])
        components = list(nx.connected_components(subgraph))
        if not components:
            return nx.Graph()
        largest_cc = max(components, key=len)
        return G.subgraph(largest_cc).copy()

    def run_lcc(self, G: nx.Graph) -> nx.Graph:
        if len(G) == 0:
            return G
        largest_cc = max(nx.connected_components(G), key=len)
        return G.subgraph(largest_cc).copy()
