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
        subgraph = G.subgraph(largest_cc).copy()
        nodes = subgraph.nodes
        result = {
            'seed_nodes': list(seed_nodes)
        }
        result['seed_nodes_module_1'] = list(nodes)
        return result
    def run_lcc_topas_style(self, G: nx.Graph, seed_nodes: set) -> nx.Graph:
        """
        Extracts the connected component of G that contains the largest number of seed nodes.
    
        Args:
            G: full graph
            seed_nodes: set of seed nodes

        Returns:
            Subgraph corresponding to the component with the most seed nodes.
        """
        if len(G) == 0:
            return G

        components = nx.connected_components(G)
        best_component = set()
        max_seed_count = 0

        for comp in components:
            comp = set(comp)
            seed_count = len(seed_nodes & comp)
            if seed_count > max_seed_count:
                max_seed_count = seed_count
                best_component = comp
    
        return G.subgraph(best_component).copy()

    def run_lcc(self, G: nx.Graph) -> nx.Graph:
        if len(G) == 0:
            return G
        largest_cc = max(nx.connected_components(G), key=len)
        return G.subgraph(largest_cc).copy()
