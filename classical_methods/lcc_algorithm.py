import networkx as nx


class LCC():
    def __init__(self) -> None:
        pass

    def run_lcc(self, G, disease_genes):
        subgraph = G.subgraph([g for g in disease_genes if g in G.nodes()])
        components = list(nx.connected_components(subgraph))
        if not components:
            return nx.Graph()
        largest_cc = max(components, key=len)
        return G.subgraph(largest_cc).copy()
