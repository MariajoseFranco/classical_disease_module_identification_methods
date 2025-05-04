import networkx as nx


class GraphPPI():
    def __init__(self):
        pass

    def create_graph(self, df_pro_pro):
        # Create a graph from the PPI data
        G_ppi = nx.from_pandas_edgelist(df_pro_pro, 'prA', 'prB')
        print(f"PPI Network: {G_ppi.number_of_nodes()} nodes, {G_ppi.number_of_edges()} edges")
        return G_ppi

    def map_dis_gen(self, df_dis_pro):
        disease_gene_map = df_dis_pro.groupby("disease_name")["protein_id"].apply(set).to_dict()
        return disease_gene_map

    def main(self, df_pro_pro, df_dis_pro):
        G_ppi = self.create_graph(df_pro_pro)
        disease_gene_map = self.map_dis_gen(df_dis_pro)
        return G_ppi, disease_gene_map
