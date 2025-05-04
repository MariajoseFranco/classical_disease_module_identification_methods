import matplotlib.pyplot as plt
import networkx as nx


class VisualizationModule():
    def __init__(self):
        pass

    def visualize_module(self, G, node_list, title="Module", seed_nodes=None):
        subG = G.subgraph(node_list).copy()
        if len(subG.nodes) == 0:
            print("Nothing to visualize: subgraph is empty.")
            return

        pos = nx.spring_layout(subG, seed=42, k=0.5)  # Increase node spacing

        plt.figure(figsize=(10, 8))

        if seed_nodes:
            seed_nodes = set(seed_nodes) & set(subG.nodes)

            connector_nodes = list(set(subG.nodes) - set(seed_nodes))

            nx.draw_networkx_nodes(
                subG,
                pos,
                nodelist=connector_nodes,
                node_color='orange',
                label='Connectors',
                node_size=300
            )
            nx.draw_networkx_nodes(
                subG,
                pos,
                nodelist=list(seed_nodes),
                node_color='red',
                label='Seed genes',
                node_size=300
            )
        else:
            nx.draw_networkx_nodes(
                subG,
                pos,
                node_color="skyblue",
                node_size=300
            )

        nx.draw_networkx_edges(
            subG,
            pos,
            alpha=0.4
        )
        nx.draw_networkx_labels(
            subG,
            pos,
            font_size=7
        )

        handles, labels = plt.gca().get_legend_handles_labels()
        if handles:
            plt.legend()
        plt.title(title)
        plt.axis('off')
        plt.show()

    def visualize_diamond_module(self, G, diamond_genes, disease_genes, disease):
        # Extract nodes in the final disease module
        module_nodes = set(gene for gene, _ in diamond_genes) | set(disease_genes)
        G_module = G.subgraph(module_nodes)

        # Filter disease_genes to those present in the current subgraph
        disease_genes_in_module = [gene for gene in disease_genes if gene in G_module.nodes]

        pos = nx.spring_layout(G_module)

        # Separate DIAMOnD predictions
        diamond_only_nodes = [
            node for node in G_module.nodes if node not in disease_genes_in_module
        ]

        # Draw everything
        plt.figure(figsize=(10, 8))
        nx.draw_networkx_nodes(
            G_module,
            pos,
            nodelist=diamond_only_nodes,
            node_size=80,
            node_color='orange',
            label='Predicted'
        )
        nx.draw_networkx_nodes(
            G_module,
            pos,
            nodelist=disease_genes_in_module,
            node_size=80,
            node_color='red',
            label='Seed Genes'
        )
        nx.draw_networkx_edges(
            G_module,
            pos,
            alpha=0.5
        )
        nx.draw_networkx_labels(
            G_module,
            pos,
            font_size=6
        )

        plt.title(f"DIAMOnD Disease Module for {disease}")
        plt.legend()
        plt.axis("off")
        plt.show()

    def visualize_domino_module(self, G, domino_modules, disease_genes, disease):
        filtered_disease_genes = [g for g in disease_genes if g in G.nodes()]
        try:
            mod = domino_modules[0]
            subG = G.subgraph(mod)
            seed_nodes = set(filtered_disease_genes)
            predicted_nodes = set(mod) - seed_nodes

            plt.figure(figsize=(10, 8))
            pos = nx.spring_layout(subG, seed=42)

            nx.draw_networkx_nodes(
                subG,
                pos,
                nodelist=list(predicted_nodes),
                node_color='orange',
                label="Predicted"
            )
            nx.draw_networkx_nodes(
                subG,
                pos,
                nodelist=list(seed_nodes & set(subG.nodes())),
                node_color='red',
                label="Seed"
            )
            nx.draw_networkx_edges(
                subG,
                pos,
                alpha=0.3
            )
            nx.draw_networkx_labels(
                subG,
                pos,
                font_size=6
            )

            plt.title(f"DOMINO: Module 1 for {disease}")
            plt.legend()
            plt.axis('off')
            plt.show()
        except Exception as ex:
            print(f"No DOMINO modules found: {ex}")

    def visualize_seed_gene_subgraph(self, disease, G_ppi, disease_gene_map):
        disease_seeds = set(disease_gene_map[disease])
        seeds_in_graph = disease_seeds & set(G_ppi.nodes())
        print(f"Total {disease.title()} seed genes: {len(disease_seeds)}")
        print(f"Seed genes present in G_ppi: {len(seeds_in_graph)}")

        G_seed_sub = G_ppi.subgraph(seeds_in_graph)
        print(
            "Seed subgraph size:",
            G_seed_sub.number_of_nodes(),
            "nodes,",
            G_seed_sub.number_of_edges(),
            "edges"
        )

        nx.draw(G_seed_sub, with_labels=True, node_color='red', node_size=300)
        plt.title(f"Seed Gene Subgraph for {disease.title()}")
        plt.show()
