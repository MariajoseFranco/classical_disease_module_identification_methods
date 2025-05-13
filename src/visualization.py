import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx


class VisualizationModule():
    def __init__(self):
        pass

    def visualize_modules(self, graph_dict, G, disease, method):
        os.makedirs("./src/outputs/imgs", exist_ok=True)

        seed_nodes = set(graph_dict.get("seed_nodes", []))

        # Build full subgraph with all nodes in all modules
        all_nodes = set()
        for k in graph_dict:
            all_nodes.update(graph_dict[k])
        subG = G.subgraph(all_nodes).copy()

        if len(subG.nodes) == 0:
            print("Nothing to visualize: combined subgraph is empty.")
            return

        # Assign a unique color to each module
        color_map = {}
        module_keys = [k for k in graph_dict if k != "seed_nodes"]
        color_list = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
        color_list = [c for c in color_list if c not in ['red', 'orange']]  # reserve red/orange

        # Assign module colors
        for i, key in enumerate(module_keys):
            color_map[key] = color_list[i % len(color_list)]

        # Set layout
        pos = nx.spring_layout(subG, seed=42, k=0.5)
        plt.figure(figsize=(12, 10))

        # Draw seed nodes
        seeds_in_subG = list(seed_nodes & subG.nodes)
        nx.draw_networkx_nodes(
            subG, pos,
            nodelist=seeds_in_subG,
            node_color='red',
            label='seed_nodes',
            node_size=300
        )

        # Draw each module
        for key, color in color_map.items():
            module_nodes = set(graph_dict[key])
            seed_in_module = module_nodes & seed_nodes
            connector_nodes = module_nodes - seed_in_module

            nx.draw_networkx_nodes(
                subG, pos,
                nodelist=list(module_nodes),
                node_color=color,
                label=key,
                node_size=300,
                alpha=0.8
            )

        # Edges and labels
        nx.draw_networkx_edges(subG, pos, alpha=0.3)
        nx.draw_networkx_labels(subG, pos, font_size=7)

        plt.title(f"{method.upper()} Modules - {disease.title()}")
        plt.axis('off')
        plt.legend()
        plt.savefig(f"./src/outputs/imgs/{method}_combined_graph_{disease}.png")
        plt.show()

    def visualize_seed_gene_subgraph(self, disease, G_ppi, disease_gene_map):
        os.makedirs("./src/outputs/imgs", exist_ok=True)
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

        # Create a figure and axes explicitly
        fig, ax = plt.subplots(figsize=(10, 8))

        # Draw using the given axes
        pos = nx.spring_layout(G_seed_sub)
        nx.draw(G_seed_sub, pos, ax=ax, with_labels=True, node_color='red', node_size=300)

        ax.set_title(f"Seed Gene Subgraph for {disease.title()}")
        fig.tight_layout()

        plt.savefig(f"./src/outputs/imgs/{disease}_seed_graph.png")
        plt.show()
