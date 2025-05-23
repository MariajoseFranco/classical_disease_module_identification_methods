import pandas as pd
from tqdm import tqdm

from classical_methods.diamond_algorithm import DIAMOND
from classical_methods.domino_algorithm import DOMINO
from classical_methods.lcc_algorithm import LCC
from classical_methods.robust_ant import ROBUST
from classical_methods.topas_algorithm import TOPAS
from data_compilation import DataCompilation
from graph_creation import GraphPPI
from visualization import VisualizationModule


class Main():
    def __init__(self, path):
        # Select the diseases to work with
        self.selected_diseases = ["Albinism", "Alcohol Use Disorder"]
        self.DC = DataCompilation(path, self.selected_diseases)
        self.GPPI = GraphPPI()
        self.V = VisualizationModule()
        self.LCC = LCC()
        self.DIAMOND = DIAMOND()
        self.DOMINO = DOMINO()
        self.ROBUST = ROBUST()
        self.TOPAS = TOPAS(expansion_steps=2, cores=4)

    def run_classical_methods(self, G_ppi, disease_pro_mapping, MIN_SEEDS=10):
        results = {}

        for disease, all_seeds in tqdm(disease_pro_mapping.items()):
            print(f"Processing: {disease} ({len(all_seeds)} raw seeds)")

            with open(f'./src/inputs/seed_nodes_{disease}.txt', 'w') as f:
                for seed in all_seeds:
                    f.write(f"{seed}\n")

            seed_nodes = [node for node in all_seeds if node in G_ppi]
            if len(seed_nodes) < MIN_SEEDS:
                print("Skipped — not enough seeds in PPI")
                continue

            results[disease] = {}

            try:
                print("\nRunning LCC...")
                results[disease]["lcc"] = self.LCC.run_lcc_per_disease(
                    G_ppi,
                    seed_nodes
                )
                print("...done")
            except Exception as e:
                print("LCC failed:", e)

            try:
                print("\nRunning TOPAS...")
                results[disease]["topas"] = self.TOPAS.run(
                    "./src/inputs/PPI.txt",
                    f"./src/inputs/seed_nodes_{disease}.txt"
                )
                print("...done")
            except Exception as e:
                print("TOPAS failed:", e)

            try:
                print("\nRunning DIAMOND...")
                results[disease]["diamond"] = self.DIAMOND.run_diamond(
                    "./src/inputs/PPI.txt",
                    f"./src/inputs/seed_nodes_{disease}.txt",
                    200
                )
                print("...done")
            except Exception as e:
                print("DIAMOnD failed:", e)

            try:
                print("\nRunning DOMINO...")
                results[disease]["domino"] = self.DOMINO.run_domino(
                    "./src/inputs/PPI.txt",
                    f"./src/inputs/seed_nodes_{disease}.txt"
                )
                print("...done")
            except Exception as e:
                print("DOMINO failed:", e)

            try:
                print("\nRunning ROBUST...")
                out_csv = f"./src/outputs/robust_{disease}.csv"
                results[disease]["robust"] = self.ROBUST.run_robust(
                    self.df_gen_gen,
                    seed_nodes,
                    out_csv
                )
                print("...done")
            except Exception as e:
                print("ROBUST failed:", e)
        return results

    def save_classical_methods_results(self, results):
        all_modules = []

        for disease, methods in results.items():
            for method, module_dict in methods.items():
                for module, genes in module_dict.items():
                    all_modules.append({
                        "disease": disease,
                        "method": method,
                        "module": module,
                        "nodes": genes
                    })

        pd.DataFrame(all_modules).to_csv("./src/outputs/multi_disease_modules.csv", index=False)
        print("Saved all modules to 'multi_disease_modules.csv'")

    def visualize_disease_results(
            self, disease, G_ppi, disease_pro_mapping, results
    ):
        # Seed Gene Subgraph
        self.V.visualize_seed_gene_subgraph(
            disease, G_ppi, disease_pro_mapping
        )

        # Methods
        self.V.visualize_modules(results[disease]["lcc"], G_ppi, disease, "lcc")
        self.V.visualize_modules(results[disease]["diamond"], G_ppi, disease, "diamond")
        self.V.visualize_modules(results[disease]["domino"], G_ppi, disease, "domino")
        self.V.visualize_modules(results[disease]["topas"], G_ppi, disease, "topas")

    def main(self):
        # Classical Methods
        df_pro_pro, df_gen_pro, df_dis_gen, df_dis_pro, df_gen_gen = self.DC.main()
        self.df_gen_gen = df_gen_gen
        df_pro_pro.to_csv("./src/inputs/PPI.txt", sep="\t", index=False)
        G_ppi, disease_pro_mapping = self.GPPI.main(df_pro_pro, df_dis_pro)
        results_classical_methods = self.run_classical_methods(G_ppi, disease_pro_mapping)
        self.save_classical_methods_results(results_classical_methods)
        for disease in self.selected_diseases:
            self.visualize_disease_results(
                disease, G_ppi, disease_pro_mapping, results_classical_methods
            )


if __name__ == "__main__":
    path = "./src/data/"
    # path = "/app/data/"
    Main(path).main()
