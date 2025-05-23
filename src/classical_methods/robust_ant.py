import os
import subprocess
import sys
import tempfile


class ROBUST:
    def __init__(self,
                 alpha=0.25, beta=0.5, n=50, tau=0.5,
                 study_bias_scores='BAIT_USAGE', gamma=1.0,
                 namespace='ENTREZ',
                 robust_script_path='./state_of_art_repos/robust_bias_aware/robust.py'):
        self.alpha = alpha
        self.beta = beta
        self.n = n
        self.tau = tau
        self.study_bias_scores = study_bias_scores
        self.gamma = gamma
        self.namespace = namespace
        self.robust_script = os.path.abspath(robust_script_path)

    def run_robust(self, df_gen_gen, seed_nodes, out_csv):
        """
        df_gen_gen : pd.DataFrame with columns ['geneA','geneB']
        seed_nodes : list de ENTREZ IDs (strings)
        out_csv    : path to save the .csv output
        """

        # Creation of temporary files for df_gen_gen and the seed genes
        with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as tmp_ppi:
            df_gen_gen.to_csv(tmp_ppi.name, sep='\t', index=False, header=False)
            ppi_path = tmp_ppi.name

        with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as tmp_seeds:
            for s in seed_nodes:
                tmp_seeds.write(f"{s}\n")
            seeds_path = tmp_seeds.name

        # Creation of the ROBUST command
        cmd = [
            sys.executable, os.path.abspath(self.robust_script),
            seeds_path, out_csv,  # los dos argumentos posicionales requeridos
            '--network',    ppi_path,
            '--namespace',  self.namespace,
            '--alpha',      str(self.alpha),
            '--beta',       str(self.beta),
            '--n',          str(self.n),
            '--tau',        str(self.tau),
            '--study-bias-scores', self.study_bias_scores,
            '--gamma',      str(self.gamma)
        ]
        subprocess.run(cmd, check=True)

        # Removal of temporary files
        os.remove(ppi_path)
        os.remove(seeds_path)
        return out_csv
