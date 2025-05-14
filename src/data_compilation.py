import pandas as pd
import numpy as np

class DataCompilation():
    def __init__(self, path, selected_diseases) -> None:
        self.path = path
        self.selected_diseases = selected_diseases

    def get_data(self):
        # Protein - Protein Interaction
        df_pro_pro = pd.read_csv(f'{self.path}pro_pro.tsv', sep='\t')
        df_pro_pro = df_pro_pro[df_pro_pro['prA'] != df_pro_pro['prB']]
        # Gen - Protein Interaction
        df_gen_pro = pd.read_csv(f'{self.path}gen_pro.tsv', sep='\t')
        # Disease - Gen Interaction
        df_dis_gen = pd.read_csv(f'{self.path}dis_gen.tsv', sep='\t')
        # Gene - Gene interaction, based on PPI
        
        return df_pro_pro, df_gen_pro, df_dis_gen

    def get_gen_gen_PPI(self, df_pro_pro, df_gen_pro) -> pd.DataFrame:

        genA = df_gen_pro.rename(columns={'protein_id':'prA', 'gene_id':'geneA'})
        genB = df_gen_pro.rename(columns={'protein_id':'prB', 'gene_id':'geneB'})
        
        merged = (
            df_pro_pro
            .merge(genA, on='prA', how='inner')
            .merge(genB, on='prB', how='inner')
        )
        
        df_gen_gen = merged[['geneA', 'geneB']].drop_duplicates()
        df_gen_gen = df_gen_gen[df_gen_gen['geneA'] != df_gen_gen['geneB']]
        
        df_gen_gen = (
            pd.DataFrame(
                np.sort(df_gen_gen.values, axis=1),
                columns=['geneA', 'geneB']
            )
            .drop_duplicates()
        )
        
        return df_gen_gen
    
    def get_dis_pro_data(self, df_dis_gen, df_gen_pro):
        df_dis_pro = df_dis_gen.merge(
            df_gen_pro, how='left', on='gene_id', indicator=True
        )
        df_dis_pro = df_dis_pro[df_dis_pro['_merge'] == 'both'].drop(
            ['_merge'], axis=1
        )
        return df_dis_pro

    def main(self):
        df_pro_pro, df_gen_pro, df_dis_gen = self.get_data()
        df_dis_pro = self.get_dis_pro_data(df_dis_gen, df_gen_pro)
        df_dis_pro = df_dis_pro[df_dis_pro['disease_name'].isin(
            self.selected_diseases
        )]
        df_gen_gen = self.get_gen_gen_PPI(df_pro_pro, df_gen_pro)
        return df_pro_pro, df_gen_pro, df_dis_gen, df_dis_pro, df_gen_gen
