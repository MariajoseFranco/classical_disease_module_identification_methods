import community as community_louvain
import networkx as nx
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


class DOMINO():
    def __init__(self) -> None:
        pass

    def domino(self, G, active_genes, min_size=5, fdr_threshold=0.3, bonferroni_threshold=0.05):
        """
        DOMINO algorithm implementation for active module identification.

        Parameters:
            G (networkx.Graph): PPI or gene interaction network
            active_genes (set or list): Set of 'active' gene IDs (e.g. significant genes)
            min_size (int): Minimum slice size to consider
            fdr_threshold (float): FDR q-value threshold for relevant slices
            bonferroni_threshold (float): Final q-value threshold for module significance

        Returns:
            modules (list of sets): Final DOMINO modules with enriched active genes
        """
        active_genes = set(active_genes)
        total_genes = set(G.nodes())
        active_in_network = active_genes.intersection(total_genes)

        # Step 0: Louvain community detection
        partition = community_louvain.best_partition(G)
        communities = {}
        for node, comm_id in partition.items():
            communities.setdefault(comm_id, set()).add(node)

        # Step 1: Identify relevant slices via hypergeometric test (FDR q < 0.3)
        relevant_slices = []
        pvals = []
        slices = []

        M = len(total_genes)  # population size
        K = len(active_in_network)  # number of active genes in network

        for nodes in communities.values():
            N = len(nodes)  # slice size
            x = len(active_in_network.intersection(nodes))  # active genes in this slice

            if N < min_size or x == 0:
                continue

            pval = hypergeom.sf(x - 1, M, K, N)
            pvals.append(pval)
            slices.append(nodes)

        # FDR correction
        if not pvals:
            print("No relevant slices found.")
            return []

        _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        for i, q in enumerate(qvals):
            if q < fdr_threshold:
                relevant_slices.append(slices[i])

        # Step 2: Partition relevant slices into modules (optional: split big ones)
        putative_modules = []
        for slc in relevant_slices:
            subG = G.subgraph(slc).copy()
            for cc in nx.connected_components(subG):
                if len(cc) >= min_size:
                    putative_modules.append(set(cc))

        # Step 3: Final module filtering via Bonferroni-corrected HG test
        final_modules = []
        pvals_final = []

        for module in putative_modules:
            N = len(module)
            x = len(active_in_network.intersection(module))
            pval = hypergeom.sf(x - 1, M, K, N)
            pvals_final.append(pval)

        _, qvals_final, _, _ = multipletests(pvals_final, method='bonferroni')
        for i, q in enumerate(qvals_final):
            if q < bonferroni_threshold:
                final_modules.append(putative_modules[i])

        for i, mod in enumerate(final_modules):
            print(f"Module {i+1}: {len(mod)} genes")

        return final_modules
