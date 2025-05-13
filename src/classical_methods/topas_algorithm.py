from functools import partial
from multiprocessing import Pool, cpu_count

import networkx as nx
import numpy as np
import pandas as pd

from classical_methods.lcc_algorithm import LCC


class TOPAS():
    def __init__(self, expansion_steps=2, cores=None):
        self.expansion_steps = expansion_steps
        self.cores = cores if cores is not None else cpu_count()
        self.lcc = LCC()

    def sp_compute(self, source_v, dist_matrix, seeds, graph, expansion_steps):
        d = dist_matrix.get(source_v, {})
        dest_v = [v for v, dist in d.items() if 1 < dist <= expansion_steps + 1 and v in seeds]
        paths = []
        for v in dest_v:
            try:
                path = nx.shortest_path(graph, source=source_v, target=v)
                paths.extend(path[1:-1])
            except nx.NetworkXNoPath:
                continue
        return paths

    def _compute_connectors(self, G_lcc, seeds):
        dist_matrix = dict(
            nx.all_pairs_shortest_path_length(G_lcc, cutoff=self.expansion_steps + 1)
        )

        with Pool(processes=self.cores) as pool:
            func = partial(
                self.sp_compute,
                dist_matrix=dist_matrix,
                seeds=seeds,
                graph=G_lcc,
                expansion_steps=self.expansion_steps
            )
            connector_lists = pool.map(func, seeds)

        connectors = set()
        for conn in connector_lists:
            connectors.update(conn)
        return connectors

    def _random_walk_prune(self, subgraph, seeds):
        A = nx.to_numpy_array(subgraph)
        M = A / A.sum(axis=1, keepdims=True)
        M = np.nan_to_num(M)

        r = np.array([1.0 if node in seeds else 0.0 for node in subgraph.nodes])
        r = r / r.sum()

        p = r.copy()
        for _ in range(100):
            p = (1 - 0.75) * np.matmul(M.T, p) + 0.75 * r

        df = pd.DataFrame({
            "node": list(subgraph.nodes),
            "val": [1 if node in seeds else 0 for node in subgraph.nodes],
            "prob": p
        }).sort_values("prob")

        for _, row in df[df["val"] == 0].iterrows():
            v = row["node"]
            if v not in subgraph:
                continue
            temp_subgraph = subgraph.copy()
            temp_subgraph.remove_node(v)
            if nx.number_connected_components(temp_subgraph) > 1:
                lcc_temp = self.lcc.run_lcc(temp_subgraph)
                if len(seeds & set(lcc_temp.nodes)) == len(seeds):
                    subgraph = lcc_temp
            else:
                subgraph = temp_subgraph

        return subgraph

    def read(self, network_file, seeds_file):
        network_df = pd.read_csv(network_file, sep="\t", header=None)
        seeds = pd.read_csv(seeds_file, sep="\t", header=None).iloc[:, 0].values
        return network_df, seeds

    def run(self, network_file, seeds_file):
        network_df, seeds = self.read(network_file, seeds_file)
        if network_df is None or network_df.shape[1] < 2:
            raise ValueError("Network file is missing.")
        if seeds is None or len(seeds) == 0:
            raise ValueError("Seeds are missing.")
        result = {
            'seed_nodes': list(seeds)
        }

        G = nx.Graph()
        G.add_edges_from(network_df.iloc[:, :2].values)
        seeds = set(map(str, seeds))

        G_lcc = self.lcc.run_lcc(G)
        seeds = seeds & set(G_lcc.nodes)
        connectors = self._compute_connectors(G_lcc, seeds)

        sub_nodes = seeds.union(connectors)
        subgraph = G_lcc.subgraph(sub_nodes).copy()

        if subgraph.number_of_edges() == 0:
            print("No module found!")
            return None

        subgraph = self.lcc.run_lcc(subgraph)
        seeds = seeds & set(subgraph.nodes)
        pruned_graph = self._random_walk_prune(subgraph, seeds)

        res = nx.to_pandas_edgelist(pruned_graph)
        sources = res['source']
        targets = res['target']
        all_nodes = pd.concat([sources, targets])
        unique_nodes = all_nodes.unique()

        result['seed_nodes_module_1'] = list(unique_nodes)

        return result
