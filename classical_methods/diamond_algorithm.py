from scipy.stats import hypergeom


class DIAMOND():
    def __init__(self) -> None:
        pass

    def compute_p_value(self, k, ks, N, s):
        """
        Compute the hypergeometric p-value of a node having ks links to s seed nodes
        in a network of N total nodes and degree k.
        """
        return hypergeom.sf(ks - 1, N, s, k)

    def run_diamond(self, G, seed_nodes, num_iterations=100):
        """
        Run the DIAMOnD algorithm on graph G starting from seed_nodes.
        Returns a list of nodes ranked by disease relevance.
        """
        N = G.number_of_nodes()
        all_nodes = set(G.nodes())
        seed_set = set(seed_nodes)
        diamond_nodes = []

        for _ in range(num_iterations):
            candidate_nodes = list(all_nodes - seed_set)
            p_values = {}

            for node in candidate_nodes:
                neighbors = set(G.neighbors(node))
                ks = len(neighbors & seed_set)
                k = len(neighbors)

                if ks == 0:
                    continue  # skip if no connection to seed

                p_val = self.compute_p_value(k, ks, N, len(seed_set))
                p_values[node] = p_val

            if not p_values:
                print("No more candidates with ks > 0.")
                break

            # Add the node with lowest p-value to the module
            next_node = min(p_values, key=p_values.get)
            seed_set.add(next_node)
            diamond_nodes.append((next_node, p_values[next_node]))
        return diamond_nodes
