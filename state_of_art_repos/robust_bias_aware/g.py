import matplotlib.pyplot as plt
import networkx as nx

G = nx.read_graphml("albinism.graphml")
nx.draw(G, with_labels=True)
plt.show()
