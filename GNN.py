import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from torch_geometric.utils import from_networkx


class GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super(GCN, self).__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index

        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)

        return F.log_softmax(x, dim=1)


def create_pyg_data(G, seed_nodes):
    # Add dummy node features if none exist
    for node in G.nodes():
        G.nodes[node]['x'] = [1.0]  # simple 1-feature vector

    # Add labels: 1 for seed nodes, 0 otherwise
    for node in G.nodes():
        G.nodes[node]['y'] = 1 if node in seed_nodes else 0

    data = from_networkx(G)
    data.x = torch.tensor([G.nodes[n]['x'] for n in G.nodes()], dtype=torch.float)
    data.y = torch.tensor([G.nodes[n]['y'] for n in G.nodes()], dtype=torch.long)

    return data


def train(model, data, epochs=200):
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

    model.train()
    for epoch in range(epochs):
        optimizer.zero_grad()
        out = model(data)
        loss = F.nll_loss(out, data.y)
        loss.backward()
        optimizer.step()

        if epoch % 20 == 0:
            pred = out.argmax(dim=1)
            acc = (pred == data.y).sum().item() / data.num_nodes
            print(f'Epoch {epoch} | Loss: {loss.item():.4f} | Accuracy: {acc:.4f}')
