import torch


class SplitNet(torch.nn.Module):
    def __init__(self, num_nodes, num_embedding=64, device='gpu'):
        super().__init__()
        self.num_nodes = num_nodes
        self.device = device
        self.num_embedding = num_embedding
        self.w = torch.nn.Parameter(torch.randn([num_nodes, num_embedding]))
        self.zero_vector = torch.nn.Parameter(torch.randn([num_embedding]))
        self.lstm = torch.nn.LSTM(num_embedding, num_embedding, 2, batch_first=True)
        self.linear_1 = torch.nn.Linear(2 * num_embedding, num_embedding)
        self.activation = torch.nn.ReLU()
        self.linear_2 = torch.nn.Linear(num_embedding, 1)

    def forward(self, splits, maps):
        encoding_splits = splits @ self.w
        with torch.no_grad():
            embeddings = torch.nn.Embedding(len(splits)+1, self.num_embedding, padding_idx=0)
            embeddings.weight = torch.nn.Parameter(torch.cat([torch.unsqueeze(self.zero_vector, dim=0), encoding_splits]))
            embeddings.weight.requires_grad = False  # set to non-trainable
            maps = maps * (torch.arange(len(splits)+1)[1:].type(torch.long).unsqueeze(dim=-1)).to('cuda')
            x = embeddings(maps.type(torch.long))
        out, (state, cell) = self.lstm(x)
        state = state.permute((1,2,0)).reshape((-1, 2*self.num_embedding))
        avg_state = state.mean(axis=0)
        # print('avg_state', avg_state.shape)
        output = self.linear_1(avg_state)
        output = self.activation(output)
        output = self.linear_2(output)
        return output

