from SIDISH.VAE_ARCHITECTURE import ARCHITECTURE as architecture
import numpy as np
import torch
from torch.utils.data import TensorDataset, DataLoader
import pyro
import torch.optim as optim
from torch_geometric.loader import NeighborLoader
from torch_geometric.data import Data
import os

class VAE():
    def __init__(self, epochs, adata, z_dim, layer_dims, lr, dropout, device, seed,
                 n_covet_genes=None, gat_hidden=None, n_heads=4):
        super(VAE, self).__init__()

        self.seed = seed

        torch.manual_seed(self.seed)
        torch.cuda.manual_seed(self.seed)
        torch.backends.cudnn.deterministic = True
        np.random.seed(self.seed)
        os.environ["PYTHONHASHSEED"] = str(self.seed)

        ## Model parameters
        self.epochs = epochs
        self.imput_dim = adata.X.shape[1]
        self.lr = lr
        self.z_dim = z_dim
        self.layer_dims = layer_dims
        self.dropout = dropout
        self.device = device

        ## Spatial / COVET parameters
        self.n_covet_genes = n_covet_genes
        self.gat_hidden = gat_hidden
        self.n_heads = n_heads
        self.use_spatial = (n_covet_genes is not None) and (gat_hidden is not None)

    def initialize(self, adata, W=None, batch_size=1024, type="Normal", num_workers=8,
                   spatial_graph=None, num_neighbors=5, covet_sqrt=None):

        ## Initialise model
        self.adata = adata

        model = architecture(
            self.imput_dim, self.z_dim, self.layer_dims, self.seed, self.dropout,
            n_covet_genes=self.n_covet_genes if self.use_spatial else None,
            gat_hidden=self.gat_hidden if self.use_spatial else None,
            n_heads=self.n_heads,
        )

        self.model = model.to(self.device)
        self.optimizer = optim.Adam(lr=self.lr, params=self.model.parameters())

        ## Initialise weights
        if W is None:
            W = np.ones(self.adata.X.shape)
        self.W = W

        ## Get the cells by genes matrix X from the adata variable
        if type == 'Dense':
            X_arr = np.array(self.adata.X.todense())
        else:
            X_arr = np.array(self.adata.X)

        ## Prep single cell data for training in VAE -- cell by genes matrix X with gene weight matrix W
        self.X_tensor = torch.from_numpy(np.array(X_arr)).float()
        self.W_tensor = torch.from_numpy(np.array(self.W)).float()

        dataset = TensorDataset(self.X_tensor, self.W_tensor)
        kwargs = {'num_workers': num_workers, 'pin_memory': True}
        self.train_loader = DataLoader(dataset, batch_size=batch_size, **kwargs, drop_last=True)

        ## Setup mini-batching for the full dataset (prediction)
        self.total_loader = DataLoader(dataset, batch_size=self.adata.X.shape[0], **kwargs)

        ## Setup spatial graph and neighbor sampling
        self.use_neighbor_sampling = False
        self.covet_tensor = None

        if self.use_spatial and spatial_graph is not None:
            self.edge_index, self.edge_weight = spatial_graph

            if self.edge_weight is None:
                self.edge_weight = torch.ones(self.edge_index.size(1))

            # Attach COVET targets to the Data object so NeighborLoader subsets them
            # per-batch alongside x / y.
            data_kwargs = dict(
                x=self.X_tensor,
                edge_index=self.edge_index,
                edge_attr=self.edge_weight,
                y=self.W_tensor,
            )
            if covet_sqrt is not None:
                self.covet_tensor = torch.from_numpy(np.array(covet_sqrt)).float()
                data_kwargs['covet'] = self.covet_tensor

            self.graph_data = Data(**data_kwargs)

            self.use_neighbor_sampling = True
            self.train_loader_spatial = NeighborLoader(
                self.graph_data, num_neighbors=[num_neighbors],
                batch_size=batch_size, shuffle=True, num_workers=num_workers,
            )

        if self.use_spatial and self.use_neighbor_sampling:
            return self.model, self.train_loader_spatial
        else:
            return self.model, self.train_loader


    def train(self, lambda_env=1.0, beta_expr=1.0, beta_niche=0.5):
        """
        Training loop.

        - Non-spatial (bulk / no COVET): single-decoder ZINB + KL, identical to original SIDISH.
        - Spatial: dual-decoder loss
              L = L_exp + lambda_env * L_env + beta_expr * KL(z_expr) + beta_niche * KL(z_niche).
        """
        self.loss = []

        for epoch in range(self.epochs):
            epoch_loss = 0.
            samples_processed = 0

            if self.use_spatial and self.use_neighbor_sampling:
                for batch in self.train_loader_spatial:
                    batch = batch.to(self.device)
                    batch_size = batch.x.size(0)

                    self.optimizer.zero_grad()

                    mu_dec, drop_logits, mu_e, lv_e, mu_n, lv_n, covet_hat = self.model(
                        batch.x, batch.edge_index, batch.edge_attr
                    )

                    if self.covet_tensor is not None and hasattr(batch, 'covet'):
                        loss = self.model.loss_function_spatial(
                            batch.x, batch.y, batch.covet,
                            mu_dec, drop_logits,
                            mu_e, lv_e, mu_n, lv_n, covet_hat,
                            lambda_env=lambda_env,
                            beta_expr=beta_expr,
                            beta_niche=beta_niche,
                        )
                    else:
                        # Spatial encoder but no COVET targets -> fall back to expression-only loss.
                        loss = self.model.loss_function(
                            batch.x, batch.y, mu_dec, drop_logits, mu_e, lv_e
                        )

                    loss.backward()
                    epoch_loss += loss.item() * batch_size
                    samples_processed += batch_size
                    self.optimizer.step()

            else:
                for x, y in self.train_loader:
                    x = x.to(self.device, non_blocking=True)
                    y_ = y.to(self.device, non_blocking=True)
                    batch_size = x.size(0)
                    self.optimizer.zero_grad()
                    mu_dec, drop_logits, mu_e, lv_e, _, _, _ = self.model(x)
                    loss = self.model.loss_function(x, y_, mu_dec, drop_logits, mu_e, lv_e)
                    loss.backward()
                    epoch_loss += loss.item() * batch_size
                    samples_processed += batch_size
                    self.optimizer.step()

            total_epoch_loss_train = epoch_loss / samples_processed
            self.epoch_loss = total_epoch_loss_train
            self.loss.append(-total_epoch_loss_train)
            print("[epoch %03d]  average training loss: %.4f" % (epoch, self.epoch_loss))

    def getBaseEncoder(self):
        return self.model.encoder.trunk

    def getLoss(self):
        self.error = np.array([i for i in self.loss])
        return self.error

    def getEmbedding(self, clustering=True):
        """
        Get z_expr (the bulk-compatible latent) for every cell.

        z_expr depends only on the MLP trunk and the expression head, so this is
        identical in spatial and non-spatial training modes and needs no graph.
        """
        self.model.eval()
        self.TZ = []
        with torch.no_grad():
            for x, y in self.total_loader:
                pyro.clear_param_store()
                x = x.to(self.device, non_blocking=True)
                z = self.model.get_latent_representation(x)
                self.TZ += z.cpu().detach().numpy().tolist()

        if clustering:
            self.adata.obsm['latent'] = np.array(self.TZ).astype(np.float32)
        return self.adata

    def getBaseEmbedding(self, clustering=True):
        """
        Alias of getEmbedding for backwards compatibility. In the new dual-encoder,
        the "base" (Cox-compatible) latent IS z_expr.
        """
        self.base_TZ = []
        self.model.eval()
        with torch.no_grad():
            for x, y in self.total_loader:
                pyro.clear_param_store()
                x = x.to(self.device, non_blocking=True)
                z = self.model.get_latent_representation(x)
                self.base_TZ += z.cpu().detach().numpy().tolist()

        if clustering:
            self.adata.obsm['base_latent'] = np.array(self.base_TZ).astype(np.float32)
        return self.adata
