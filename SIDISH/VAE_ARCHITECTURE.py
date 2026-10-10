################################## IMPORT LIBRARY ############################
import os
import torch
import torch.nn as nn
import numpy as np

import pyro.contrib.examples.util
from pyro.distributions.zero_inflated import ZeroInflatedNegativeBinomial

import torch.nn.functional as F
from torch.distributions import kl_divergence as KL
from torch.distributions import Normal

from torch_geometric.nn import GATConv

pyro.distributions.enable_validation(False)

################################ BUILDING BLOCKS #########################################

class GATConvLayer(nn.Module):
    """
    Single-layer Graph Attention wrapper.

    Uses concat=False so output dim equals out_features (averaged across heads),
    keeping the rest of the network agnostic to the head count.
    """
    def __init__(self, in_features, out_features, heads=4, dropout=0.0):
        super(GATConvLayer, self).__init__()
        self.gat = GATConv(in_features, out_features, heads=heads,
                           concat=False, dropout=dropout, add_self_loops=True)

    def forward(self, x, edge_index, edge_weight=None):
        # GATConv accepts edge_attr for attention bias; pass distances if available.
        return self.gat(x, edge_index, edge_attr=edge_weight)


class MLPBlock(nn.Module):
    """
    Linear -> (Dropout) -> Softplus stack, identical to the original encoder/decoder MLP.
    Factored out so DualEncoder and Decoder share the same building block.
    """
    def __init__(self, input_dim, layer_dims, dropout=0.0):
        super(MLPBlock, self).__init__()
        layers = []
        in_dim = input_dim
        for dim in layer_dims:
            layers.append(nn.Linear(in_dim, dim))
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            layers.append(nn.Softplus())
            in_dim = dim
        self.net = nn.Sequential(*layers)
        self.out_dim = in_dim

    def forward(self, x):
        return self.net(x)


################################ DECODERS #########################################

class Decoder(nn.Module):
    """
    Expression decoder (DecExp). Identical behaviour to the original SIDISH Decoder:
    z -> MLP -> two heads (mu, dropout_logits) for the ZINB distribution.
    """
    def __init__(self, input_dim, z_dim, layer_dims, dropout=0):
        super(Decoder, self).__init__()
        self.before_last_layer = MLPBlock(z_dim, layer_dims, dropout)
        in_dim = self.before_last_layer.out_dim
        self.last_layer_1 = nn.Linear(in_dim, input_dim)
        self.last_layer_2 = nn.Linear(in_dim, input_dim)

    def forward(self, z):
        bll = self.before_last_layer(z)
        mu = self.last_layer_1(bll)
        dropout_logits = self.last_layer_2(bll)
        return torch.exp(mu), dropout_logits


class DecEnv(nn.Module):
    """
    Environment decoder (COVET reconstruction head).

    z_niche -> MLP -> chol_size raw values -> lower-triangular L (softplus on diag)
            -> L @ L.T (guaranteed PSD)

    Output matches the pre-computed sqrt(COVET) target [n_covet_genes, n_covet_genes]
    under a Frobenius loss.
    """
    def __init__(self, z_dim, n_covet_genes, hidden_dims=(128, 256)):
        super(DecEnv, self).__init__()
        self.n_covet_genes = n_covet_genes
        self.chol_size = n_covet_genes * (n_covet_genes + 1) // 2
        self.mlp = MLPBlock(z_dim, list(hidden_dims), dropout=0.0)
        self.head = nn.Linear(self.mlp.out_dim, self.chol_size)

        tril_idx = torch.tril_indices(n_covet_genes, n_covet_genes)
        self.register_buffer('tril_row', tril_idx[0])
        self.register_buffer('tril_col', tril_idx[1])

        diag_mask = (tril_idx[0] == tril_idx[1])
        self.register_buffer('diag_mask', diag_mask)

    def forward(self, z_niche):
        B = z_niche.size(0)
        n = self.n_covet_genes

        raw = self.head(self.mlp(z_niche))                                    # [B, chol_size]
        diag = F.softplus(raw[:, self.diag_mask]) + 1e-6
        offdiag = raw[:, ~self.diag_mask]

        L = torch.zeros(B, n, n, device=z_niche.device, dtype=z_niche.dtype)
        L[:, self.tril_row[self.diag_mask], self.tril_col[self.diag_mask]] = diag
        L[:, self.tril_row[~self.diag_mask], self.tril_col[~self.diag_mask]] = offdiag

        return L @ L.transpose(-1, -2)                                        # [B, n, n], PSD


################################ ENCODER #########################################

class DualEncoder(nn.Module):
    """
    Shared MLP trunk producing h_hidden, with two heads:

      * Expression head (always active):
            h_hidden -> fc_expr_mean / fc_expr_logvar -> z_expr
        This is the only path used at bulk inference and feeds DEEP_COX.

      * Niche head (active only when edge_index is provided and gat_hidden is set):
            h_hidden -> GATConv(gat_hidden) -> fc_niche_mean / fc_niche_logvar -> z_niche
        This is what the COVET decoder consumes during training.

    Returns (mu_expr, lv_expr, mu_niche, lv_niche) with niche tensors None if
    the niche branch is inactive.
    """
    def __init__(self, input_dim, z_dim, layer_dims, dropout=0.0,
                 gat_hidden=None, n_heads=4, gat_dropout=0.0):
        super(DualEncoder, self).__init__()

        self.trunk = MLPBlock(input_dim, layer_dims, dropout)
        # Back-compat alias: DEEP_COX_ARCHITECTURE looks up `before_last_layer` on
        # the inner encoder to wire its transfer-learning MLP.
        self.before_last_layer = self.trunk
        h_dim = self.trunk.out_dim

        # Expression head (always present)
        self.fc_expr_mean = nn.Linear(h_dim, z_dim)
        self.fc_expr_logvar = nn.Linear(h_dim, z_dim)

        # Niche head (optional)
        self.use_niche = gat_hidden is not None
        if self.use_niche:
            self.gat = GATConvLayer(h_dim, gat_hidden, heads=n_heads, dropout=gat_dropout)
            self.fc_niche_mean = nn.Linear(gat_hidden, z_dim)
            self.fc_niche_logvar = nn.Linear(gat_hidden, z_dim)

    def forward(self, x, edge_index=None, edge_weight=None):
        h = self.trunk(x)

        mu_expr = self.fc_expr_mean(h)
        lv_expr = self.fc_expr_logvar(h)

        if self.use_niche and edge_index is not None:
            g = F.relu(self.gat(h, edge_index, edge_weight))
            mu_niche = self.fc_niche_mean(g)
            lv_niche = self.fc_niche_logvar(g)
        else:
            mu_niche = None
            lv_niche = None

        return mu_expr, lv_expr, mu_niche, lv_niche


################################ FULL MODEL #########################################

class ARCHITECTURE(nn.Module):
    """
    SIDISH-COVET VAE.

    Non-spatial (bulk) path:        x -> DualEncoder.trunk -> z_expr -> Decoder (ZINB)
    Spatial (training) path adds:    z_niche -> DecEnv (sqrt(COVET) reconstruction)

    n_covet_genes=None disables the COVET branch entirely (DecEnv is not built,
    GAT is not built); the model collapses to the original SIDISH non-spatial VAE.
    """
    def __init__(self, input_dim, z_dim, layer_dims, seed, dropout=0.5,
                 n_covet_genes=None, gat_hidden=None, n_heads=4, use_cuda=False):
        super(ARCHITECTURE, self).__init__()

        self.seed = seed
        torch.manual_seed(self.seed)
        torch.cuda.manual_seed(self.seed)
        torch.backends.cudnn.deterministic = True
        np.random.seed(self.seed)
        os.environ["PYTHONHASHSEED"] = str(self.seed)

        self.input_dim = input_dim
        self.z_dim = z_dim
        self.layer_dims = layer_dims
        self.use_cuda = use_cuda
        self.n_covet_genes = n_covet_genes
        self.use_spatial = (n_covet_genes is not None) and (gat_hidden is not None)

        self.log_theta = nn.Parameter(torch.randn(self.input_dim))

        self.encoder = DualEncoder(input_dim, z_dim, layer_dims, dropout,
                                   gat_hidden=gat_hidden if self.use_spatial else None,
                                   n_heads=n_heads)
        self.decoder = Decoder(input_dim, z_dim, layer_dims[::-1], dropout)

        if self.use_spatial:
            self.dec_env = DecEnv(z_dim, n_covet_genes)

        # Kept for back-compat with VAE.getBaseEncoder() and old checkpoints.
        self.before_last_layer = self.encoder.trunk

    # ---------- VAE helpers ----------

    def reparameterize(self, mu, logvar):
        return Normal(mu, logvar).rsample()

    def forward(self, x, edge_index=None, edge_weight=None):
        x = torch.log(x + 1)
        mu_expr, lv_expr, mu_niche, lv_niche = self.encoder(
            x.view(-1, self.input_dim), edge_index, edge_weight
        )
        z_expr = self.reparameterize(mu_expr, lv_expr)
        mu_decoder, dropout_logits = self.decoder(z_expr)

        if self.use_spatial and mu_niche is not None:
            z_niche = self.reparameterize(mu_niche, lv_niche)
            covet_hat = self.dec_env(z_niche)
        else:
            covet_hat = None

        return mu_decoder, dropout_logits, mu_expr, lv_expr, mu_niche, lv_niche, covet_hat

    def get_latent_representation(self, x, edge_index=None, edge_weight=None):
        """
        Returns z_expr (the expression latent) regardless of spatial state.
        This is the latent consumed by DEEP_COX, so it MUST be the bulk-compatible one.
        """
        x = torch.log(x + 1)
        mu_expr, lv_expr, _, _ = self.encoder(
            x.view(-1, self.input_dim), edge_index, edge_weight
        )
        return mu_expr + torch.exp(0.5 * lv_expr)

    # ---------- Loss components ----------

    def kl_d(self, mu, logvar):
        z_loc = torch.zeros_like(mu)
        z_scale = torch.ones_like(logvar)
        return KL(Normal(mu, logvar), Normal(z_loc, z_scale)).sum(dim=1)

    def reconstruction_loss(self, x, mu, dropout_logits, w):
        theta = F.softplus(self.log_theta)
        nb_logits = (mu + 1e-5).log() - (theta + 1e-5).log()
        distribution = ZeroInflatedNegativeBinomial(
            total_count=theta, logits=nb_logits,
            gate_logits=dropout_logits, validate_args=False
        )
        log_prob = distribution.log_prob(x) * w
        return -log_prob.sum(dim=-1)

    def env_loss(self, covet_hat, covet_target):
        """
        Frobenius^2 between predicted and target sqrt(COVET), averaged over the
        n*n matrix entries so its scale is comparable to per-gene losses.
        """
        diff = covet_hat - covet_target
        return (diff * diff).flatten(1).mean(dim=1)

    def loss_function(self, x, w, mu_decoder, dropout_logits, mu_expr, lv_expr):
        """Non-spatial loss — identical to original SIDISH."""
        rec = self.reconstruction_loss(x, mu_decoder, dropout_logits, w)
        kl = self.kl_d(mu_expr, lv_expr)
        return torch.mean(rec, dim=0) + torch.mean(kl, dim=0)

    def loss_function_spatial(self, x, w, covet_target,
                              mu_decoder, dropout_logits,
                              mu_expr, lv_expr, mu_niche, lv_niche, covet_hat,
                              lambda_env=1.0, beta_expr=1.0, beta_niche=0.5):
        """
        L = L_exp + lambda * L_env + beta_expr * KL(z_expr) + beta_niche * KL(z_niche)
        """
        rec = self.reconstruction_loss(x, mu_decoder, dropout_logits, w)
        kl_expr = self.kl_d(mu_expr, lv_expr)
        kl_niche = self.kl_d(mu_niche, lv_niche)
        env = self.env_loss(covet_hat, covet_target)

        return (
            torch.mean(rec, dim=0)
            + lambda_env * torch.mean(env, dim=0)
            + beta_expr * torch.mean(kl_expr, dim=0)
            + beta_niche * torch.mean(kl_niche, dim=0)
        )

######################################## END ############################################
