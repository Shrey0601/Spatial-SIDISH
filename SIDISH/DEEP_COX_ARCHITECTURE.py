import torch.nn as nn

class DEEPCOX_ARCHITECTURE(nn.Module):
    """
    Deep Cox architecture used in SIDISH for survival prediction.
    This network integrates pretrained encoder representations (from the VAE)
    with a Cox proportional hazards regression layer for modeling survival risk.
    """
    
    def __init__(self, hidden, encoder, dropout):
        """
        Args:
            hidden (int): Number of hidden units in the intermediate layer.
            encoder (nn.Module): Pretrained encoder (e.g., from SIDISH VAE).
            dropout (float): Dropout rate for regularization.
        """
        
        super(DEEPCOX_ARCHITECTURE, self).__init__()

        # Cox uses the shared MLP trunk of the VAE encoder. In the dual-encoder
        # design (post-COVET), `before_last_layer` is the trunk in both spatial
        # and non-spatial modes; its output feeds either the expression head
        # (always) or, in spatial mode, the GAT/niche head as well.
        enc_module = encoder.model.encoder
        self.encoder_layer = enc_module.before_last_layer
        if hasattr(self.encoder_layer, 'out_dim'):
            last_out = self.encoder_layer.out_dim
        else:
            # Fall back to the expression-mean layer's in_features, which equals
            # the trunk's output dim.
            fc_mean = getattr(enc_module, 'fc_expr_mean',
                              getattr(enc_module, 'fc_mean', None))
            last_out = fc_mean.in_features
        self.af1 = nn.Tanh()
        self.dr1 = nn.Dropout(dropout)

        # Add a new fully connected layer for risk prediction
        self.new_layer = nn.Linear(last_out, hidden)
        self.dr2 = nn.Dropout(dropout)
        self.af2 = nn.Tanh()
        self.final_layer = nn.Linear(hidden, 1, bias=False)

    def forward(self, x):
        x_ = self.af1(self.dr1(self.encoder_layer(x)))
        x__ = self.af2(self.dr2(self.new_layer(x_)))
        final_x = self.final_layer(x__)
        return final_x
