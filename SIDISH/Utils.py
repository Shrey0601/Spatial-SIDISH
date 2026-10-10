import torch
import numpy as np
from scipy.stats import weibull_min
from imblearn.over_sampling import RandomOverSampler
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from tqdm import tqdm
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.utils.class_weight import compute_class_weight
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from imblearn.under_sampling import RandomUnderSampler
import copy
import shap
import scanpy as sc
import pandas as pd
from scipy import stats
from scipy.linalg import sqrtm
from typing import Literal
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import kneighbors_graph, NearestNeighbors
import torch_geometric.utils as pyg_utils
from torch_geometric.data import Data
from scipy.stats import binomtest
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

def extractFeature(adata,type='Normal'):
    feature_acc = []

    C = adata.to_df().values
    risk = np.array(adata.obs.SIDISH_value.values)

    train_X, test_X, train_y, test_y = train_test_split(C, risk, test_size=0.2, stratify=risk, random_state=42)
    ss_train = StandardScaler()
    train_X = ss_train.fit_transform(train_X)
    test_X = ss_train.transform(test_X)

    classes = np.unique(train_y)
    class_weights = compute_class_weight(class_weight="balanced", classes=np.array([0, 1]), y=train_y)
    class_weight_dict = dict(zip(classes, class_weights))

    ros = RandomUnderSampler(random_state=42)
    train_X, train_y = ros.fit_resample(train_X, train_y)

    #model = LogisticRegression(penalty="l1", solver="liblinear", max_iter=10000,class_weight=class_weight_dict, n_jobs=-1, random_state=42).fit(train_X, train_y)
    model = RandomForestClassifier(criterion='log_loss',n_jobs=-1, class_weight=class_weight_dict,random_state=42).fit(train_X, train_y)
    return model, train_X, test_X


def sigmoid(x, a=100, b=0):
    value = 1 / (1 + np.exp(-a * (x - b)))
    return value

def create_spatial_graph(spatial_coords, method='knn', k=5, radius=None, include_self=False):
    """
    Create a graph from spatial coordinates with limited neighborhood size
    
    Parameters:
    -----------
    spatial_coords : numpy.ndarray
        Array of shape [n_cells, n_dims] containing spatial coordinates
    method : str
        Method to create the graph, either 'knn' or 'radius'
    k : int
        Number of nearest neighbors for KNN graph (default: 5)
    radius : float
        Radius for radius graph
    include_self : bool
        Whether to include self-loops
    
    Returns:
    --------
    edge_index : torch.LongTensor
        Edge indices in COO format [2, num_edges]
    edge_weight : torch.FloatTensor
        Edge weights based on distance [num_edges]
    """
    n_cells = spatial_coords.shape[0]
    
    if method == 'knn':
        # Create KNN graph with limited neighborhood size
        adjacency = kneighbors_graph(
            spatial_coords, k, mode='distance', include_self=include_self
        )
        
    elif method == 'radius':
        if radius is None:
            raise ValueError("Radius must be provided for radius graph method")
        
        # Calculate pairwise distances
        dist_matrix = squareform(pdist(spatial_coords))
        
        # Create binary adjacency matrix based on radius
        adjacency = (dist_matrix <= radius).astype(np.float32)
        
        # Remove self-loops if needed
        if not include_self:
            np.fill_diagonal(adjacency, 0)
    else:
        raise ValueError(f"Unknown graph creation method: {method}")
    
    # Convert to edge_index and edge_weight
    adjacency_coo = adjacency.tocoo()
    edge_index = torch.LongTensor(np.vstack((adjacency_coo.row, adjacency_coo.col)))
    edge_weight = torch.FloatTensor(adjacency_coo.data)
    
    # Convert weights to similarity (1 / (distance + epsilon))
    edge_weight = 1.0 / (edge_weight + 1e-6)
    
    return edge_index, edge_weight

def get_spatial_graph_from_adata(adata, spatial_key='spatial', method='knn', k=5, radius=None, include_self=False):
    """
    Create a graph from spatial coordinates stored in AnnData object
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing spatial coordinates
    spatial_key : str
        Key in adata.obsm where spatial coordinates are stored
    method : str
        Method to create the graph, either 'knn' or 'radius'
    k : int
        Number of nearest neighbors for KNN graph (default: 5)
    radius : float
        Radius for radius graph
    include_self : bool
        Whether to include self-loops
        
    Returns:
    --------
    edge_index : torch.LongTensor
        Edge indices in COO format [2, num_edges]
    edge_weight : torch.FloatTensor
        Edge weights based on distance [num_edges]
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")
    
    spatial_coords = adata.obsm[spatial_key]
    return create_spatial_graph(spatial_coords, method, k, radius, include_self)

def create_pytorch_geometric_data(adata, spatial_key='spatial', method='knn', k=5):
    """
    Create a PyTorch Geometric Data object from AnnData with spatial information
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing gene expression and spatial data
    spatial_key : str
        Key in adata.obsm where spatial coordinates are stored
    method : str
        Method to create the graph, either 'knn' or 'radius'
    k : int
        Number of nearest neighbors for KNN graph (default: 5)
        
    Returns:
    --------
    data : torch_geometric.data.Data
        PyTorch Geometric Data object containing the graph
    """
    # Extract expression data
    if isinstance(adata.X, np.ndarray):
        x = torch.from_numpy(adata.X).float()
    else:
        x = torch.from_numpy(adata.X.todense()).float()
    
    # Create spatial graph
    edge_index, edge_attr = get_spatial_graph_from_adata(
        adata, 
        spatial_key=spatial_key, 
        method=method, 
        k=k
    )
    
    # Create Data object
    data = Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr
    )
    
    return data


def compute_covet(adata, n_hvg=50, k=5, spatial_key='spatial', eps=1e-6, X_precomputed=None):
    """
    Compute per-cell shifted-covariance niche representation (COVET) and its
    matrix square root, following Haviv et al. (Nature Biotechnology, 2025).

    For each cell i with k spatial nearest neighbors:
        COVET_i = (1/k) * (E_i - X_bar)^T (E_i - X_bar)
    where E_i is the [k, n_hvg] expression of neighbors and X_bar is the
    GLOBAL mean over all cells (not the local niche mean). The global-mean
    shift is what makes niches directly comparable across the dataset.

    The matrix square root is returned because AOT distance between niches
    is L2^2 between their square roots, making it the natural decoder target
    for an MSE/Frobenius loss.

    Parameters
    ----------
    adata : AnnData
        AnnData with .X (cells x genes) and .obsm[spatial_key] (cells x 2/3).
    n_hvg : int
        Number of highly variable genes to use for COVET. Ignored if
        X_precomputed is given (its column count sets the dimension instead).
    k : int
        Number of spatial nearest neighbors per cell (self excluded).
    spatial_key : str
        Key in adata.obsm holding spatial coordinates.
    eps : float
        Ridge added to the diagonal of COVET before sqrtm for numerical
        stability when matrices are near-singular.
    X_precomputed : np.ndarray, optional
        [N, D] matrix to use directly as the per-cell niche features instead
        of the internal dispersion-based HVG selection (e.g. PCA components
        computed over the full gene panel). When given, the internal gene
        selection is skipped entirely and hvg_idx is returned as None.

    Returns
    -------
    covet_sqrt : np.ndarray
        Array of shape [N, D, D] containing sqrt(COVET_i) for each cell,
        where D = n_hvg or X_precomputed.shape[1].
    hvg_idx : np.ndarray or None
        Indices into adata.var of the selected HVGs (length n_hvg), or None
        when X_precomputed was used.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")

    if X_precomputed is not None:
        X_hvg = np.asarray(X_precomputed, dtype=np.float32)
        hvg_idx = None
    else:
        X = adata.X
        if not isinstance(X, np.ndarray):
            X = X.toarray()
        X = X.astype(np.float32)

        # Highly variable gene selection by dispersion (Fano factor: var / mean).
        # Works on raw counts, no extra dependencies. If 'highly_variable' is already
        # set in adata.var, defer to that selection instead.
        if 'highly_variable' in adata.var.columns:
            hvg_idx = np.where(adata.var['highly_variable'].values)[0][:n_hvg]
            if len(hvg_idx) < n_hvg:
                mean = X.mean(axis=0)
                var = X.var(axis=0)
                dispersion = var / (mean + 1e-8)
                hvg_idx = np.argsort(-dispersion)[:n_hvg]
        else:
            mean = X.mean(axis=0)
            var = X.var(axis=0)
            dispersion = var / (mean + 1e-8)
            hvg_idx = np.argsort(-dispersion)[:n_hvg]

        X_hvg = X[:, hvg_idx]

    n_hvg = X_hvg.shape[1]
    X_bar = X_hvg.mean(axis=0, keepdims=True)

    coords = adata.obsm[spatial_key]
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, knn_idx = nn.kneighbors(coords)
    knn_idx = knn_idx[:, 1:]

    N = X_hvg.shape[0]
    covet_sqrt = np.zeros((N, n_hvg, n_hvg), dtype=np.float32)
    ridge = eps * np.eye(n_hvg, dtype=np.float32)
    for i in tqdm(range(N), desc='Computing COVET'):
        E_i = X_hvg[knn_idx[i]]
        D = E_i - X_bar
        C = (D.T @ D) / float(k)
        S = sqrtm(C + ridge)
        covet_sqrt[i] = np.real(S).astype(np.float32)

    return covet_sqrt, hvg_idx


def r_squared(sample, dist, params):
    x_sorted  = np.sort(sample)
    n         = len(sample)
    ecdf      = np.arange(1, n + 1) / n
    model_cdf = dist.cdf(x_sorted, *params)
    sse = np.sum((ecdf - model_cdf) ** 2)
    sst = np.sum((ecdf - ecdf.mean()) ** 2)
    return 1 - sse / sst


class Utils():

    def getWeightVector(patients, adata, model, percentile: float, device: str, distribution: Literal["default", "fitted"] = "default", dist = None):

        labels = []
        p = patients
        all_p = p.to(device, non_blocking=True)
        val = model(all_p).detach().cpu().flatten()

        X = adata.to_df().values

        # X = np.log(1+X) #MinMaxScaler().fit_transform(np.log(1 + X))
        all_X = torch.from_numpy(X).type(torch.float32)
        all_X = all_X.to(device, non_blocking=True)
        val_sc = model(all_X).detach().cpu().flatten()

        ##### Cell high risk labeling #####
        cell_hazard = val_sc
        patient_hazard = val
        
        DISTROS = { "Weibull": stats.weibull_min, "Exponential": stats.expon, "Gamma" : stats.gamma }
            
        returned_dist =None
        if distribution == "default":
            params = weibull_min.fit(cell_hazard)
            params_patients = weibull_min.fit(patient_hazard)

            percentile_cells = weibull_min.ppf(percentile, *params)
            percentile_patients = weibull_min.ppf(percentile, *params_patients)
        
        else:
            if dist is not None:
                params = DISTROS[dist].fit(cell_hazard)
                params_patients = DISTROS[dist].fit(patient_hazard)

                percentile_cells = DISTROS[dist].ppf(percentile, *params)
                percentile_patients = DISTROS[dist].ppf(percentile, *params_patients)
                
            else:
                
                exp_cell_hazard = np.exp(cell_hazard)
                
                results = []
                for idx, (label, dist) in enumerate(DISTROS.items()):
                    # 4.1  Fit parameters
                    if label != "Weibull":
                        params = dist.fit(exp_cell_hazard, floc=0)
                    else:
                        params = dist.fit(exp_cell_hazard)
                    
                    # 4.2  Metrics
                    logL = np.sum(dist.logpdf(exp_cell_hazard, *params))
                    k, n = len(params), len(exp_cell_hazard)
                    aic  = 2 * k - 2 * logL
                    bic  = k * np.log(n) - 2 * logL
                    ks_D, ks_p = stats.kstest(exp_cell_hazard, dist.name, args=params)
                    r2 = r_squared(exp_cell_hazard, dist, params)

                    results.append(dict(distribution=label, logL=logL, AIC=aic, BIC=bic, KS_D=ks_D, KS_p=ks_p, R2=r2, params=params))
                df = pd.DataFrame(results).round(6)
                    
                # Hybrid rank-sum (lower is better)
                ranked = df.assign(
                    r_aic = df["AIC"].rank(method="min", ascending=True),
                    r_bic = df["BIC"].rank(method="min", ascending=True),
                    r_ks  = df["KS_D"].rank(method="min", ascending=True),
                    r_r2  = df["R2"].rank(method="min", ascending=False),
                )
                ranked["rank_sum"] = ranked[["r_aic","r_bic","r_ks","r_r2"]].sum(axis=1)

                # Pick best by rank_sum, tie-break by AIC
                best_idx = ranked.sort_values(["rank_sum","AIC"], ascending=[True, True]).index[0]
                best_hybrid = df.loc[best_idx, "distribution"]
                returned_dist = best_hybrid
                print(best_hybrid)
                print("Best Distribution: ".format(best_hybrid))
                    
                params = DISTROS[best_hybrid].fit(cell_hazard)
                params_patients = DISTROS[best_hybrid].fit(patient_hazard)
                    
                percentile_cells = DISTROS[best_hybrid].ppf(percentile, *params)
                percentile_patients = DISTROS[best_hybrid].ppf(percentile, *params_patients)

        high_risk_cells = cell_hazard >= percentile_cells
        high_risk_cells_ = high_risk_cells.type(torch.int)
        mask_patients = patient_hazard >= percentile_patients

        adata.obs["SIDISH_value"] = high_risk_cells_
        adata.obs["risk_value"] = val_sc

        for i in adata.obs.SIDISH_value:
            if i == 1:
                labels.append("h")
            else:
                labels.append("b")

        adata.obs["SIDISH"] = labels

        val_p = torch.sigmoid(val)
        val_p[~mask_patients] = 0.0

        return (torch.FloatTensor(val_p), adata, percentile_cells, percentile_cells.max(), percentile_cells.min(), returned_dist)
    
    def annotateCells(adata, model,percentile_cells, device, percentile, mode, perturbation=False):

        labels = []
        # X = np.log(1 + X)
        X = adata.to_df().values
        all_X = torch.from_numpy(X).type(torch.float32)
        all_X = all_X.to(device, non_blocking=True)
        val_sc = model(all_X).detach().cpu().flatten()

        ##### Cell high risk labeling #####
        cell_hazard = val_sc
        if mode == 'test':
            params = weibull_min.fit(cell_hazard)
            percentile_cells = weibull_min.ppf(percentile, *params)
        else:
            percentile_cells = percentile_cells

        high_risk_cells = cell_hazard >= percentile_cells
        high_risk_cells_ = high_risk_cells.type(torch.int)
        
        if perturbation:
            old_risk_score = adata.obs["risk_value"].values
            
            delta_score = old_risk_score - val_sc.numpy()
            adata.obs['perturbation_score'] = delta_score
            adata.obs['perturbation_score_abs'] = np.abs(delta_score)
            
            adata.obs["SIDISH_value"] = high_risk_cells_
            adata.obs["risk_value"] = val_sc
        
        else: 
            adata.obs["SIDISH_value"] = high_risk_cells_
            adata.obs["risk_value"] = val_sc

        for i in adata.obs.SIDISH_value:
            if i == 1:
                labels.append("h")
            else:
                labels.append("b")

        adata.obs["SIDISH"] = labels
        return  adata

    def getWeightMatrix(adata, seed, steepness=100,type='Normal'):
        model, train_X, test_X = extractFeature(adata, type)
        explainer = shap.Explainer(model, train_X, seed=seed, feature_names=adata.to_df().columns)
        shap_values = explainer(test_X, check_additivity=False)
        shap_values_class_1 = shap_values.values[:, :, 0]
        shap_values = pd.DataFrame(copy.deepcopy(shap_values_class_1), columns=adata.to_df().columns)
        q = shap_values.mean(axis=0)

        W = []
        for i in adata.obs.SIDISH:
            if i in "h":
                x = q.values
                weights_ = sigmoid(x, steepness)
                W.append(weights_)

            else:
                W.append(np.zeros(adata.shape[1]))
        W = np.array(W)

        return W

    def get_threshold(adata, model, percentile, device):

        X = adata.to_df().values
        # X = np.log(1+X) #MinMaxScaler().fit_transform(np.log(1 + X))
        all_X = torch.from_numpy(X).type(torch.float32)
        all_X = all_X.to(device, non_blocking=True)
        val_sc = model(all_X).detach().cpu().flatten()

        ##### Cell high risk labeling #####
        cell_hazard = val_sc

        params = weibull_min.fit(cell_hazard)

        percentile_cells = weibull_min.ppf(percentile, *params)
        print(percentile_cells)

        return percentile_cells


def cox_weighted_external_validation(bulk_train, gene_list, expr_external, surv_external,
                                      penalizer_grid=(0.1, 0.5, 1, 2, 5, 10, 20),
                                      n_folds=5, seed=0):
    """
    Score an independent bulk cohort using Cox-weighted marker gene expression,
    instead of a naive unweighted mean.

    A ridge-penalized CoxPHFitter is fit on `bulk_train` (the SIDISH training
    cohort) to get signed, per-gene coefficients. The ridge penalizer is chosen
    by k-fold cross-validation on `bulk_train` ALONE (mean held-out concordance
    index), never touching `expr_external` — this avoids picking the penalizer
    against the validation set. The final model (refit on all of bulk_train
    with the CV-selected penalizer) then scores `expr_external` patients via
    risk = expression @ coef, splits at the median, and runs a log-rank test
    against `surv_external`.

    An unweighted mean can null out real signal when marker genes have
    opposing hazard directions (some raise risk, some lower it) -- averaging
    them cancels rather than reinforces. Cox weighting fixes this by signing
    and scaling each gene's contribution before combining.

    Parameters
    ----------
    bulk_train : pd.DataFrame
        Training cohort with 'duration', 'event' columns plus one column per
        gene in `gene_list`.
    gene_list : list of str
        Marker genes to use as Cox covariates. Only genes present in both
        `bulk_train` and `expr_external` are used.
    expr_external : pd.DataFrame
        External cohort expression, patients x genes.
    surv_external : pd.DataFrame
        External cohort survival, indexed to match `expr_external`, with
        'duration' and 'event' columns.
    penalizer_grid : tuple of float
        Candidate ridge penalizer strengths to cross-validate over.
    n_folds : int
        Number of CV folds for penalizer selection.
    seed : int
        Random seed for the fold split.

    Returns
    -------
    dict with keys: penalizer, coef (pd.Series), log_rank_p, high_group
    (pd.Series of bool, indexed like expr_external), risk_scores (np.ndarray)
    """
    avail = [g for g in gene_list if g in bulk_train.columns and g in expr_external.columns]
    if len(avail) < 2:
        raise ValueError(f"Need >=2 genes present in both cohorts, got {len(avail)}")

    cols = avail + ['duration', 'event']
    train_df = bulk_train[cols].copy()

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    mean_ci = {}
    for pen in penalizer_grid:
        fold_cis = []
        for train_idx, val_idx in skf.split(train_df, train_df['event']):
            fold_train = train_df.iloc[train_idx]
            fold_val = train_df.iloc[val_idx]
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(fold_train, duration_col='duration', event_col='event')
            fold_cis.append(cph.score(fold_val, scoring_method='concordance_index'))
        mean_ci[pen] = np.mean(fold_cis)

    best_penalizer = max(mean_ci, key=mean_ci.get)

    cph_final = CoxPHFitter(penalizer=best_penalizer)
    cph_final.fit(train_df, duration_col='duration', event_col='event')
    coef = cph_final.summary['coef'].loc[avail]

    risk_scores = expr_external[avail].values @ coef.values
    high_group = pd.Series(risk_scores >= np.median(risk_scores), index=expr_external.index)

    lr = logrank_test(
        surv_external.loc[high_group, 'duration'], surv_external.loc[~high_group, 'duration'],
        event_observed_A=surv_external.loc[high_group, 'event'],
        event_observed_B=surv_external.loc[~high_group, 'event'],
    )

    return dict(
        penalizer=best_penalizer,
        cv_concordance=mean_ci,
        coef=coef,
        log_rank_p=lr.p_value,
        high_group=high_group,
        risk_scores=risk_scores,
    )


def paper_matched_insample_validation(gene_list, expr_external, surv_external, penalizer=0.0):
    """
    DIAGNOSTIC ONLY -- NOT a test of generalization. Do not report this as external
    validation on its own; always pair it with cox_weighted_external_validation()'s
    blind (TCGA-fit) result and label both clearly.

    Matches the SIDISH paper's own described Kaplan-Meier methodology (Methods,
    "Kaplan-Meier survival analysis"): a Cox model is fit DIRECTLY on the same
    cohort it is then tested on -- i.e. `expr_external`/`surv_external` supplies
    both the coefficients AND the test set. This lets the model see the external
    cohort's own survival outcomes before scoring it, so a small p-value here
    only shows the gene list CAN be tuned to fit this cohort's own outcomes, not
    that it generalizes to unseen patients. Use it to check whether a poor blind
    result is due to weak marker genes (this would also look poor) or genuinely
    just to test-strictness (this would look strong even when the blind result
    doesn't).

    Parameters
    ----------
    gene_list : list of str
        Marker genes to use as Cox covariates. Only genes present in
        `expr_external` are used.
    expr_external : pd.DataFrame
        External cohort expression, patients x genes -- supplies BOTH the
        Cox-fitting data and the test set (paper's described approach).
    surv_external : pd.DataFrame
        External cohort survival, indexed to match `expr_external`, with
        'duration' and 'event' columns.
    penalizer : float
        Ridge penalizer for numerical stability only (the paper describes a
        "simple Cox regression model" with no mention of regularization;
        0.0 is lifelines' own default and matches that description literally
        -- raised only if the fit fails to converge).

    Returns
    -------
    dict with keys: penalizer, coef (pd.Series), log_rank_p, high_group,
    risk_scores (np.ndarray), converged (bool), note (str)
    """
    avail = [g for g in gene_list if g in expr_external.columns]
    if len(avail) < 2:
        raise ValueError(f"Need >=2 genes present in expr_external, got {len(avail)}")

    cols = avail + ['duration', 'event']
    df = expr_external.join(surv_external[['duration', 'event']])[cols].dropna()

    converged = True
    pen = penalizer
    cph = CoxPHFitter(penalizer=pen)
    try:
        cph.fit(df, duration_col='duration', event_col='event')
    except Exception:
        converged = False
        pen = max(penalizer, 0.1)
        cph = CoxPHFitter(penalizer=pen)
        cph.fit(df, duration_col='duration', event_col='event')

    coef = cph.summary['coef'].loc[avail]
    risk_scores = df[avail].values @ coef.values
    high_group = pd.Series(risk_scores >= np.median(risk_scores), index=df.index)

    lr = logrank_test(
        df.loc[high_group, 'duration'], df.loc[~high_group, 'duration'],
        event_observed_A=df.loc[high_group, 'event'],
        event_observed_B=df.loc[~high_group, 'event'],
    )

    return dict(
        penalizer=pen,
        coef=coef,
        log_rank_p=lr.p_value,
        high_group=high_group,
        risk_scores=risk_scores,
        converged=converged,
        note="DIAGNOSTIC ONLY -- Cox fit directly on the same cohort tested (paper-matched "
             "methodology). Not a generalization test; pair with the blind TCGA-fit result.",
    )