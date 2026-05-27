"""Batch-adjusted linear models with limma-style empirical Bayes moderation.

This is an independent implementation of a subset of limma's ordinary eBayes
model, not a call to Bioconductor limma. Equal residual degrees of freedom and
complete data are required. It should be cross-validated in R before clinical use.
"""
import numpy as np
import pandas as pd
from scipy import optimize, special, stats
from statsmodels.stats.multitest import multipletests

def fit_moderated(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    if list(expression.columns) != list(metadata.sample_id):
        raise ValueError("Expression / metadata misalignment")
    if set(metadata.subtype) != {"Luminal A", "Luminal B"}:
        raise ValueError("Requires both luminal groups")
    if expression.isna().any().any():
        raise ValueError("Missing expression")
    batch = pd.get_dummies(metadata.batch.astype(str), drop_first=True, dtype=float)
    design = np.column_stack([np.ones(len(metadata)),
                              (metadata.subtype == "Luminal B").astype(float), batch.to_numpy()])
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Confounded design")
    y = expression.to_numpy(dtype=float).T
    xtx_inv = np.linalg.inv(design.T@design)
    beta = xtx_inv@design.T@y
    residual = y-design@beta
    df = len(metadata)-design.shape[1]
    s2 = np.maximum(np.sum(residual**2, axis=0)/df, np.finfo(float).tiny)
    log_s2 = np.log(s2)
    excess = max(np.var(log_s2, ddof=1)-special.polygamma(1, df/2), 1e-8)
    df0 = 2*optimize.brentq(lambda x: special.polygamma(1, x)-excess, 1e-5, 1e9)
    log_scale = (log_s2.mean()-special.digamma(df/2)+np.log(df/2)
                 +special.digamma(df0/2)-np.log(df0/2))
    scale = np.exp(log_scale)
    posterior = (df0*scale+df*s2)/(df0+df)
    se = np.sqrt(posterior*xtx_inv[1, 1])
    t = beta[1]/se
    p = 2*stats.t.sf(np.abs(t), df+df0)
    q = multipletests(p, method="fdr_bh")[1]
    result = pd.DataFrame({"gene_symbol": expression.index, "log2_fold_change": beta[1],
                           "effect_size": beta[1], "t_moderated": t, "p_value": p,
                           "adjusted_p_value": q, "significant_fdr05": q < 0.05,
                           "mean_log2_intensity": y.mean(axis=0)})
    return result.sort_values("adjusted_p_value").reset_index(drop=True)
