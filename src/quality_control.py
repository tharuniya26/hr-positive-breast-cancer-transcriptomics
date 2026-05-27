"""Distribution and PCA diagnostics without modifying source values."""
import numpy as np
import pandas as pd

def qc_metrics(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    vals = expression.to_numpy(dtype=float)
    if vals.shape[1] != len(metadata):
        raise ValueError("Expression / metadata mismatch")
    if np.isnan(vals).any():
        raise ValueError("Missing expression")
    out = metadata.copy()
    out["median_intensity"] = np.median(vals, axis=0)
    out["iqr_intensity"] = np.percentile(vals, 75, axis=0)-np.percentile(vals, 25, axis=0)
    out["missing_fraction"] = np.mean(~np.isfinite(vals), axis=0)
    return out

def pca(expression: pd.DataFrame, metadata: pd.DataFrame, top_variable: int = 3000):
    if list(expression.columns) != list(metadata.sample_id):
        raise ValueError("Expression columns not aligned with metadata")
    x = expression.to_numpy(dtype=float)
    variance = x.var(axis=1)
    idx = np.argsort(variance)[-min(top_variable, len(variance)):]
    z = x[idx].T
    z = z-z.mean(axis=0)
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    points = u[:, :2]*s[:2]
    # SVD sign is arbitrary; orient reproducibly.
    for j in range(2):
        if points[np.argmax(np.abs(points[:, j])), j] < 0:
            points[:, j] *= -1
    explained = (s[:2]**2)/(s**2).sum()
    return pd.DataFrame({"sample_id": metadata.sample_id, "subtype": metadata.subtype,
                         "batch": metadata.batch, "PC1": points[:, 0], "PC2": points[:, 1]}), explained
