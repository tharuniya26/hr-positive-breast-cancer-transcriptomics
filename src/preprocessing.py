"""Select source-labelled luminal tumours and map unambiguous probes to genes."""
import numpy as np
import pandas as pd

GROUPS = ("Luminal A", "Luminal B")

def select_samples(expression: pd.DataFrame, metadata: pd.DataFrame):
    if metadata.sample_id.duplicated().any() or expression.index.duplicated().any():
        raise ValueError("Duplicate sample or probe identifier")
    selected = metadata.loc[metadata.subtype.isin(GROUPS)].copy()
    if selected.shape[0] != 59 or selected.subtype.value_counts().to_dict() != {"Luminal B": 30, "Luminal A": 29}:
        raise ValueError("Unexpected eligible sample counts")
    if not selected.source.str.contains("Tumor Sample").all():
        raise ValueError("Non-tumour among eligible samples")
    if not set(selected.sample_id).issubset(expression.columns):
        raise ValueError("Missing expression columns")
    selected = selected.sort_values("sample_id").reset_index(drop=True)
    data = expression.loc[:, selected.sample_id].copy()
    if data.isna().any().any() or not np.isfinite(data.to_numpy(dtype=float)).all():
        raise ValueError("Missing or non-finite expression values")
    return data, selected

def gene_matrix(expression: pd.DataFrame, annotation: pd.DataFrame):
    common = expression.index.intersection(annotation.index)
    mapping = annotation.loc[common].copy()
    # Ambiguous symbols in GPL annotations use /// delimiters. Drop them.
    mapping = mapping.loc[mapping.gene_symbol.notna()]
    mapping = mapping.loc[~mapping.gene_symbol.astype(str).str.contains(r"///|//|;|,|\s", regex=True)]
    mapping = mapping.loc[mapping.gene_symbol.astype(str).str.fullmatch(r"[A-Za-z0-9_.-]+")]
    measured = expression.loc[mapping.index]
    mapping = mapping.loc[measured.index]
    mapping["mean_intensity"] = measured.mean(axis=1)
    # One probe per symbol, chosen without inspecting group contrasts.
    mapping = mapping.reset_index().sort_values(["gene_symbol", "mean_intensity", "probe_id"],
                                                ascending=[True, False, True]).set_index("probe_id")
    mapping = mapping.drop_duplicates("gene_symbol")
    selected = measured.loc[mapping.index].copy()
    selected.index = mapping.gene_symbol.to_numpy()
    selected.index.name = "gene_symbol"
    mapping = mapping.reset_index().set_index("gene_symbol")
    return selected, mapping
