"""Reactome over-representation analysis using the measured-gene universe."""
from pathlib import Path
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def read_reactome(path: Path) -> dict[str, set[str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    sets = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        cells = line.split("\t")
        if len(cells) >= 4:
            # Reactome v97 uses name, stable ID, genes. Older exports inserted
            # a literal third field, "Reactome Pathway".
            start = 3 if cells[2] == "Reactome Pathway" else 2
            sets[cells[0]+" ["+cells[1]+"]"] = set(cells[start:])
    return sets

def overrepresentation(de: pd.DataFrame, sets: dict[str, set[str]], min_size=10, max_size=500):
    universe = set(de.gene_symbol)
    hits = set(de.loc[(de.adjusted_p_value < 0.05) & (de.log2_fold_change.abs() >= 0.5), "gene_symbol"])
    rows = []
    for name, members in sets.items():
        genes = members & universe
        if not min_size <= len(genes) <= max_size:
            continue
        overlap = genes & hits
        p = hypergeom.sf(len(overlap)-1, len(universe), len(genes), len(hits))
        rows.append({"pathway": name, "universe_size": len(universe), "set_size": len(genes),
                     "hit_list_size": len(hits), "overlap_size": len(overlap),
                     "overlap_genes": ";".join(sorted(overlap)), "p_value": p})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["adjusted_p_value"] = multipletests(out.p_value, method="fdr_bh")[1]
    return out.sort_values("adjusted_p_value").reset_index(drop=True)
