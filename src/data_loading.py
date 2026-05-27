"""Read the GEO matrix and GPL570 annotation with strict alignment checks."""
from __future__ import annotations
import csv
import gzip
import hashlib
from pathlib import Path
import pandas as pd

MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45827/matrix/GSE45827_series_matrix.txt.gz"
ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_geo(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    fields: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                row = next(csv.reader([line], delimiter="\t"))
                fields.setdefault(row[0], []).append(row[1:])
        else:
            raise ValueError("GEO table marker missing")
        matrix = pd.read_csv(stream, sep="\t", index_col=0, low_memory=False)
    matrix = matrix.loc[matrix.index != "!series_matrix_table_end"]
    matrix.index.name = "probe_id"
    matrix.columns = matrix.columns.str.strip('"')
    ids = fields["!Sample_geo_accession"][0]
    if len(ids) != matrix.shape[1] or list(matrix.columns) != ids or len(ids) != len(set(ids)):
        raise ValueError("GEO sample identifiers do not align")
    char = fields["!Sample_characteristics_ch1"]
    records = []
    for i, accession in enumerate(ids):
        label = char[1][i].removeprefix("tumor subtype: ")
        batch = char[2][i].removeprefix("batch: ")
        records.append({"sample_id": accession, "title": fields["!Sample_title"][0][i],
                        "source": fields["!Sample_source_name_ch1"][0][i],
                        "subtype": label, "batch": batch,
                        "hybridization": char[3][i].removeprefix("hybridation: ")})
    metadata = pd.DataFrame(records)
    return matrix.apply(pd.to_numeric, errors="raise"), metadata

def load_annotation(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        table = pd.read_csv(stream, sep="\t", skiprows=27, comment="!", low_memory=False)
    table = table.rename(columns={"ID": "probe_id", "Gene symbol": "gene_symbol", "Gene ID": "gene_id"})
    required = {"probe_id", "gene_symbol", "gene_id"}
    if not required.issubset(table):
        raise ValueError(f"Missing annotation fields: {required-set(table)}")
    table = table[["probe_id", "gene_symbol", "gene_id"]].copy()
    table["probe_id"] = table["probe_id"].astype(str)
    return table.drop_duplicates("probe_id").set_index("probe_id")
