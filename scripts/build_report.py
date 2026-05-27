"""Generate tables, compile the XeLaTeX report, and package standalone sources."""
from pathlib import Path
import csv
import json
import re
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
TABLES = REPORTS / "tables"
STEM = "molecular_transcriptomic_analysis"


def escape(value):
    mapping = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
               "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}"}
    return "".join(mapping.get(c, c) for c in str(value))


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def probability(value):
    if not value:
        return "---"
    number = float(value)
    if number >= 0.001:
        return f"{number:.4f}"
    mantissa, exponent = f"{number:.2e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def write_table(name, text):
    (TABLES / name).write_text(text + "\n", encoding="utf-8")


def generate_tables():
    TABLES.mkdir(exist_ok=True)
    metadata = rows(ROOT / "data/metadata/sample_metadata.csv")
    batch_rows = []
    for batch in sorted({r["batch"] for r in metadata}, key=int):
        a = sum(r["batch"] == batch and r["subtype"] == "Luminal A" for r in metadata)
        b = sum(r["batch"] == batch and r["subtype"] == "Luminal B" for r in metadata)
        batch_rows.append(f"{batch} & {a} & {b} & {a+b}" + r"\\")
    write_table("batch.tex", r"""\begin{table}[!htbp]
\centering
\caption{Distribution of analysed samples across recorded batches. Each batch contains both source-labelled groups.}
\label{tab:batch}
\begin{tabular}{@{}lrrr@{}}\toprule
Batch & Luminal A & Luminal B & Total\\\midrule
""" + "\n".join(batch_rows) + r"""
\midrule Total & 29 & 30 & 59\\\bottomrule
\end{tabular}
\end{table}""")
    markers = rows(ROOT / "results/marker_summary.csv")
    main = {"ESR1", "PGR", "ERBB2", "MKI67", "CCND1", "FOXO1", "WNT4", "TNFSF11"}
    marker_rows, extended_rows = [], []
    for row in markers:
        measured = bool(row["log2_fold_change"])
        effect = f"{float(row['log2_fold_change']):+.3f}" if measured else "---"
        q = probability(row["adjusted_p_value"])
        gene = r"\gene{" + escape(row["gene_symbol"]) + "}"
        if row["gene_symbol"] in main:
            marker_rows.append(f"{gene} & {escape(row['probe_id'])} & ${effect}$ & {q}" + r"\\")
        status = "Not retained"
        if measured:
            status = "Below FDR threshold" if float(row["adjusted_p_value"]) < .05 else "Not significant"
        extended_rows.append(f"{gene} & {effect} & {q} & {status}" + r"\\")
    write_table("markers.tex", r"""\begin{table}[!htbp]
\centering\small
\caption{Principal marker results. Coefficients are adjusted Luminal B minus Luminal A log$_2$ differences. Probabilities are BH-adjusted across all 13,823 tested genes.}
\label{tab:markers}
\begin{tabular}{@{}llrr@{}}\toprule
Gene & Retained probe & Log$_2$ difference & Adjusted $p$\\\midrule
""" + "\n".join(marker_rows) + r"""
\bottomrule\end{tabular}
\end{table}""")
    write_table("extended_markers.tex", r"""\begin{longtable}{@{}lrrl@{}}
\caption{Complete marker panel from the saved analysis. ``Below FDR threshold'' refers to adjusted $p<0.05$, regardless of effect magnitude.}\label{tab:extended}\\
\toprule Gene & Log$_2$ difference & Adjusted $p$ & Measurement / test result\\\midrule
\endfirsthead
\toprule Gene & Log$_2$ difference & Adjusted $p$ & Measurement / test result\\\midrule
\endhead
""" + "\n".join(extended_rows) + r"""
\bottomrule
\end{longtable}""")
    pathway_rows = []
    for row in rows(ROOT / "results/pathway_enrichment.csv")[:8]:
        name = row["pathway"].split(" [")[0]
        fold = int(row["overlap_size"]) / (int(row["set_size"]) * int(row["hit_list_size"]) / int(row["universe_size"]))
        pathway_rows.append(f"{escape(name)} & {row['set_size']} & {row['overlap_size']} & {fold:.2f} & {probability(row['adjusted_p_value'])}" + r"\\")
    write_table("pathways.tex", r"""\begin{table}[!htbp]
\centering\small
\caption{Eight leading Reactome sets, ranked by adjusted probability. Size is the number of members in the measured universe; overlap is the number in the foreground. Enrichment is observed overlap divided by its hypergeometric expectation.}
\label{tab:pathways}
\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash}Xrrrr@{}}\toprule
Pathway & Size & Overlap & Enrichment & Adjusted $p$\\\midrule
""" + "\n".join(pathway_rows) + r"""
\bottomrule\end{tabularx}
\end{table}""")
    summary = json.loads((ROOT / "results/summary.json").read_text())
    hashes = []
    for title, key in [("GSE45827 series matrix", "geo_sha256"), ("GPL570 annotation", "annotation_sha256"), ("Reactome release 97 GMT", "reactome_sha256")]:
        digest = summary[key]
        hashes.append(r"\noindent\textbf{" + title + r"}\par{\small\ttfamily " + digest[:32] + r"\newline " + digest[32:] + r"\par}\medskip")
    write_table("hashes.tex", "\n".join(hashes))


def author_names(value):
    names = []
    for name in value.split(" and "):
        if name == "others":
            names.append("et al.")
        elif "," in name:
            surname, given = name.split(",", 1)
            initials = "".join(part[0] + "." for part in given.split())
            names.append(escape(surname + " " + initials))
        else:
            names.append(escape(name))
    return ", ".join(names)


def generate_bibliography():
    # The project bibliography uses one simple brace-delimited record per line.
    entries = []
    for line in (ROOT / "references/bibliography.bib").read_text(encoding="utf-8").splitlines():
        match = re.match(r"@\w+\{([^,]+),", line)
        if not match:
            continue
        key = match.group(1)
        fields = dict(re.findall(r"(\w+)\s*=\s*\{([^{}]*)\}", line))
        item = rf"\bibitem{{{key}}} "
        if "author" in fields:
            item += author_names(fields["author"]) + " "
        item += "(" + fields.get("year", "n.d.") + "). " + escape(fields["title"]) + ". "
        if "journal" in fields:
            item += r"\textit{" + escape(fields["journal"]) + "}"
            if "volume" in fields:
                item += " " + fields["volume"]
            if "pages" in fields:
                item += ":" + fields["pages"]
            item += ". "
        if "doi" in fields:
            item += r"\url{https://doi.org/" + fields["doi"] + "}."
        elif "url" in fields:
            item += r"\url{" + fields["url"] + "}. Accessed 16 October 2026."
        entries.append((key, item))
    cited = []
    for keys in re.findall(r"\\cite\{([^}]+)\}", (REPORTS / f"{STEM}.tex").read_text(encoding="utf-8")):
        for key in keys.split(","):
            if key not in cited:
                cited.append(key)
    lookup = dict(entries)
    missing = set(cited) - set(lookup)
    if missing:
        raise ValueError(f"Missing bibliography keys: {missing}")
    (REPORTS / "bibliography.tex").write_text(
        "\\begin{thebibliography}{99}\n\\small\\raggedright\n" + "\n\n".join(lookup[k] for k in cited) +
        "\n\\end{thebibliography}\n", encoding="utf-8")


def build():
    generate_tables()
    generate_bibliography()
    engine = shutil.which("xelatex")
    if not engine:
        raise SystemExit("XeLaTeX is required. Install MiKTeX/TeX Live or compile the source with XeLaTeX in Overleaf.")
    build_dir = REPORTS / ".build"
    build_dir.mkdir(exist_ok=True)
    version = subprocess.run([engine, "--version"], capture_output=True, text=True, check=True).stdout
    options = ["-interaction=nonstopmode", "-halt-on-error", "-file-line-error", "-no-shell-escape"]
    if "MiKTeX" in version:
        options.append("--disable-installer")
    command = [engine, *options, f"-output-directory={build_dir}", f"{STEM}.tex"]
    for iteration in range(3):
        result = subprocess.run(command, cwd=REPORTS, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=180)
        (build_dir / f"pass_{iteration+1}.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"XeLaTeX failed. See {build_dir / f'pass_{iteration+1}.txt'}\n{result.stdout[-3500:]}")
    log = (build_dir / f"{STEM}.log").read_text(encoding="utf-8", errors="replace")
    if "undefined references" in log or "Citation" in log and "undefined" in log:
        raise RuntimeError("Unresolved references remain in the LaTeX log.")
    shutil.copy2(build_dir / f"{STEM}.pdf", REPORTS / f"{STEM}.pdf")
    archive = REPORTS / "latex_source.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in [REPORTS / f"{STEM}.tex", REPORTS / "bibliography.tex", REPORTS / "LATEX_README.md"]:
            bundle.write(path, path.name)
        for path in sorted(TABLES.glob("*.tex")):
            bundle.write(path, "tables/" + path.name)
        for path in sorted((ROOT / "figures").glob("*.png")):
            bundle.write(path, "figures/" + path.name)
        bundle.write(ROOT / "references/bibliography.bib", "bibliography.bib")
    print(f"PDF: {REPORTS / f'{STEM}.pdf'}")
    print(f"LaTeX source archive: {archive}")
    warnings = [line for line in log.splitlines() if "Overfull" in line or "Warning:" in line]
    print("Layout warnings:", "\n".join(warnings) if warnings else "none")


if __name__ == "__main__":
    build()
