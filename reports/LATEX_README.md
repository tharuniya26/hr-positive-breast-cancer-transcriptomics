# Detailed LaTeX research report

The detailed report is `molecular_transcriptomic_analysis.tex`.
The PDF is compiled with XeLaTeX. The original `.qmd` is retained as an earlier,
shorter narrative and is no longer the input to the report builder.

## Typography and layout

A4 paper, 25 mm margins, 11 pt Times New Roman body text, restrained headings,
numbered equations, booktabs tables, numbered citations, linked references,
and page headers. Where Times New Roman is unavailable (including Overleaf),
the source uses TeX Gyre Termes, a compatible professional serif family.
Figures are the original analysis outputs; no scientific results are simulated.

## Build from the repository

Install MiKTeX or TeX Live with XeLaTeX and these standard packages:
`fontspec`, `geometry`, `amsmath`, `amssymb`, `graphicx`, `booktabs`, `tabularx`,
`longtable`, `caption`, `fancyhdr`, `xcolor`, and `hyperref`.
Then run from the repository root:

```sh
python -m scripts.build_report
```

This updates `tables/*.tex` from the saved result files, updates
`bibliography.tex` from `references/bibliography.bib`, compiles three passes,
and creates `latex_source.zip`. Temporary files remain in the ignored `.build/`.
The bibliography uses numbered `thebibliography` entries; BibTeX/Biber is not
required. Edit the simple, one-record-per-line `.bib` file and rebuild when
changing references.

## Overleaf or standalone editing

Upload `latex_source.zip` as a new Overleaf project. Set the compiler to
**XeLaTeX** and the main document to `molecular_transcriptomic_analysis.tex`.
The archive includes figures, generated tables and bibliography; it does not
require Python or private/raw data to compile. Locally, extract the archive and
run XeLaTeX on the main document three times to resolve contents and citations.

The report preserves the actual analysis date, 16 October 2026. Commit timestamps
are independent of the dates of data retrieval, analysis and report preparation.
