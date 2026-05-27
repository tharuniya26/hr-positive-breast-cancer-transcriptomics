# Project quality audit — 16 October 2026

Status terms: **passed** means directly checked in this workspace; **not verified** means the evidence is incomplete; **failed** identifies an unmet target.

| Check | Status | Evidence or remaining issue |
|---|---|---|
| Distinct computational project | Passed | GEO matrix, new code, executed analyses and report; prior coursework used only for background. |
| Original documents excluded from Git | Passed | Both in `private_sources/`; `git check-ignore` confirmed exclusion. Public reports contain no cover sheets, student numbers or declaration pages. |
| Presentation attribution | Passed | Private inventory lists all five group contributors; public source review acknowledges group provenance. No third-party source figure copied. |
| Source figure inventory | Passed | 11 DOCX media objects and 27 PDF image instances extracted into the private archive; visible figure captions and reuse restrictions recorded. No source image published. |
| Real dataset and provenance | Passed | GSE45827, GPL570 and downloaded file hashes recorded in `results/summary.json` and `docs/data_provenance.md`. |
| Sample annotations and group counts | Passed | Matrix alignment checks; 29 Luminal A and 30 Luminal B source-labelled tumours. No independent receptor IHC or PAM50 labels claimed. |
| Platform-appropriate preprocessing | Passed | GEO-reported GC-RMA log2 values used without RNA-seq normalisation; ambiguous probe mappings removed. |
| Batch-adjusted differential expression | Passed | Fixed-effects model, variance moderation and BH FDR executed for 13,823 genes. |
| Numerical equivalence to Bioconductor limma | Not verified | Python empirical-Bayes implementation has unit tests but no external limma cross-validation. This limits confidence in exact p-values. |
| Multiple-testing and measured-gene ORA universe | Passed | BH correction in code; 13,823 tested genes form the pathway universe. |
| Exact Reactome release | Passed | Official V97 GMT downloaded, extracted and SHA-pinned; checksum also matches an independent V97 manifest. Enrichment remains exploratory. |
| Figures based on actual results | Passed | Figures 2–6 generated from executed tables/matrix. Figures 1, 7 and 8 clearly identified as conceptual. PNG and SVG generated. |
| Notebook execution | Passed | Six notebooks executed sequentially with `nbclient` from repository root on 16 October 2026, including within the isolated environment. All code cells have execution counts and zero error outputs. Windows emitted a nonfatal ZMQ teardown warning after the last kernel. |
| Automated tests | Passed | Clean-environment `python -m pytest -q`: 8 passed. Tests cover real input dimensions, duplicated IDs, missing data, PCA determinism, differential-model structure, both Reactome GMT formats, ORA universe, explicit missing-marker status and figure output. |
| Report PDF and editable source | Passed | Eight-page PDF rendered from `.qmd` by `scripts/build_report.py`; text extraction and page count checked. |
| Bibliographic authenticity | Passed for cited core references | Dataset article, subtype, PR, RANKL, model and resistance DOIs checked against GEO/publisher or PubMed records; inconsistent source bibliography entries omitted. Not every historical-source citation was exhaustively audited because it was not reused. |
| Fresh isolated environment | Passed | A new `.venv` was installed from `requirements.txt`; download verification, full pipeline, PDF build, eight tests and six notebooks completed there. The check caught and resolved a missing `statsmodels` dependency. Tested with Python 3.14.2, NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1, statsmodels 0.15.0 and Matplotlib 3.11.2. |
| Independent biological/clinical validation | Not verified | One observational cohort; no IHC, outcome, treatment or external validation dataset. No clinical claims made. |
| Data redistribution permissions | Not verified | Public GEO data remain download-only and Git-ignored; exact upstream redistribution terms were not established. |

## Reviewer synthesis

- **Molecular oncology:** PGR, MKI67, ERBB2 and cell-cycle interpretations are plausible as cohort associations; no mechanistic or treatment-response inference is supported.
- **Bioinformatics:** Input validation, hashes, probe mapping, batch covariates and reusable modules are present. The probe annotation is from 2016 and one-probe selection omits alternative signals.
- **Transcriptomics:** Source labels support a luminal contrast, but a classifier or individual receptor assay is unavailable. The original study had a different research aim.
- **Biostatistics:** BH and a measured-gene universe are implemented. The independent moderated model needs a limma numerical comparison; unmeasured confounding and cohort selection remain.
- **Software and admissions review:** The repository has executable commands, six notebooks, tests, report source and honest limitations. It is a portfolio analysis, not a peer-reviewed or clinically validated discovery.

## Next validation priorities

1. Cross-validate the moderated coefficients, degrees of freedom and p-values against Bioconductor limma on the same matrix and design.
2. Preserve the official V97 GMT checksum and rerun enrichment if the Reactome release is changed.
3. Replicate marker and cell-cycle patterns in an independent cohort with documented IHC or PAM50 labels and appropriate covariates.
