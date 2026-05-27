"""Generate six sequential, readable notebooks from shared analysis functions."""
from pathlib import Path
import nbformat as nbf

root=Path("notebooks");root.mkdir(exist_ok=True)
intro="Run from the repository root after `python -m scripts.download_data` and `python -m scripts.run_analysis`. These notebooks use reusable functions from `src/`. Observations come from GSE45827, not from the historical coursework."

plans=[
("01_dataset_exploration","Dataset exploration","Which GEO columns are suitable for a Luminal A/B comparison?",
 "Inspect source labels, eligibility, and matrix scale. GEO reports GC-RMA normalized log2 signals.",
 "from pathlib import Path\nfrom src.data_loading import parse_geo\nfrom src.preprocessing import select_samples\nexpression, all_metadata = parse_geo(Path('data/raw/GSE45827_series_matrix.txt.gz'))\nselected, metadata = select_samples(expression, all_metadata)\nprint(expression.shape, selected.shape)\nprint(all_metadata.subtype.value_counts())\nprint(metadata.groupby(['subtype', 'batch']).size())\nprint(selected.stack().describe())",
 "The eligible set contains 29 source-labelled Luminal A and 30 Luminal B primary tumours. Labels are source annotations, not independently verified receptor status.",
 "The matrix lacks individual IHC, PAM50 method and outcome covariates."),
("02_quality_control","Quality control","Do the selected intensities show missingness, unusual distributions, or broad variation?",
 "Map probes to genes without using group differences; inspect per-sample medians and PCA.",
 "from pathlib import Path\nfrom src.data_loading import parse_geo, load_annotation\nfrom src.preprocessing import select_samples, gene_matrix\nfrom src.quality_control import qc_metrics, pca\nfrom src.visualization import pca_plot\nexpression, all_metadata = parse_geo(Path('data/raw/GSE45827_series_matrix.txt.gz'))\nselected, metadata = select_samples(expression, all_metadata)\ngenes, mapping = gene_matrix(selected, load_annotation(Path('data/raw/GPL570.annot.gz')))\nmetrics = qc_metrics(genes, metadata)\npoints, variance = pca(genes, metadata)\nprint(metrics[['subtype','batch','median_intensity','iqr_intensity','missing_fraction']].describe(include='all'))\nprint('PCA variance:', variance)\npca_plot(points, variance, Path('figures'))\nfrom IPython.display import display, Image\ndisplay(Image(filename='figures/02_pca.png'))",
 "PCA displays broad sample variation; its axes are not proof of biological subtyping. No sample is excluded after viewing the plot.",
 "Potential batch and tumour-purity variation remain; source preprocessing is not fully auditable."),
("03_differential_expression","Differential expression","Which mapped genes differ between source-labelled groups after batch adjustment?",
 "Use a fixed-effects subtype-plus-batch model with empirical-Bayes variance moderation and BH FDR.",
 "from pathlib import Path\nfrom src.data_loading import parse_geo, load_annotation\nfrom src.preprocessing import select_samples, gene_matrix\nfrom src.differential_expression import fit_moderated\nfrom src.visualization import volcano, heatmap\nexpression, all_metadata = parse_geo(Path('data/raw/GSE45827_series_matrix.txt.gz'))\nselected, metadata = select_samples(expression, all_metadata)\ngenes, mapping = gene_matrix(selected, load_annotation(Path('data/raw/GPL570.annot.gz')))\nde = fit_moderated(genes, metadata)\nprint('Genes tested:', len(de), 'FDR<0.05:', sum(de.significant_fdr05))\nprint(de.head(10).to_string(index=False))\nvolcano(de, Path('figures'))\nheatmap(genes, metadata, de, Path('figures'))\nfrom IPython.display import display, Image\ndisplay(Image(filename='figures/03_volcano.png'))\ndisplay(Image(filename='figures/04_heatmap.png'))",
 "Positive coefficients indicate higher log2 signal in labelled Luminal B. The heatmap is scaled within each gene.",
 "This independent limma-style implementation is not formally cross-validated against Bioconductor limma; unmeasured covariates remain."),
("04_pathway_analysis","Pathway analysis","Are descriptive differential genes over-represented in curated pathways?",
 "Use Reactome gene sets, a tested-gene universe, hypergeometric p-values and BH correction.",
 "from pathlib import Path\nimport pandas as pd\nfrom src.pathway_analysis import read_reactome, overrepresentation\nfrom src.visualization import pathway_plot\nde = pd.read_csv('results/differential_expression.csv')\nsets = read_reactome(Path('data/raw/ReactomePathways.gmt'))\nenrichment = overrepresentation(de, sets)\nprint(enrichment[['pathway','overlap_size','adjusted_p_value']].head(12).to_string(index=False))\npathway_plot(enrichment, Path('figures'))\nfrom IPython.display import display, Image\ndisplay(Image(filename='figures/06_pathways.png'))",
 "Cell-cycle sets are prominent; many Reactome sets overlap and the result has no directionality.",
 "Reactome release 97 sets overlap extensively, and ORA cannot measure pathway activation."),
("05_hormone_receptor_gene_analysis","Hormone receptor and candidate genes","How do selected endocrine and signalling transcripts vary across groups?",
 "Display individual normalized signals and link each candidate to adjusted statistics.",
 "from pathlib import Path\nimport pandas as pd\nfrom src.visualization import marker_plot\ngenes = pd.read_csv('data/processed/gene_expression.csv.gz', index_col=0)\nmetadata = pd.read_csv('data/metadata/sample_metadata.csv')\nde = pd.read_csv('results/differential_expression.csv')\ncandidates = ['ESR1','PGR','ERBB2','MKI67','CCND1','CDK4','CDK6','RB1','FOXO1','PIK3CA','AKT1','MTOR','MAPK1','TNFSF11','TNFRSF11A','WNT4','PGRMC1','CYP19A1']\nprint(de.loc[de.gene_symbol.isin(candidates),['gene_symbol','log2_fold_change','adjusted_p_value']].sort_values('gene_symbol').to_string(index=False))\nprint('Not mapped:', sorted(set(candidates)-set(genes.index)))\nmarker_plot(genes, metadata, de, Path('figures'))\nfrom IPython.display import display, Image\ndisplay(Image(filename='figures/05_markers.png'))",
 "PGR marks PR biology, ESR1 marks ERα, ERBB2 marks HER2, MKI67 marks proliferation, CCND1/CDK4/6/RB1 mark cell-cycle control, FOXO1 and PI3K/AKT/MTOR mark growth signalling, and TNFSF11/TNFRSF11A/WNT4 mark possible paracrine programmes. Report each observed coefficient and FDR rather than assuming candidates are differential.",
 "Microarray transcript abundance does not measure protein abundance, receptor activity, phosphorylation, cell type or drug response. Missing platform genes are not biologically absent."),
("06_results_interpretation","Results interpretation","What conclusions follow from the observed comparison?",
 "Read executed summary and leading tables. Separate cohort findings from literature mechanisms.",
 "import json, pandas as pd\nfrom pathlib import Path\nsummary=json.loads(Path('results/summary.json').read_text())\nde=pd.read_csv('results/differential_expression.csv')\npathways=pd.read_csv('results/pathway_enrichment.csv')\nprint(summary)\nprint(de.loc[de.gene_symbol.isin(['ESR1','PGR','ERBB2','MKI67','CCND1','FOXO1','WNT4','TNFSF11']),['gene_symbol','log2_fold_change','adjusted_p_value']].to_string(index=False))\nprint(pathways[['pathway','overlap_size','adjusted_p_value']].head(5).to_string(index=False))",
 "In this source-labelled cohort, PGR is lower and MKI67/ERBB2 are higher in B, while ESR1 shows no significant contrast. Mitotic pathway over-representation is compatible with proliferation-associated differences. These are observational associations, not a validated biomarker or mechanistic activity assay.",
 "Independent subtype verification, external cohort replication, clinical covariates and numerical comparison with limma are outstanding.")]

for stem,title,question,methods,code,interpretation,limitations in plans:
    nb=nbf.v4.new_notebook(metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
    nb.cells=[nbf.v4.new_markdown_cell(f"# {title}\n\n{intro}"),
              nbf.v4.new_markdown_cell(f"## Research question\n\n{question}\n\n## Data and methods\n\n{methods}"),
              nbf.v4.new_code_cell(code),
              nbf.v4.new_markdown_cell(f"## Biological interpretation\n\n{interpretation}\n\n## Limitations\n\n{limitations}")]
    nbf.write(nb,root/(stem+".ipynb"))
print("Generated",len(plans),"notebooks")
