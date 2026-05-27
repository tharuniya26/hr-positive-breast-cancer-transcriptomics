"""Run the full observed-data analysis and generate result tables and figures."""
from pathlib import Path
import json
import pandas as pd
from src.data_loading import parse_geo,load_annotation,sha256
from src.preprocessing import select_samples,gene_matrix
from src.quality_control import qc_metrics,pca
from src.differential_expression import fit_moderated
from src.pathway_analysis import read_reactome,overrepresentation
from src import visualization as viz

raw=Path("data/raw");out=Path("results");out.mkdir(exist_ok=True)
Path("data/processed").mkdir(parents=True,exist_ok=True)
Path("data/metadata").mkdir(parents=True,exist_ok=True)
matrix,all_meta=parse_geo(raw/"GSE45827_series_matrix.txt.gz")
probe,meta=select_samples(matrix,all_meta)
annot=load_annotation(raw/"GPL570.annot.gz")
genes,mapping=gene_matrix(probe,annot)
meta.to_csv("data/metadata/sample_metadata.csv",index=False)
genes.to_csv("data/processed/gene_expression.csv.gz",compression="gzip")
mapping.to_csv("data/processed/gene_probe_mapping.csv")
qc=qc_metrics(genes,meta);qc.to_csv(out/"quality_metrics.csv",index=False)
points,explained=pca(genes,meta);points.to_csv(out/"pca_scores.csv",index=False)
de=fit_moderated(genes,meta)
de=de.merge(mapping[["probe_id","gene_id"]],left_on="gene_symbol",right_index=True,how="left")
de.to_csv(out/"differential_expression.csv",index=False)
candidate_roles={
    "ESR1":("ERα receptor","Oestrogen signalling"),"PGR":("Nuclear progesterone receptor","Progesterone and ER cross-talk"),
    "ERBB2":("HER2 receptor tyrosine kinase","Growth-factor signalling"),"MKI67":("Proliferation-associated marker","Cell cycle"),
    "CCND1":("Cyclin D1","CDK4/6–RB–E2F"),"CDK4":("Cyclin-dependent kinase 4","CDK4/6–RB–E2F"),
    "CDK6":("Cyclin-dependent kinase 6","CDK4/6–RB–E2F"),"RB1":("Retinoblastoma protein","CDK4/6–RB–E2F"),
    "FOXO1":("Forkhead transcription factor","PI3K/AKT and apoptosis"),"PIK3CA":("PI3K catalytic subunit","PI3K/AKT/mTOR"),
    "AKT1":("AKT kinase","PI3K/AKT/mTOR"),"MTOR":("mTOR kinase","PI3K/AKT/mTOR"),
    "MAPK1":("ERK2 kinase","MAPK/ERK"),"TNFSF11":("RANKL ligand","RANK/NF-κB"),
    "TNFRSF11A":("RANK receptor","RANK/NF-κB"),"WNT4":("WNT ligand","WNT/β-catenin"),
    "PGRMC1":("Progesterone receptor membrane component","Non-genomic progesterone-associated signalling"),
    "CYP19A1":("Aromatase","Oestrogen synthesis")}
markers=pd.DataFrame([{"gene_symbol":g,"biological_role":role,"pathway":pathway}
                      for g,(role,pathway) in candidate_roles.items()])
markers=markers.merge(de,on="gene_symbol",how="left")
markers["measurement_status"]=markers.probe_id.notna().map({True:"measured",False:"not in processed matrix"})
markers.to_csv(out/"marker_summary.csv",index=False)
sets=read_reactome(raw/"ReactomePathways.gmt")
enrich=overrepresentation(de,sets)
enrich.to_csv(out/"pathway_enrichment.csv",index=False)
fig=Path("figures")
viz.workflow(fig);viz.pca_plot(points,explained,fig);viz.volcano(de,fig)
viz.heatmap(genes,meta,de,fig);viz.marker_plot(genes,meta,de,fig)
viz.pathway_plot(enrich,fig);viz.conceptual(fig)
summary={"samples":len(meta),"luminal_a":int((meta.subtype=="Luminal A").sum()),
         "luminal_b":int((meta.subtype=="Luminal B").sum()),"genes_tested":len(de),
         "significant_fdr05":int((de.adjusted_p_value<.05).sum()),
         "significant_fdr05_abs_log2fc_ge_0_5":int(((de.adjusted_p_value<.05)&(de.log2_fold_change.abs()>=.5)).sum()),
         "pathways_tested":len(enrich),"pathways_fdr05":int((enrich.adjusted_p_value<.05).sum()),
         "pca_variance_pc1":float(explained[0]),"pca_variance_pc2":float(explained[1]),
         "geo_sha256":sha256(raw/"GSE45827_series_matrix.txt.gz"),
         "annotation_sha256":sha256(raw/"GPL570.annot.gz"),
         "reactome_sha256":sha256(raw/"ReactomePathways.gmt")}
(out/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
print(de.loc[de.gene_symbol.isin(["ESR1","PGR","ERBB2","MKI67","CCND1","FOXO1","WNT4","TNFSF11"]),["gene_symbol","log2_fold_change","adjusted_p_value"]].to_string(index=False))
