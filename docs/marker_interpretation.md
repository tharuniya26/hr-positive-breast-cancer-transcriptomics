# Candidate marker interpretation

All effects below are batch-adjusted **Luminal B minus Luminal A** log2 signal differences in GSE45827. FDR is Benjamini–Hochberg across 13,823 genes. Gene IDs and selected probes come from the GPL570 annotation; exact rows are in [`results/marker_summary.csv`](../results/marker_summary.csv). Candidates were chosen from biological background before interpreting differential-expression results. A missing row does not mean absent expression in breast tissue.

| Gene | Biological role and pathway | Log2 difference | FDR | Reading of this cohort |
|---|---|---:|---:|---|
| *ESR1* | ERα; oestrogen transcription | −0.228 | 0.673 | No significant difference; protein activity unknown |
| *PGR* | PR; progesterone and ER cross-talk | −4.993 | 1.57×10⁻⁷ | Lower transcript in B; PR IHC unavailable |
| *ERBB2* | HER2; receptor tyrosine kinase | +3.350 | 3.65×10⁻¹³ | Higher transcript in B; HER2 IHC unavailable |
| *MKI67* | Proliferation marker; cell cycle | +2.716 | 7.82×10⁻¹¹ | Higher transcript in B; Ki-67 protein unknown |
| *CCND1* | Cyclin D1; CDK4/6–RB–E2F | +0.165 | 0.640 | No significant difference |
| *CDK4* | CDK4 kinase; cell cycle | +0.216 | 0.166 | No significant difference |
| *CDK6* | CDK6 kinase; cell cycle | −0.255 | 0.559 | No significant difference |
| *RB1* | Retinoblastoma protein; E2F control | −0.058 | 0.755 | No significant difference; phosphorylation unknown |
| *FOXO1* | Forkhead factor; AKT/apoptosis | −0.453 | 0.0409 | Modest lower transcript; nuclear localisation unknown |
| *PIK3CA* | PI3K catalytic subunit | −0.276 | 0.367 | No significant difference; mutation/activity unknown |
| *AKT1* | AKT kinase; PI3K/AKT/mTOR | +0.609 | 0.00543 | Higher transcript; kinase phosphorylation unknown |
| *MTOR* | mTOR kinase | +0.109 | 0.292 | No significant difference; complex activity unknown |
| *MAPK1* | ERK2 kinase; MAPK/ERK | +0.028 | 0.867 | No significant difference; phosphorylation unknown |
| *TNFSF11* | RANKL ligand; RANK/NF-κB | −0.149 | 0.769 | No significant difference; paracrine secretion unknown |
| *TNFRSF11A* | RANK receptor; RANK/NF-κB | −0.058 | 0.923 | No significant difference; cell composition unknown |
| *WNT4* | WNT ligand; WNT/β-catenin | +0.548 | 0.404 | No significant difference; pathway activity unknown |
| *PGRMC1* | Membrane-associated progesterone component | −0.126 | 0.525 | No significant difference; receptor function unknown |
| *CYP19A1* | Aromatase; oestrogen synthesis | — | — | GPL570 has probes, but none remain in the processed GSE45827 matrix, consistent with the depositor's low-signal filter |

The last explanation is an inference from the GPL570 annotation, the series matrix and GEO's documented low-signal filtering; the exact reason each CYP19A1 probe was omitted cannot be proved from the processed file alone. No highlighted gene is a validated new biomarker here. Bulk-tumour mRNA may reflect cancer cells, stroma and infiltrating cells and cannot establish protein abundance, activity, endocrine response or causal mechanisms. See [molecular background](molecular_background.md) for the literature basis and [limitations](limitations.md) for study-level uncertainties.
