# Limitations

1. GEO source labels are not independent receptor measurements or reproducible PAM50 calls. The apparent Luminal B pattern may partly reflect how samples were labelled.
2. The original study sought a different biological question, so this contrast is secondary and exploratory.
3. No age, grade, tumour purity, treatment, HER2 IHC, patient pairing or outcome covariates were available in the selected matrix; confounding and heterogeneity remain.
4. GEO reports GC-RMA normalization and a mixed-model procedure, but precise transformations of stored values are not fully auditable here. Batch is modelled as a fixed effect; residual batch effects can persist.
5. GPL570 annotation is dated 2016; gene identifiers and probe specificity can change. Selecting one high-intensity probe per symbol loses alternative isoform information.
6. The independent empirical-Bayes implementation needs formal numerical cross-validation against Bioconductor limma. It is appropriate for exploratory software research, not a validated clinical tool.
7. Reactome release 97 gene sets overlap extensively and can change in future releases. ORA does not infer pathway activity or direction.
8. mRNA abundance does not measure ER/PR protein, phosphorylation, kinase activity, ESR1 mutations, drug sensitivity or prognosis.
9. The observations come from one microarray cohort. There is no independent validation cohort or prospective clinical assessment.
