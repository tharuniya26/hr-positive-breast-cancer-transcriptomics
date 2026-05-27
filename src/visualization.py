"""Figures from observed data and separate conceptual schematics."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {"Luminal A": "#0072B2", "Luminal B": "#D55E00"}

def save(fig, folder: Path, stem: str):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder/f"{stem}.png", dpi=250, bbox_inches="tight")
    fig.savefig(folder/f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)

def workflow(folder):
    stages = ["GSE45827\npublic matrix", "Sample\nvalidation", "QC & PCA",
              "Batch-adjusted\nDE", "Reactome\nenrichment", "Molecular\ninterpretation"]
    fig, ax = plt.subplots(figsize=(13, 2.5))
    for i, label in enumerate(stages):
        ax.text(i, 0.5, label, ha="center", va="center", fontsize=10,
                bbox=dict(boxstyle="round,pad=.5", fc="#E8F1F8", ec="#487899"))
        if i < len(stages)-1:
            ax.annotate("", xy=(i+.8, .5), xytext=(i+.38, .5),
                        arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.set_xlim(-.5, 5.5); ax.set_ylim(0, 1); ax.axis("off")
    save(fig, folder, "01_workflow")

def pca_plot(points, explained, folder):
    fig, ax = plt.subplots(figsize=(7, 5))
    for label, group in points.groupby("subtype"):
        ax.scatter(group.PC1, group.PC2, label=f"{label} (n={len(group)})",
                   c=COLORS[label], s=38, alpha=.8, edgecolor="white", linewidth=.3)
    ax.set(xlabel=f"PC1 ({explained[0]*100:.1f}% variance)",
           ylabel=f"PC2 ({explained[1]*100:.1f}% variance)",
           title="PCA of 3,000 most variable mapped genes")
    ax.legend(frameon=False); fig.tight_layout(); save(fig, folder, "02_pca")

def volcano(de, folder):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    x=de.log2_fold_change.to_numpy(); y=-np.log10(de.adjusted_p_value.clip(lower=1e-300))
    sig=de.significant_fdr05.to_numpy() & (np.abs(x)>=.5)
    ax.scatter(x[~sig], y[~sig], s=4, color="#AAAAAA", alpha=.45, rasterized=True)
    ax.scatter(x[sig], y[sig], s=6, color="#7B3294", alpha=.65, rasterized=True)
    for gene in ["ESR1", "PGR", "MKI67", "CCND1", "FOXO1", "ERBB2"]:
        q=de.loc[de.gene_symbol==gene]
        if len(q):
            row=q.iloc[0]; ax.annotate(gene, (row.log2_fold_change, -np.log10(max(row.adjusted_p_value,1e-300))), fontsize=8)
    ax.axhline(-np.log10(.05), color="#555555", lw=.7, ls="--")
    ax.axvline(-.5, color="#555555", lw=.7, ls="--");ax.axvline(.5, color="#555555", lw=.7, ls="--")
    ax.set(xlabel="Adjusted log2 difference: Luminal B − Luminal A", ylabel="−log10(FDR)",
           title="Gene-level differential expression")
    fig.tight_layout();save(fig,folder,"03_volcano")

def heatmap(expression, metadata, de, folder):
    genes=de.loc[(de.adjusted_p_value<.05)&(de.log2_fold_change.abs()>=.5)].head(30).gene_symbol
    if len(genes)<2:return
    ordered=metadata.sort_values(["subtype","sample_id"])
    z=expression.loc[genes,ordered.sample_id].to_numpy()
    z=(z-z.mean(axis=1,keepdims=True))/np.maximum(z.std(axis=1,keepdims=True),1e-10)
    fig,ax=plt.subplots(figsize=(10,7));im=ax.imshow(z,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
    ax.set_yticks(range(len(genes)),genes,fontsize=7);ax.set_xticks([])
    ax.axvline((ordered.subtype=="Luminal A").sum()-.5,color="black",lw=1)
    ax.set(xlabel="Samples ordered by source subtype",title="Top 30 DE genes; row z scores")
    fig.colorbar(im,ax=ax,label="Within-gene z score",shrink=.7)
    fig.tight_layout();save(fig,folder,"04_heatmap")

def marker_plot(expression, metadata, de, folder):
    genes=[x for x in ["ESR1","PGR","ERBB2","MKI67","CCND1","FOXO1","WNT4","TNFSF11"] if x in expression.index]
    fig,axes=plt.subplots(2,4,figsize=(12,6));rng=np.random.default_rng(24)
    for ax,gene in zip(axes.flat,genes):
        for i,label in enumerate(COLORS):
            ids=metadata.loc[metadata.subtype==label,"sample_id"]
            vals=expression.loc[gene,ids].to_numpy()
            ax.scatter(i+rng.uniform(-.13,.13,len(vals)),vals,color=COLORS[label],s=11,alpha=.6)
            ax.errorbar(i,np.mean(vals),yerr=1.96*np.std(vals,ddof=1)/np.sqrt(len(vals)),
                        fmt="_",color="black",capsize=4,markersize=14)
        row=de.loc[de.gene_symbol==gene].iloc[0]
        ax.set(title=f"{gene} | FDR {row.adjusted_p_value:.2g}",xticks=[0,1],xticklabels=["A","B"])
    for ax in axes.flat[len(genes):]:ax.axis("off")
    fig.suptitle("GCRMA log2 signal; black bars show means and 95% normal CIs")
    fig.tight_layout();save(fig,folder,"05_markers")

def pathway_plot(enrichment, folder):
    view=enrichment.loc[enrichment.adjusted_p_value<.05].head(12).iloc[::-1]
    fig,ax=plt.subplots(figsize=(9,5))
    if len(view):
        labels=view.pathway.str.replace(r" \[.*", "", regex=True).str.slice(0,65)
        ax.barh(labels,-np.log10(view.adjusted_p_value.clip(lower=1e-300)),color="#009E73")
    else:ax.text(.5,.5,"No pathways met FDR < 0.05",ha="center",transform=ax.transAxes)
    ax.set(xlabel="−log10(FDR)",title="Reactome over-representation; non-directional")
    fig.tight_layout();save(fig,folder,"06_pathways")

def conceptual(folder):
    fig,ax=plt.subplots(figsize=(9,5))
    nodes={"Oestrogen / ERα":(.15,.78),"Progesterone / PR":(.15,.42),"RTK / PI3K–AKT":(.49,.85),
           "MAPK–ERK":(.49,.6),"RANKL / WNT4":(.49,.3),"CCND1–CDK4/6":(.77,.63),"RB–E2F cycle":(.77,.3)}
    for name,(x,y) in nodes.items():ax.text(x,y,name,ha="center",va="center",fontsize=9,bbox=dict(boxstyle="round",fc="#E8F1F8",ec="#487899"))
    for a,b in [("Oestrogen / ERα","CCND1–CDK4/6"),("Progesterone / PR","RANKL / WNT4"),("RTK / PI3K–AKT","CCND1–CDK4/6"),("MAPK–ERK","CCND1–CDK4/6"),("CCND1–CDK4/6","RB–E2F cycle")]:
        ax.annotate("",xy=nodes[b],xytext=nodes[a],
                    arrowprops=dict(arrowstyle="->",color="#555555",lw=1.2,
                                    shrinkA=38,shrinkB=38))
    ax.text(.5,.03,"Literature-based schematic; arrows are context-dependent and do not show measured activity.",ha="center",fontsize=9)
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");save(fig,folder,"07_signalling")
    fig,ax=plt.subplots(figsize=(8,4.5))
    drugs=[("Tamoxifen","ER modulation","Approved"),("Fulvestrant","ER degradation","Approved"),
           ("Anastrozole","Aromatase inhibition","Approved"),("CDK4/6 inhibitors","CDK4/6 inhibition","Approved"),
           ("Docetaxel","Microtubule stabilization","Approved"),("Onapristone","PR antagonism","Investigational")]
    for i,(drug,target,status) in enumerate(drugs):
        y=1-i/7; color="#0072B2" if status=="Approved" else "#D55E00"
        ax.text(.05,y,drug,color=color,fontsize=11);ax.annotate("",xy=(.49,y),xytext=(.36,y),arrowprops=dict(arrowstyle="->",color=color));ax.text(.52,y,target,fontsize=10)
    ax.text(.05,.02,"Blue: approved drug class or agent; orange: investigational in breast cancer.",fontsize=9)
    ax.set_xlim(0,1);ax.set_ylim(0,1.08);ax.axis("off");save(fig,folder,"08_therapies")
