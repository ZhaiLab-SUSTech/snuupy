import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
import numpy as np
import scanpy as sc
import tempTools

# compare nanopore mtx with illumina mtx
nanoporeExpressionPath = "generateMtx/NanoporeMultiMat"
illuminaExpressionPath = "filtered_feature_bc_matrix.h5"

illuminaAdata = sc.read_10x_h5(illuminaExpressionPath)
nanoporeAdata = sc.read_10x_mtx(nanoporeExpressionPath)
nanoporeAdata = nanoporeAdata[:, nanoporeAdata.var.index.str.find("_") == -1]

sc.pp.filter_genes(illuminaAdata, min_cells=0)
sc.pp.filter_cells(illuminaAdata, min_genes=0)
sc.pp.filter_cells(illuminaAdata, min_counts=0)
sc.pp.filter_genes(nanoporeAdata, min_cells=0)
sc.pp.filter_cells(nanoporeAdata, min_genes=0)
sc.pp.filter_cells(nanoporeAdata, min_counts=0)

mergedAdataExpression = illuminaAdata.obs.merge(
    nanoporeAdata.obs,
    left_index=True,
    right_index=True,
    suffixes=["_illumina", "_nanopore"],
)

fig, ax = plt.subplots(figsize=(2.5, 2.5))  # Figure 3c
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
sns.scatterplot(
    "n_counts_nanopore",
    "n_counts_illumina",
    data=mergedAdataExpression,
    alpha=0.05,
    color="#000000",
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(2000))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2000))
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.set_ylabel("Illumina", size=12, fontweight="bold")
ax.set_xlabel("Nanopore", size=12, fontweight="bold")
ax.set_title("UMI counts per nucleus", fontweight="bold", size=12)

fig, ax = plt.subplots(figsize=(2.5, 2.5))  # Figure 3c
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
sns.scatterplot(
    "n_genes_nanopore",
    "n_genes_illumina",
    data=mergedAdataExpression,
    alpha=0.05,
    color="#000000",
    ax=ax,
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.set_ylabel("Illumina", size=12, fontweight="bold")
ax.set_xlabel("Nanopore", size=12, fontweight="bold")
ax.set_title("Gene counts per nucleus", fontweight="bold", size=12)

## get combine umap
INPUT_H5 = "filtered_feature_bc_matrix.h5"
nucleiAdata = sc.read_10x_h5(INPUT_H5, genome=None, gex_only=True)
tempTools.plotCellScatter(nucleiAdata)  # plot basic information
nucleiAdata = tempTools.filterAd(nucleiAdata, 2300, 350, 0.01)

nanoporeAdata = sc.read_10x_mtx("generateMtx/NanoporeMultiMat")
tempTools.plotCellScatter(nanoporeAdata)
nanoporeAdata = nanoporeAdata[
    nanoporeAdata.obs.index.isin(nucleiAdata.obs.index),
    ~nanoporeAdata.var.index.str.contains("_"),
]  # only use abundance info
combineAdata = tempTools.combineAdataUseScanorama(
    [nucleiAdata, nanoporeAdata], "batch", ["Illumina", "Nanopore"]
)

fig, ax = plt.subplots(figsize=(5.8, 3.5))  # figure 3e
sc.pl.umap(
    combineAdata,
    color="batch",
    title="",
    show=False,
    legend_fontsize=10,
    ax=ax,
    alpha=0.4,
)
ax.set_xlabel("UMAP 1", size=12, weight="bold")
ax.set_ylabel("UMAP 2", size=12, weight="bold")

# umap by nanopore abundance
## read processed illuminaAd
illuminaProcessedAd = sc.read_h5ad("processedIlluminaAd.h5ad") # the annotated h5ad

## get nanopore umap and transfer lable

adata = sc.read_10x_mtx("generateMtx/NanoporeMultiMat")
adata = adata[:, ~adata.var.index.str.contains("_")]
tempTools.plotCellScatter(adata)
adata = adata[adata.obs.index.isin(illuminaProcessedAd.obs.index)]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.012, max_mean=3, min_disp=1.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out(adata, ["n_count", "percent_ct"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
sc.pl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata, random_state=0)
adata.obs["illuminaLouvain"] = illuminaProcessedAd.obs["louvain"]

fig, ax = plt.subplots(figsize=(5.8, 3.5)) # figure 3d
sc.pl.umap(
    adata,
    color="illuminaLouvain",
    title="",
    show=False,
    legend_fontsize=10,
    legend_loc="on data",
    ax=ax,
)
ax.set_xlabel("UMAP 1", size=12, weight="bold")
ax.set_ylabel("UMAP 2", size=12, weight="bold")


# compare with sicelore
SICELORE_FOUND = "addBarcodeUmi_umifound_.bam"
siceloreFound = pysam.AlignmentFile(SICELORE_FOUND)
siceloreList = []
siceloreBcDict = {}
for read in siceloreFound:
    name = read.qname.split("_")[0]
    siceloreList.append(name)
    BC = read.get_tag("BC") + "_" + read.get_tag("U8")
    siceloreBcDict[name] = BC
siceloreBcDt = pd.DataFrame.from_dict(siceloreBcDict, "index")
siceloreBcDt.columns = ["siceloreBc"]

snuupyFoundPath = "parseMismatchResult/nanoporeReadWithBarcode.feather"
snuupyFound = pd.read_feather(snuupyFoundPath, columns=[0, 1])
snuupyFound.set_index("name", inplace=True)
snuupyFound.columns = ["snuupyBc"]

bothBcDt = pd.concat([siceloreBcDt, snuupyFound], axis=1)
bothBcDt.dropna(inplace=True)
bothDetectedNum, bothDetectedSameNum = len(bothBcDt), len(
    bothBcDt.query("siceloreBc == snuupyBc")
)
siceloreBcDetectedNum, snuupyDetectedNum = len(siceloreBcDt), len(planCFound)

from matplotlib_venn import venn2

v = venn2(
    subsets=(
        snuupyDetectedNum - bothDetectedNum,
        siceloreBcDetectedNum - bothDetectedNum,
        bothDetectedNum,
    ),
    set_labels=("snuupy", "Sicelore"),
    set_colors=("#EE7785", "#84B1ED"),
    alpha=1.0,
)  # figure 3b
[x.set_text(f"{int(x.get_text()):,d}") for x in v.subset_labels]


# multilayer clustering
INPUT_MTX = "generateMtx/IlluminaMultiMat"
adata = sc.read_10x_mtx(INPUT_MTX)
tempTools.plotCellScatter(adata)  # plot basic information
adata = tempTools.filterAd(adata, 5000, 550, 0.01)
geneAdata = adata[:, adata.var.index.str.find("_") == -1]
adata.X = adata.X.A
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=1.5)
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
sc.pl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
sc.pp.neighbors(adata, n_pcs=35)
sc.tl.umap(adata, random_state=0)
sc.tl.louvain(
    adata,
    resolution=1.2,
)
adata.obs["louvain"] = adata.obs["louvain"].astype(int) + 1
adata.obs["louvain"] = adata.obs["louvain"].astype("str").astype("category")
orderLouvain = [str(x) for x in range(1, len(adata.obs.louvain.unique()) + 1)]
adata.obs.louvain.cat.reorder_categories(orderLouvain, inplace=True)
setAdataColor = tempTools.setAdataColor
getAdataColor = tempTools.getAdataColor
typeColor = {
    "3": "#ffb482",
    "None": "#c0c0c0",
    "12": "#ff6600",
    "13": "#55a868",
    "7": "#8de5a1",
}
setAdataColor(
    adata,
    "louvain",
    {
        "1": "#1f77b4",
        "2": "#ffb442",
        "3": "#ffb482",
        "4": "#279e68",
        "5": "#aa40fc",
        "6": "#8c564b",
        "7": "#1E90FF",
        "8": "#b5bd61",
        "9": "#17becf",
        "10": "#aec7e8",
        "11": "#ffbb78",
        "12": "#ff6600",
        "13": "#a1c9f4",
        "14": "#c5b0d5",
        "15": "#c49c94",
        "16": "#f7b6d2",
    },
)
ax = sc.pl.umap(
    adata,
    color="louvain",
    title="Hybrid",
    show=False,
    legend_fontsize=10,
    groups=["3", "12", "13", "7"],
)  # figure 4a
ax.set_xlabel("UMAP 1", size=12, weight="bold")
ax.set_ylabel("UMAP 2", size=12, weight="bold")

subLouvainMap = {"3": "2.1", "12": "2.2", "7": "10.2", "13": "10.1", "None": "None"}
adata.obs["louvain"] = adata.obs["louvain"].map(subLouvainMap)
subAdata = adata[adata.obs.louvain.isin(["2.1", "2.2"])].copy()
rawNanoporeAdata = sc.read_10x_mtx("generateMtx/NanoporeMultiMat")

ax = sc.pl.umap(
    subAdata, color=["AT3G19010"], title="", show=False, cmap="Reds", size=200
)  # figure 4d
plt.ylim(1, 12)
plt.xlim(1, 12)

ax = sc.pl.umap(
    subAdata,
    color=["AT3G19010_False_fullySpliced"],
    title="",
    show=False,
    cmap="Reds",
    size=200,
)  # figure 4d
plt.ylim(1, 12)
plt.xlim(1, 12)

ax = sc.pl.umap(
    subAdata,
    color=["AT3G19010_True_fullySpliced"],
    title="",
    show=False,
    cmap="Reds",
    size=200,
)  # figure 4d
plt.ylim(1, 12)
plt.xlim(1, 12)

subAdata = subAdata[:, ~subAdata.var.index.str.contains("_")]
sc.tl.rank_genes_groups(
    subAdata, "louvain", method="wilcoxon", use_raw=False, groups=["2.2"]
)
sc.pl.rank_genes_groups(subAdata, n_genes=25, sharey=False, ax=ax)  # figure 4e
sc.pl.umap(
    subAdata,
    color="AT2G34600",
    show=False,
    cmap="Reds",
    size=200,
)  # figure 4e
