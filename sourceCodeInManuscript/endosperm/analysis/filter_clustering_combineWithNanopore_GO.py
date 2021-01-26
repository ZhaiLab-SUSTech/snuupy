import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import copy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
import tempTools

# read marker
markerDf = pd.read_excel(
    "PC2016-00845-CR1_Supplemental_Data_Set_1.xlsx", header=None
).rename(
    {0: "stage", 1: "tissue", 2: "locus"}, axis=1
)  # obtained from doi/10.1105/tpc.16.00845
nameMappingDt = {
    x: y
    for x, y in zip(
        ["EP", "SUS", "MCE", "PEN", "CZE", "CZSC", "GSC"],
        [
            "Embryo proper",
            "Suspensor",
            "Micropylar endosperm",
            "Peripheral endosperm",
            "Chalazal endosperm",
            "Chalazal seed coat",
            "General seed coat",
        ],
    )
}
markerDf = markerDf.assign(tissue=lambda df: df["tissue"].map(nameMappingDt))
htMarkerDt = (
    markerDf.query("stage == 'Heart'").groupby("tissue")["locus"].apply(list).to_dict()
)
orderLs = [
    "Embryo proper",
    "Suspensor",
    "Micropylar endosperm",
    "Peripheral endosperm",
    "Chalazal endosperm",
    "Chalazal seed coat",
    "General seed coat",
]
markerDt = (
    markerDf.query('stage in ["Globular", "Heart", "Linear Cotyledon"]')
    .groupby("tissue")["locus"]
    .agg(set)
    .to_dict()
)

# process illumina data
INPUT_H5 = "filtered_feature_bc_matrix.h5"
adata = sc.read_10x_h5(INPUT_H5, genome=None, gex_only=True)
tempTools.plotCellScatter(adata)
adata = adata[adata.obs["n_genes"] < 3000, :]
adata = adata[adata.obs["n_genes"] > 400, :]
adata.obs["n_count"] = adata.X.sum(axis=1)

## combine with nanopore data
snAd = sc.read_10x_mtx("generateMtx/NanoporeMultiMat")[adata.obs.index]
snCountAd = snAd[:, ~snAd.var.index.str.contains("_")]
tempTools.plotCellScatter(snCountAd)
illuminaWithSn_ObsDf = adata.obs.merge(
    snCountAd.obs,
    left_index=True,
    right_index=True,
    suffixes=("_illumina", "_Nanopore"),
)

fig, ax = plt.subplots(figsize=(3.5, 3.5))  # figure S7a
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
sns.scatterplot(
    "n_count_Nanopore",
    "n_count_illumina",
    data=illuminaWithSn_ObsDf,
    alpha=0.05,
    color="#000000",
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.set_ylabel("Illumina", size=12, fontweight="bold")
ax.set_xlabel("Nanopore", size=12, fontweight="bold")
ax.set_title("UMI counts per nucleus", fontweight="bold", size=12)

fig, ax = plt.subplots(figsize=(3.5, 3.5))  # figure S7a
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
sns.scatterplot(
    "n_genes_Nanopore",
    "n_genes_illumina",
    data=illuminaWithSn_ObsDf,
    alpha=0.05,
    color="#000000",
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
ax.yaxis.set_major_locator(ticker.MultipleLocator(500))
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.set_ylabel("Illumina", size=12, fontweight="bold")
ax.set_xlabel("Nanopore", size=12, fontweight="bold")
ax.set_title("Gene counts per nucleus", fontweight="bold", size=12)

illuminaWithSnAd = tempTools.combineAdataUseScanorama(
    [adata, snCountAd], "batch", ["Illumina", "Nanopore"]
)

fig, ax = plt.subplots(figsize=(6.4, 3.5))  #  figure S7b
sc.pl.umap(
    illuminaWithSnAd,
    color="batch",
    title="",
    show=False,
    legend_fontsize=10,
    ax=ax,
    alpha=0.5,
)
ax.set_xlabel("UMAP 1", size=12, weight="bold")
ax.set_ylabel("UMAP 2", size=12, weight="bold")

# clustering
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=1.5)
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
sc.tl.leiden(
    adata,
    resolution=0.8,
)
adata.obs["leiden"] = adata.obs["leiden"].astype(int) + 1
adata.obs["leiden"] = adata.obs["leiden"].astype(str).astype("category")
orderLouvain = [str(x) for x in range(1, len(adata.obs.leiden.unique()) + 1)]
adata.obs.leiden.cat.reorder_categories(orderLouvain, inplace=True)
adata.obs = adata.obs.assign(n_count_log=lambda df: np.log10(df["n_count"] + 1))

fig, ax = plt.subplots(figsize=(8, 3.5))  # Figure 5a
sc.pl.umap(
    adata,
    color="leiden",
    title="",
    show=False,
    legend_fontsize=14,
    legend_loc="on data",
    ax=ax,
)
ax.set_xlabel("UMAP 1", size=12, weight="bold")
ax.set_ylabel("UMAP 2", size=12, weight="bold")

# annotation
detectedGeneSet = set(adata.var.index)
orderTissueLs = [
    "Embryo proper",
    "Suspensor",
    "Micropylar endosperm",
    "Peripheral endosperm",
    "Chalazal endosperm",
    "Chalazal seed coat",
    "General seed coat",
]
markerDf["tissue"] = (
    markerDf["tissue"].astype("category").cat.reorder_categories(orderTissueLs)
)
heartStageMarkerDf = markerDf.query('locus in @detectedGeneSet and stage in ["Heart"]')
heartStageMarkerDt = heartStageMarkerDf.pipe(
    lambda df: df.groupby("tissue")["locus"].agg(list).to_dict()
)
heartCellScoreByGenesDf, heartClusterTypeRatio = tempTools.cellTypeAnnoByCellScore(
    adata, heartStageMarkerDt, "leiden"
)
adata.obs = adata.obs.assign(
    htClusterAnnoByCS=lambda df: df["leiden"].map(
        heartClusterTypeRatio.idxmax(1).to_dict()
    ),
    htCellAnnoByCS=lambda df: df.index.map(
        heartCellScoreByGenesDf["typeName"].to_dict()
    ),
)
heartClusterTypeRatio = heartClusterTypeRatio.reindex(
    [*orderTissueLs, "Unknown"], axis=1
)
fig, ax = plt.subplots(figsize=[3.5, 3.15])  # figure 5d
sns.heatmap(
    heartClusterTypeRatio.T.fillna(0),
    cmap="Reds",
    annot=heartClusterTypeRatio.T.fillna(0)
    .mul(adata.obs["leiden"].value_counts(), axis=1)
    .astype(int),
    cbar=False,
    fmt="d",
    ax=ax,
)
cb = ax.figure.colorbar(
    ax.collections[0],
)  # 显示colorbar
cb.ax.tick_params(labelsize=12)  # 设置colorbar刻度字体大小。
[x.set_fontsize(12) for x in ax.texts]
plt.ylabel("")
plt.xlabel("")

plt.yticks(fontsize=12, rotation=0)
plt.xticks(fontsize=12)

# calculate Incompletely read ratio
nanoporeAd = sc.read_10x_mtx("generateMtx/NanoporeMultiMat")[adata.obs.index]
nanoporeSpliceInfoAd = tempTools.getSpliceInfoFromSnuupyAd(nanoporeAd)
adata.obs = adata.obs.assign(
    irRatio=nanoporeSpliceInfoAd.to_df("unspliced").sum(1)
    / nanoporeSpliceInfoAd.to_df().sum(1),
)

fig, ax = plt.subplots(figsize=[6.29, 3.15])  # figure 5b
sc.pl.umap(
    adata, color="irRatio", title="", show=False, legend_fontsize=10, cmap="Reds", ax=ax
)
ax.set_xlabel("UMAP 1", size=12, weight="bold")
ax.set_ylabel("UMAP 2", size=12, weight="bold")

fig, ax = plt.subplots(figsize=[4, 3])  # figure 5c
plotDf = adata.obs.copy()
sns.boxplot(data=plotDf, x="leiden", y="irRatio", ax=ax, fliersize=0, width=0.8)
for i, box in enumerate(ax.artists):
    box.set_edgecolor(sns.color_palette()[i])
    box.set_facecolor("white")
    for j in range(6 * i, 6 * (i + 1)):
        ax.lines[j].set_color(sns.color_palette()[i])
        if j % 6 == 4:
            ax.lines[j].set_linewidth(3)
            ax.lines[j].set_xdata(
                np.array(ax.lines[j].get_xdata()) + np.array([0.05, -0.05])
            )

plt.xticks(fontsize=12)
plt.xlim(-1, 6)
sns.despine(top=True, right=True)
plt.xlabel("")

plt.yticks(fontsize=12)
plt.ylabel("Incompletely spliced ratio", fontsize=14)

# identify up-regulated gene
endospermAd = adata[adata.obs.query("leiden in ['1','2', '3', '4','5','6']").index]
sc.tl.rank_genes_groups(endospermAd, "leiden", method="wilcoxon", use_raw=True)
sc.pl.rank_genes_groups(endospermAd, n_genes=25, sharey=False, ax=ax)
sc.tl.filter_rank_genes_groups(
    endospermAd, max_out_group_fraction=0.25, min_fold_change=1.5
)
tempTools.getEnrichedGensDf(
    endospermAd, "4"
)  # then gene within this returned dataframe will be used fo GO analysis performed on website AGRIGO2

# GO analysis visualization
import glob
import re


def readGoResult(path):
    goDf = pd.read_table(path).assign(cluster=re.search(r"cluster(\w+).txt", path)[1])
    return goDf


goDf = pd.concat([readGoResult(x) for x in glob.glob("AGRIGO2_results*.txt")])
termTypeDt = dict(
    P="Biological Process", C="Cellular Component", F="Molecular Function"
)
goDf = goDf.query("FDR <= 0.05").assign(
    logTransformedFDR=lambda df: -np.log10(df["FDR"]),
    GO_Term=lambda df: df["GO_acc"] + ": " + df["Term"],
    term_type=lambda df: df["term_type"].map(termTypeDt),
)
goDf["-log(FDR)"] = goDf["logTransformedFDR"]
goDf.sort_values(["term_type", "cluster"], inplace=True)
plotDf = goDf.query("cluster == '4' and term_type == 'Cellular Component'")

fig, ax = plt.subplots(figsize=[4, 3.5])  # figure 5e
sns.barplot(
    data=plotDf,
    x="-log(FDR)",
    y="Term",
    ax=ax,
    dodge=False,
    color=sns.color_palette()[0],
)

for rec in ax.patches:
    rec.set_xy((rec.get_xy()[0], rec.get_xy()[1] + 0.2))
    rec.set_height(0.4)
plt.legend(loc=10, bbox_to_anchor=[-0.35, -0.1])

plt.xlabel("-Log$_{10}$(FDR)", fontsize=12)
plt.ylabel("")
sns.despine(top=True, right=True)