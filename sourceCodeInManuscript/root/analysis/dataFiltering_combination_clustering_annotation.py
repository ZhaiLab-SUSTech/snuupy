import copy
import scanorama
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
import tempTools


# load snRNA-seq data
INPUT_H5 = "filtered_feature_bc_matrix.h5"
nucleiAdata = sc.read_10x_h5(INPUT_H5, genome=None, gex_only=True)
tempTools.plotCellScatter(nucleiAdata)  # plot basic information
nucleiAdata = tempTools.filterAd(nucleiAdata, 2300, 350, 0.01)

# combine and calculate Alignment score
# all need data can be obtained in the corresponding GEO database

## load mp data
mpAdata = sc.read_10x_h5("mp_filtered_feature_bc_matrix.h5")
tempTools.plotCellScatter(mpAdata)
mpAdata = mpAdata[mpAdata.obs.eval("500 <= n_genes <= 5000")]

## load science data
scienceAdata = sc.read_10x_h5("science_filtered_gene_bc_matrices.h5")
tempTools.plotCellScatter(scienceAdata)
scienceAdata = scienceAdata[:, ~scienceAdata.var.index.duplicated()]
scienceAdata.var.index = scienceAdata.var.gene_ids

## load dc data
dcAdata = sc.read_text("dc_Root_single_cell_wt_datamatrix.csv", ",")
dcAdata = dcAdata.T
tempTools.plotCellScatter(dcAdata)

## load pp data
ppAdata = sc.read_text("pp_5way_merge_raw.tsv", "\t")
ppAdata = ppAdata.T
ppuseBarcodeLs = list(
    ppAdata.obs.index.str.split("_").map(
        lambda x: x[2] + "-1-" + str(int(x[1][-1]) - 1)
    )
)
ppRawAdatas = [
    sc.read_10x_h5(x)
    for x in [
        "pp_1_filtered_feature_bc_matrix.h5",
        "pp_2_filtered_feature_bc_matrix.h5",
        "pp_3_filtered_feature_bc_matrix.h5",
    ]
]
ppRawAdata = ppRawAdatas[0].concatenate(ppRawAdatas[1:])
ppAdata = ppRawAdata[ppRawAdata.obs.index.isin(ppuseBarcodeLs)]

## load pc data
pcUseCellLs = list(
    pd.read_csv("pc_controlAnno.tsv", sep="\s").index.str.split("_").str[0]
)
pcAdata = sc.read_10x_h5("pc_filtered_feature_bc_matrix.h5")
pcAdata = pcAdata[pcAdata.obs.index.isin(pcUseCellLs)]

## start combine

allAdataLs = [nucleiAdata, mpAdata, scienceAdata, dcAdata, ppAdata, pcAdata]
allAdataNameLs = [
    "Nuclei",
    "Molecular Plant",
    "Science",
    "Developmental Cell",
    "Plant Physiology",
    "Plant Cell",
]
allCombineAdata = tempTools.combineAdataUseScanorama(
    allAdataLs, "batch", allAdataNameLs, True, 1186
)  # 1186 is the nuclei counts of the snRNA-seq data

## plot combined results
combineAsDf = pd.DataFrame(
    [
        [0.0, 0.538, 0.361, 0.198, 0.077, 0.17],
        [0.0, 0.0, 0.659, 0.358, 0.327, 0.213],
        [0.0, 0.0, 0.0, 0.862, 0.323, 0.547],
        [0.0, 0.0, 0.0, 0.0, 0.572, 0.808],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.367],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ],
)  # this matrix is obtained from the Standard output of the function combineAdataUseScanorama
combineAsDf = combineAsDf.T + combineAsDf
dfNameLs = [
    "Day 10, Nuclei\n(0.5cm primary root tips)",
    "Day 10, Zhang et al.\n(0.5cm primary root tips)",
    "Day 6, Wendrich et al.\n(Primary root tips)",
    "Day 6, Denyer et al.\n(1 cm primary root tips)",
    "Day 5, Ryu et al.\n(Primary root tips)",
    "Day 8, Jean-Baptiste et al.\n(Whole roots)",
]
np.fill_diagonal(combineAsDf.values, 1)
combineAsDf.index = dfNameLs
combineAsDf.columns = dfNameLs
orderedNameLs = [
    "Day 10, Nuclei\n(0.5cm primary root tips)",
    "Day 10, Zhang et al.\n(0.5cm primary root tips)",
    "Day 8, Jean-Baptiste et al.\n(Whole roots)",
    "Day 6, Denyer et al.\n(1 cm primary root tips)",
    "Day 6, Wendrich et al.\n(Primary root tips)",
    "Day 5, Ryu et al.\n(Primary root tips)",
]
combineAsDf = combineAsDf.reindex(orderedNameLs).reindex(orderedNameLs, axis=1)
combineAsDf = combineAsDf.iloc[::-1]

fig, ax = plt.subplots(figsize=[12, 10])  # figure 2a
sns.heatmap(combineAsDf, cmap="Reds", vmin=0, vmax=1.0, ax=ax, cbar=False)
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
cb = ax.figure.colorbar(
    ax.collections[0],
)
cb.ax.tick_params(labelsize=24)
plt.xticks(rotation=90, ma="right", fontsize=24)
plt.yticks(fontsize=24)

allAdataLs = [nucleiAdata, mpAdata]
allCombineAdata = tempTools.combineAdataUseScanorama(
    allAdataLs, "batch", allAdataNameLs, True, 1186
)  # 1186 is the nuclei counts of the snRNA-seq data

fig, ax = plt.subplots()  # figure 2b
ax = sc.pl.umap(allCombineAdata, color="batch", title="", vmax=1500, alpha=0.5, ax=ax)

# clustering
adata = nucleiAdata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=1.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
sc.pl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata, random_state=0)
sc.tl.louvain(
    adata,
    resolution=1.2,
)
adata.obs["louvain"] = adata.obs["louvain"].astype(int) + 1
adata.obs["louvain"] = adata.obs["louvain"].astype(str).astype("category")
orderLouvain = [str(x) for x in range(1, len(adata.obs.louvain.unique()) + 1)]
adata.obs.louvain.cat.reorder_categories(orderLouvain, inplace=True)

# annotation
## read marker gene
markerGenesDf = pd.read_table(
    "SpecificGeneFromBenfey.tsv"
)  # can be obtained from DOI: 10.1101/2020.06.29.178863; the processed table were listed in supplemental
correctTypeNameDt = {
    x: y
    for x, y in zip(
        [
            "Atrichoblast",
            "Columella",
            "Cortex",
            "Endodermis",
            "Lateral Root Cap",
            "Pericycle",
            "Phloem",
            "Procambium",
            "Quiescent Center",
            "Stem Cell Niche",
            "Trichoblast",
            "Xylem",
        ],
        [
            "Non-hair",
            "Columella",
            "Cortex",
            "Endodermis",
            "Lateral Root Cap",
            "Pericycle",
            "Phloem",
            "Procambium",
            "Quiescent Center",
            "Stem Cell Niche",
            "Root Hair",
            "Xylem",
        ],
    )
}
markerGenesDf = markerGenesDf.assign(
    celltype=lambda x: x["celltype"].map(correctTypeNameDt)
).assign(cluster=lambda x: x["dev stage"] + "-" + x["celltype"])
markerGenesDf = markerGenesDf.query("gene in @adata.raw.var.index")
markerGenesDf = markerGenesDf.sort_values(["celltype", "combined_rank"])
stageMarkerGenesDf = (
    markerGenesDf.groupby("cluster").apply(lambda x: x.iloc[:30]).reset_index(drop=True)
)
stageMarkerGenesDt = (
    stageMarkerGenesDf.groupby("cluster")["gene"].agg(lambda x: list(x)).to_dict()
)
stageMarkerGenesDt.pop("Meristematic-Quiescent Center")

markerGenesDf = (
    markerGenesDf.groupby(["celltype"])
    .apply(lambda x: x.drop_duplicates("gene").iloc[:30])
    .reset_index(drop=True)
)
markerGenesDt = (
    markerGenesDf.groupby("celltype")["gene"].agg(lambda x: list(set(x))).to_dict()
)
markerGenesDt.pop("Quiescent Center")

## start annotation
cellScoreByGenesDf, clusterTypeRatio = tempTools.cellTypeAnnoByCellScore(
    adata, markerGenesDt, "louvain"
)
adata.obs["typeName"] = cellScoreByGenesDf["typeName"]
typeColorsDt = {
    "Endodermis": "#8c8c8c",
    "Stem Cell Niche": "#55a868",
    "Phloem": "#8172b3",
    "Pericycle": "#8da0cb",
    "Xylem": "#c44e52",
    "Cortex": "#a1c9f4",
    "Non-hair": "#ff9f9b",
    "Root Hair": "#fab0e4",
    "Lateral Root Cap": "#dd8452",
    "Columella": "#da8bc3",
    "Procambium": "#f5e77a",
    "Unknown": "#cfcfcf",
}
tempTools.setAdataColor(adata, "typeName", typeColorsDt)
fig, ax = plt.subplots(figsize=[5, 3])  # figure S3a
sc.pl.umap(adata, color="typeName", show=False, legend_fontsize=10, ax=ax)
plt.tight_layout()

sns.clustermap(clusterTypeRatio.T.fillna(0), cmap="Reds", figsize=[6, 6])  # figure S3b

nonHairAdata = adata[
    adata.obs["louvain"].isin(["2", "6", "14"])
].copy()  # These cluster were annotated as Non hair
stageMarkerGenesDt = {
    x: stageMarkerGenesDt[x]
    for x in ["Meristematic-Non-hair", "Elongating-Non-hair", "Mature-Non-hair"]
}
stageCellScoreByGenesDf, stageClusterTypeRatio = tempTools.cellTypeAnnoByCellScore(
    nonHairAdata, stageMarkerGenesDt, "louvain"
)
nonHairAdata.obs["stageTypeName"] = stageCellScoreByGenesDf["typeName"]
typeColorsDt = {
    "Meristematic-Non-hair": "#ffb482",
    "Elongating-Non-hair": "#ff9f9b",
    "Mature-Non-hair": "#8de5a1",
    "Unknown": "#cfcfcf",
}
tempTools.setAdataColor(nonHairAdata, "stageTypeName", typeColorsDt)
fig, ax = plt.subplots(figsize=[5, 3])  # figure S3c
ax = sc.pl.umap(nonHairAdata, color="stageTypeName", show=False, ax=ax, size=100)
plt.tight_layout()

###plot annotated Umap
annoResultDt = clusterTypeRatio.idxmax(1).to_dict()
annoResultDt.update(stageClusterTypeRatio.idxmax(1).to_dict())
annoResultDt = {x: f"{x}: {y}" for x, y in annoResultDt.items()}
adata.obs = adata.obs.assign(
    typeAnnoByCellScore=lambda x: x["louvain"].map(annoResultDt)
)
typeColorDt =  {'1: Stem Cell Niche': '#937860',
 '2: Mature-Non-hair': '#ffb482',
 '3: Pericycle': '#8da0cb',
 '4: Stem Cell Niche': '#55a868',
 '5: Endodermis': '#8c8c8c',
 '6: Mature-Non-hair': '#8de5a1',
 '7: Root Hair': '#fab0e4',
 '8: Endodermis': '#4c72b0',
 '9: Lateral Root Cap': '#dd8452',
 '10: Cortex': '#a1c9f4',
 '11: Xylem': '#c44e52',
 '12: Stem Cell Niche': '#d0bbff',
 '13: Phloem': '#8172b3',
 '14: Elongating-Non-hair': '#ff9f9b'}
 clusterColorDt = {
    '7': '#fab0e4',
    '14': '#ff9f9b',
    '2': '#ffb482',
    '6': '#8de5a1',
    '8': '#4c72b0',
    '5': '#8c8c8c',
    '10': '#a1c9f4',
    '11': '#c44e52',
    '3': '#8da0cb',
    '13': '#8172b3',
    '9': '#dd8452',
    '1': '#937860',
    '4': '#55a868',
    '12': '#d0bbff'
}
tempTools.setAdataColor(adata, "typeAnnoByCellScore", typeColorDt)
ax = sc.pl.umap(
    adata, color="typeAnnoByCellScore", title="", color_map="Reds", show=False
)  # figure 1c

sc.tl.rank_genes_groups(adata, "louvain", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.tl.filter_rank_genes_groups(
    adata, min_fold_change=1.5, max_out_group_fraction=0.25, min_in_group_fraction=0.25
)

plotMarkerGene = {
    "Root Hair": ["AT1G48930"],
    "Non-hair": ["AT1G79840"],
    "Mature Non-hair": ["AT3G01970"],
    "Cortex": ["AT5G07990"],
    "Endodermis": ["AT5G42180"],
    "Xylem": ["AT5G12870"],
    "Phloem": ["AT3G43270"],
    "Pericycle": ["AT1G32450"],
    "Stem Cell Niche": ["AT4G25630"],
    "Root Cap": ["AT1G79580"],
}
plotMarkerGeneNameDt = {
    "AT1G48930": "GH9C1",
    "AT1G79840": "GL2",
    "AT3G01970": "WRKY45",
    "AT5G07990": "TT7",
    "AT5G42180": "PER64",
    "AT5G12870": "MYB46",
    "AT3G43270": "PME32",
    "AT1G32450": "NRT1.5",
    "AT4G25630": "FIB2",
    "AT1G79580": "SMB",
}
ax = sc.pl.stacked_violin(
    adata,
    plotMarkerGene,
    groupby="louvain",
    stripplot=True,
    jitter=0.25,
    color=[clusterColorDt[str(x)] for x in range(1, 15)],
    figsize=(12, 16),
    show=False,
)  # figure 1d
mainAx = ax["mainplot_ax"]
plt.axes(mainAx)
xticksPa = plt.xticks()
plt.xticks(
    xticksPa[0],
    [plotMarkerGeneNameDt[x.get_text()] for x in xticksPa[1]],
    rotation=-45,
    ha="left",
    fontsize=20,
    fontstyle="italic",
)
plt.yticks(fontsize=20)
plt.axes(mainAx)
geneAx = ax["gene_group_ax"]
plt.axes(geneAx)
for singleText in geneAx.texts:
    singleText.update(dict(rotation=45, ha="left", fontsize=20))