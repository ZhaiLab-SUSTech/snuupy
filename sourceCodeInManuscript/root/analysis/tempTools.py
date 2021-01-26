import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
from scipy.stats import spearmanr, pearsonr, zscore
from loguru import logger
from io import StringIO
from concurrent.futures import ProcessPoolExecutor as Mtp
import sh
import h5py

#################
##一些简单的工具###
#################


def mkdir(dirPath):
    try:
        sh.mkdir(dirPath)
    except:
        logger.warning(f"{dirPath} existed!!")


def normalizeMultiAd(multiAd, removeAmbiguous=True):
    """
    对二代三代数据分开normalize, 最终获得的每个细胞有3e4的reads
    """
    multiCountAd = multiAd[:, ~multiAd.var.index.str.contains("_")]
    multiOtherAd = multiAd[:, multiAd.var.index.str.contains("_")]
    sc.pp.normalize_total(multiCountAd, target_sum=1e4)
    sc.pp.normalize_total(multiOtherAd, target_sum=2e4)
    multiAd = sc.concat([multiCountAd, multiOtherAd], axis=1)
    if removeAmbiguous:
        multiAd = multiAd[
            :,
            ~(
                multiAd.var.index.str.contains("Ambiguous")
                | multiAd.var.index.str.contains("_N_")
            ),
        ]
    return multiAd


def filterAd(adata, geneMax, geneMin, ctMax):
    adata = adata[adata.obs["n_genes"] < geneMax, :]
    adata = adata[adata.obs["n_genes"] > geneMin, :]
    adata = adata[adata.obs["percent_ct"] < ctMax, :]
    return adata


def getAdataColor(adata, label):
    """
    获得adata中label的颜色
    """
    return {
        x: y
        for x, y in zip(adata.obs[label].cat.categories, adata.uns[f"{label}_colors"])
    }


def setAdataColor(adata, label, colorDt, hex=True):
    """
    设置adata中label的颜色
    """
    adata.obs[label] = adata.obs[label].astype("category")
    if not hex:
        from matplotlib.colors import to_hex

        colorDt = {x: hex(y) for x, y in colorDt.items()}
    adata.uns[f"{label}_colors"] = [colorDt[x] for x in adata.obs[label].cat.categories]
    return adata


def onlyKeepAdCount(adata):
    """
    去除adata和adata.raw中含有`_`的index
    """
    adata = adata.copy()
    adata = adata[:, ~adata.var.index.str.contains("_")]
    adataRaw = adata.raw.to_adata()
    adataRaw = adataRaw[:, ~adataRaw.var.index.str.contains("_")]
    adata.raw = adataRaw
    return adata


def plotCellScatter(adata):
    """
    绘制anndata基本情况
    """
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=0)

    adata.obs["n_count"] = adata.X.sum(axis=1)
    ctGene = (adata.var_names.str.startswith("ATCG")) | (
        adata.var_names.str.startswith("ATMG")
    )
    adata.obs["percent_ct"] = np.sum(adata[:, ctGene].X, axis=1) / np.sum(
        adata.X, axis=1
    )
    sc.pl.violin(
        adata, ["n_count", "n_genes", "percent_ct"], multi_panel=True, jitter=0.4
    )
    # fig, axs = plt.subplots(1, 3, figsize=(15, 3))
    # for x, y in zip(["n_count", "n_genes", "percent_ct"], axs):
    #     sc.pl.violin(adata, x, jitter=0.4, ax=y, show=False)
    #     plt.subplot(y)
    #     plt.ylabel("")


def creatAnndataFromDf(df, **layerInfoDt):
    """
    dataframe转换成anndata
    df,
    layerInfoDt:
        key为layer名
        value为mtx
    均行为barcode 列为feature 维度相同
    """
    transformedAd = anndata.AnnData(
        X=df.values,
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns),
    )
    for layerName, layerMtx in layerInfoDt.items():
        transformedAd.layers[layerName] = layerMtx

    return transformedAd


def mergeAdataExpress(
    adata, groupby="louvain", useRaw=True, logTransformed=True, method="mean"
):
    """
    通过adata.obs中的<groupby>对表达量合并
    """

    def __creatAnndataFromDf(df):
        transformedAd = anndata.AnnData(
            X=df.values,
            obs=pd.DataFrame(index=df.index),
            var=pd.DataFrame(index=df.columns),
        )
        return transformedAd

    adataExpressDf = (
        pd.DataFrame(adata.raw.X.A, columns=adata.raw.var.index, index=adata.obs.index)
        if useRaw
        else adata.to_df()
    )
    adataExpressDf = np.exp(adataExpressDf) - 1 if logTransformed else adataExpressDf

    groupbyExpressDf = (
        adataExpressDf.join(adata.obs[groupby]).groupby(groupby).agg(method)
    )
    groupbyExpressDf = (
        np.log(groupbyExpressDf + 1) if logTransformed else groupbyExpressDf
    )

    groupbyExpressAd = __creatAnndataFromDf(groupbyExpressDf)

    return groupbyExpressAd


def mergeAdata(adata, groupby, mergeLayer=[], method="sum"):
    """
    通过adata.obs中的<groupby>合并X和layer
    """
    adataXDf = adata.to_df()
    groupbyXDf = adataXDf.join(adata.obs[groupby]).groupby(groupby).agg(method)

    adataLayerDfDt = {}
    for singleLayer in mergeLayer:
        adataLayerDfDt[singleLayer] = (
            adata.to_df(singleLayer)
            .join(adata.obs[groupby])
            .groupby(groupby)
            .agg(method)
        )
    return creatAnndataFromDf(groupbyXDf, **adataLayerDfDt)


def getEnrichedGensDf(adata, name, useUns="rank_genes_groups_filtered"):
    """
    从adata的uns中获得富集基因
    """
    adataGeneInfoDt = adata.uns[useUns]
    useAttr = ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]
    nameUseInfoLs = [adataGeneInfoDt[x][name] for x in useAttr]
    nameEnrichedDf = pd.DataFrame(nameUseInfoLs, index=useAttr).T.dropna()
    return nameEnrichedDf.query("pvals_adj <= 0.05")


def calculateExpressionRatio(adata, clusterby):
    """
    逐个计算adata中每个基因在每个cluster中的表达比例

    adata:
        需要含有raw
    clusterby:
        adata.obs中的某个列名
    """
    transformAdataRawToAd = lambda adata: anndata.AnnData(
        X=adata.raw.X, obs=adata.obs, var=adata.raw.var
    )
    rawAd = transformAdataRawToAd(adata)
    expressionOrNotdf = (rawAd.to_df() > 0).astype(int)
    expressionOrNotdf[clusterby] = rawAd.obs[clusterby]
    expressionRatioDf = expressionOrNotdf.groupby(clusterby).agg(
        "sum"
    ) / expressionOrNotdf.groupby(clusterby).agg("count")
    return expressionRatioDf


def calculateGeneAverageEx(expressionMtxDf, geneDt, method="mean"):
    """
    根据geneDt对expressionMtxDf计算平均值或中位数

    expressionMtxDf:
        形如adata.to_df()

    geneDt:
        形如:{
    "type1": [
        "AT5G42235",
        "AT4G00540",
        ],
    "type2": [
        "AT1G55650",
        "AT5G45980",
        ],
    }
    method:
        'mean|median'

    """
    averageExLs = []
    for typeName, geneLs in geneDt.items():
        typeAvgExpress = (
            expressionMtxDf.reindex(geneLs, axis=1).mean(1)
            if method == "mean"
            else expressionMtxDf.reindex(geneLs, axis=1).median(1)
        )
        typeAvgExpress.name = typeName
        averageExLs.append(typeAvgExpress)
    averageExDf = pd.concat(averageExLs, axis=1)

    return averageExDf


def getClusterEnrichedGene(
    adata,
    useAttri=["names", "pvals", "pvals_adj", "logfoldchanges"],
    geneAnno=False,
    geneMarker=False,
):
    """
    获得每个cluster enriched的基因
    adata:
        obs中有louvain
        执行过:
            sc.tl.rank_genes_groups(adata)
            sc.tl.filter_rank_genes_groups(adata)

    geneAnno：
        dict;
        存在两个key['curator_summary', 'Note'] 从gff文件中转换而来

    geneMarker:
        dict;
        key为gene
        value为tissue


    useAttri = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
    """
    useLouvain = adata.obs["louvain"].unique()
    useAttriMap = {x: y for x, y in zip(range(len(useAttri)), useAttri)}
    useDict = adata.uns["rank_genes_groups_filtered"]
    allLouvainUseAttri = []
    for pointLouvain in useLouvain:
        pointLouvainUseAttri = []
        for pointAttr in useAttri:
            pointLouvainUseAttri.append(pd.Series(useDict[pointAttr][pointLouvain]))
        pointLouvainUseAttri = (
            pd.concat(pointLouvainUseAttri, axis=1)
            .rename(columns=useAttriMap)
            .assign(clusters=pointLouvain)
        )
        allLouvainUseAttri.append(pointLouvainUseAttri)
    allLouvainUseAttri = (
        pd.concat(allLouvainUseAttri)
        .dropna()
        .sort_values("clusters")
        .reindex(["clusters", "names", "pvals", "pvals_adj", "logfoldchanges"], axis=1)
    )
    if geneMarker:
        allLouvainUseAttri = allLouvainUseAttri.assign(
            markerInfo=allLouvainUseAttri["names"].map(
                lambda x: geneMarker.get(x, "Unknown")
            )
        )
    if geneAnno:
        allLouvainUseAttri = allLouvainUseAttri.assign(
            curator_summary=allLouvainUseAttri["names"].map(
                geneAnno["curator_summary"]
            ),
            note=allLouvainUseAttri["names"].map(geneAnno["Note"]),
        ).reset_index(drop=True)
    #     allLouvainUseAttri.dropna(inplace=True)
    return allLouvainUseAttri.query("pvals_adj > 0.05")


def __shuffleLabel(adata, label, i):
    """
    used for getEnrichedScore
    """
    shuffleAd = adata.copy()
    #     shuffleAd.obs[label] = adata.obs[label].sample(frac=1).reset_index(drop=True)
    shuffleAd.obs[label] = adata.obs[label].sample(frac=1, random_state=i).values
    #     print(shuffleAd.obs.iloc[0])
    shuffleClusterDf = (
        mergeAdataExpress(shuffleAd, label).to_df().reset_index().assign(label=i)
    )

    return shuffleClusterDf


def getEnrichedScore(adata, label, geneLs, threads=12, times=100):
    """
    获得ES值。ES值是通过对adata.obs中的label进行重排times次，然后计算原始label的zscore获得

    adata:
        必须有raw且为log-transformed

    label:
        adata.obs中的列名

    geneLs:
        需要计算的基因

    threads:
        使用核心数

    times:
        重排的次数
    """

    geneLs = geneLs[:]
    geneLs[0:0] = [label]
    adata = adata.copy()

    allShuffleClusterExpressLs = []
    with Mtp(threads) as mtp:
        for time in range(1, times + 1):
            allShuffleClusterExpressLs.append(
                mtp.submit(__shuffleLabel, adata, label, time)
            )

    allShuffleClusterExpressLs = [x.result() for x in allShuffleClusterExpressLs]
    originalClusterDf = (
        mergeAdataExpress(adata, label).to_df().reset_index().assign(label=0)
    )
    allShuffleClusterExpressLs.append(originalClusterDf)
    allShuffleClusterExpressDf = (
        pd.concat(allShuffleClusterExpressLs).set_index("label").reindex(geneLs, axis=1)
    )
    logger.info(f"start calculate z score")
    allShuffleClusterZscoreDf = (
        allShuffleClusterExpressDf.groupby(label)
        .apply(lambda x: x.set_index(label, append=True).apply(zscore))
        .reset_index(level=0, drop=True)
    )
    #     import pdb;pdb.set_trace()
    clusterZscoreDf = (
        allShuffleClusterZscoreDf.query(f"label == 0")
        .reset_index(level=0, drop=True)
        .fillna(0)
    )
    return clusterZscoreDf


#####解析IR结果####


def getSpliceInfoOnIntronLevel(irInfoPath, useIntronPath=None):
    """
    从intron水平获得剪接情况
    irInfoPath:
        snuupy getSpliceInfo的结果
    useIntronPath:
        使用的intron列表，需要表头'intron_id'

    return:
        adata:
            X: unsplice + splice
            layer[unspliced, spliced]
    """
    irInfoDf = pd.read_table(irInfoPath)
    intronCountMtxDt = {}
    intronRetenMtxDt = {}
    # 输入 0base
    # 输出 1base
    allLinesCounts = len(irInfoDf)
    for i, line in enumerate(irInfoDf.itertuples()):
        barcode = line.Name.split("_")[0]
        lineCountMtxDt = intronCountMtxDt.get(barcode, {})
        lineRetenMtxDt = intronRetenMtxDt.get(barcode, {})

        exonOverlapInfo = [int(x) for x in line.ExonOverlapInfo.split(",")]
        minIntron = min(exonOverlapInfo)
        maxIntron = max(exonOverlapInfo)
        intronCov = list(range(minIntron, maxIntron))

        if pd.isna(line.IntronOverlapInfo):
            intronOverlapInfo = []
        else:
            intronOverlapInfo = [int(x) for x in line.IntronOverlapInfo.split(",")]

        intronCov.extend(intronOverlapInfo)
        intronCov = set(intronCov)

        for intronCovNum in intronCov:
            lineCountMtxDt[f"{line.geneId}_intron_{intronCovNum+1}"] = (
                lineCountMtxDt.get(f"{line.geneId}_intron_{intronCovNum+1}", 0) + 1
            )
        for intronRentNum in intronOverlapInfo:
            lineRetenMtxDt[f"{line.geneId}_intron_{intronRentNum+1}"] = (
                lineRetenMtxDt.get(f"{line.geneId}_intron_{intronRentNum+1}", 0) + 1
            )

        intronCountMtxDt[barcode] = lineCountMtxDt
        intronRetenMtxDt[barcode] = lineRetenMtxDt
        if i % 1e5 == 0:
            logger.info(f"{i}/{allLinesCounts}")
    intronCountMtxDf = pd.DataFrame.from_dict(intronCountMtxDt, "index")
    intronRetenMtxDf = pd.DataFrame.from_dict(intronRetenMtxDt, "index")
    if useIntronPath:
        useIntronDf = pd.read_table(useIntronPath)
        useIntronLs = list(
            useIntronDf["intron_id"].str.split(".").str[0]
            + "_intron_"
            + useIntronDf["intron_id"].str.split("intron").str[1]
        )
        intronRetenMtxDf = intronRetenMtxDf.loc[
            :, intronRetenMtxDf.columns.isin(useIntronLs)
        ]
        intronCountMtxDf = intronCountMtxDf.loc[
            :, intronCountMtxDf.columns.isin(useIntronLs)
        ]
    intronCountMtxDf.index = intronCountMtxDf.index + "-1"
    intronRetenMtxDf.index = intronRetenMtxDf.index + "-1"
    intronRetenMtxDf = intronRetenMtxDf.fillna(0)
    intronCountMtxDf = intronCountMtxDf.fillna(0)
    intronCountMtxAd = creatAnndataFromDf(intronCountMtxDf)
    intronRetenMtxAd = creatAnndataFromDf(intronRetenMtxDf)

    useIntronLs = list(intronRetenMtxAd.var.index | intronCountMtxAd.var.index)
    useCellLs = list(intronRetenMtxAd.obs.index | intronCountMtxAd.obs.index)

    intronRetenMtxDf = (
        intronRetenMtxAd.to_df()
        .reindex(useIntronLs, axis=1)
        .reindex(useCellLs)
        .fillna(0)
    )
    intronCountMtxDf = (
        intronCountMtxAd.to_df()
        .reindex(useIntronLs, axis=1)
        .reindex(useCellLs)
        .fillna(0)
    )

    return creatAnndataFromDf(
        intronCountMtxDf,
        unspliced=intronRetenMtxDf,
        spliced=intronCountMtxDf - intronRetenMtxDf,
    )


def getSpliceInfoFromSnuupyAd(nanoporeAd):
    """
    用于从snuupy crMode产生的NanoporeMtx中提取产生splice和unsplice的read

    return:
        adata:
            X: unsplice + splice
            layer[unspliced, spliced]
    """
    nanoporeCountAd = nanoporeAd[:, ~nanoporeAd.var.index.str.contains("_")]
    unsplicedAd = nanoporeAd[:, nanoporeAd.var.index.str.contains("False_fullySpliced")]
    unsplicedAd.var.index = unsplicedAd.var.index.str.split("_").str[0]
    splicedAd = nanoporeAd[:, nanoporeAd.var.index.str.contains("True_fullySpliced")]
    splicedAd.var.index = splicedAd.var.index.str.split("_").str[0]
    useGeneLs = sorted(list(set(splicedAd.var.index) | set(unsplicedAd.var.index)))
    unsplicedDf = unsplicedAd.to_df().reindex(useGeneLs, axis=1).fillna(0)
    splicedDf = splicedAd.to_df().reindex(useGeneLs, axis=1).fillna(0)
    allSpliceDf = splicedDf + unsplicedDf
    return creatAnndataFromDf(allSpliceDf, spliced=splicedDf, unspliced=unsplicedDf)


def getDiffSplicedIntron(
    snSpliceIntronInfoAd,
    groupby,
    minCount,
    minDiff=0.1,
    threads=24,
    winflatPath="/public/home/jiajb/soft/IRFinder/IRFinder-1.2.5/bin/util/winflat",
):
    from pandarallel import pandarallel
    from statsmodels.stats import multitest
    import os

    pandarallel.initialize(nb_workers=threads)

    def calcuPvalue(line):
        nonlocal winflatPath
        xUnsplice = line.iloc[2]
        yUnsplice = line.iloc[3]
        xTotal = line.iloc[0]
        yTotal = line.iloc[1]
        resultStr = (
            os.popen(
                f"{winflatPath} -xvalue {xUnsplice} -yvalue {yUnsplice} -diff {xTotal} {yTotal}"
            )
            .read()
            .strip()
        )
        if not resultStr:
            return 1.0
        resultFloat = [
            float(x)
            for x in [x.strip().split("=")[-1].strip() for x in resultStr.split("\n")]
        ][1]

        return resultFloat

    allClusterDiffDt = {}
    for singleCluster in snSpliceIntronInfoAd.obs[groupby].unique():
        snSpliceIntronInfoAd.obs = snSpliceIntronInfoAd.obs.assign(
            cate=lambda df: df[groupby].map(
                lambda x: {singleCluster: singleCluster}.get(x, f"non-{singleCluster}")
            )
        )
        clusterSpliceIntronInfoAd = mergeAdata(
            snSpliceIntronInfoAd, "cate", ["unspliced", "spliced"]
        )
        clusterSpliceIntronInfoAd = clusterSpliceIntronInfoAd[
            :, clusterSpliceIntronInfoAd.to_df().min(0) >= minCount
        ]

        clusterSpliceIntronInfoDf = pd.concat(
            [
                clusterSpliceIntronInfoAd.to_df("unspliced").T,
                clusterSpliceIntronInfoAd.to_df().T,
            ],
            axis=1,
        )
        clusterSpliceIntronInfoDf.columns = [
            "unspliced",
            "non-unspliced",
            "total",
            "non-total",
        ]

        clusterSpliceIntronInfoDf["pvalue"] = clusterSpliceIntronInfoDf.parallel_apply(
            calcuPvalue, axis=1
        )
        clusterSpliceIntronInfoDf["fdr"] = multitest.fdrcorrection(
            clusterSpliceIntronInfoDf["pvalue"], 0.05
        )[1]

        clusterSpliceIntronInfoDf = clusterSpliceIntronInfoDf.assign(
            diffRatio=lambda df: df["unspliced"] / df["total"]
            - df["non-unspliced"] / df["non-total"]
        )

        clusterSpliceIntronInfoDf = clusterSpliceIntronInfoDf.eval(
            f"significantDiff = (fdr <= 0.05) & (diffRatio >= {minDiff})"
        )
        allClusterDiffDt[singleCluster] = clusterSpliceIntronInfoDf
        logger.info(
            f"{singleCluster} processed; {len(clusterSpliceIntronInfoDf)} input; {clusterSpliceIntronInfoDf['significantDiff'].sum()} output"
        )
    return allClusterDiffDt


#################
## 细胞注释
#################


def cellTypeAnnoByCorr(
    adata,
    bulkExpressionDf,
    threads=1,
    method="pearsonr",
    reportFinalUseGeneCounts=False,
    geneCountsCutoff=0,
    logTransformed=True,
    returnR=False,
    keepZero=True,
    useRaw=True,
    reportCounts=50,
):
    """
    通过bulk数据的相关性鉴定细胞类型

    adata: log-transformed adata

    bulkExpressionDf:
                                                AT1G01010  AT1G01030  AT1G01040
        bending cotyledon : chalazal endosperm   3.018853   2.430005   8.284994
        bending cotyledon : chalazal seed coat   2.385562   2.364294   8.674318
        bending cotyledon : embryo proper        2.258559   2.249158   7.577717

    returnR:
        返回最大值还是所有结果的r

    method:
        'pearsonr' | 'spearmanr'

    geneCountsCutoff:
        CPM
    """

    def __getSpearmanRForCell(cellExpressionSr):
        nonlocal i, bulkExpressionDf, keepZero, reportCounts, threads, method, geneCountsCutoff, reportFinalUseGeneCounts
        if not keepZero:
            cellExpressionSr = cellExpressionSr.pipe(lambda x: x[x != 0])
        cellExpressionSr = cellExpressionSr.pipe(lambda x: x[x >= geneCountsCutoff])
        #         print(cellExpressionSr)
        useGeneLs = cellExpressionSr.index
        bulkExpressionDf = bulkExpressionDf.reindex(useGeneLs, axis=1).dropna(axis=1)
        useGeneLs = bulkExpressionDf.columns
        cellExpressionSr = cellExpressionSr.reindex(useGeneLs)

        if reportFinalUseGeneCounts:
            logger.info(len(cellExpressionSr))

        if len(cellExpressionSr) <= 1:
            return None

        if method == "spearmanr":
            bulkExpressionCorrDf = bulkExpressionDf.apply(
                lambda x: spearmanr(cellExpressionSr, x)[0], axis=1
            )
        elif method == "pearsonr":
            bulkExpressionCorrDf = bulkExpressionDf.apply(
                lambda x: pearsonr(cellExpressionSr, x)[0], axis=1
            )
        else:
            logger.error("Unrecognized method")
            1 / 0

        i += 1
        if i % reportCounts == 0:
            logger.info(f"{i * threads} / {cellCounts} processed")

        if not returnR:
            mostSimilarBulk = bulkExpressionCorrDf.idxmax()
            return mostSimilarBulk
        else:
            return bulkExpressionCorrDf

    i = 0
    adata = adata.copy()
    cellCounts = len(adata)
    geneCountsCutoff = np.log(geneCountsCutoff + 1)

    adataExpressDf = (
        pd.DataFrame(adata.raw.X.A, columns=adata.raw.var.index, index=adata.obs.index)
        if useRaw
        else adata.to_df()
    )
    adataExpressDf = np.exp(adataExpressDf) - 1 if logTransformed else adataExpressDf
    adataExpressDf = adataExpressDf.div(adataExpressDf.sum(1), axis=0) * 1000000
    adataExpressDf = np.log(adataExpressDf + 1) if logTransformed else adataExpressDf
    #     print(adataExpressDf)

    if threads == 1:
        cellAnnotatedType = adataExpressDf.apply(__getSpearmanRForCell, axis=1)
    else:
        from pandarallel import pandarallel

        pandarallel.initialize(nb_workers=threads)
        cellAnnotatedType = adataExpressDf.parallel_apply(__getSpearmanRForCell, axis=1)
    return cellAnnotatedType


def cellTypeAnnoByMarker(adata, allMarkerUse, label="louvain", method="mean"):
    """
    通过marker基因表达量鉴定细胞类型

    adata:
        adata.obs中有louvain
        通过adata.raw.var来判断哪些基因表达
        存在raw, log-transformed

    allMarkerUse:
         {
            'Columella root cap':
                ['AT4G27400','AT3G18250', 'AT5G20045']
         }

    method = mean|median

    return:
        df: markerExpressCount(not logtransformed), expressRatio
    """
    #     import pdb; pdb.set_trace()
    markerRevDt = {z: x for x, y in allMarkerUse.items() for z in y}
    rawCountMtx = np.exp(adata.raw.to_adata().to_df()) - 1
    rawCountMtxWithLabel = rawCountMtx.join(adata.obs[label])

    clusterExMtx = rawCountMtxWithLabel.groupby(label).agg(method)
    #     return clusterExMtx.T
    clusterExMtxTr = clusterExMtx.T
    clusterExMtxTr.columns = clusterExMtxTr.columns.astype("str")
    clusterExMtxTr = clusterExMtxTr.assign(
        cellType=lambda df: df.index.map(lambda x: markerRevDt.get(x, "Others"))
    )
    clusterTypeExpMtx = clusterExMtxTr.groupby("cellType").agg(method).T

    cellExRatioMtxTr = rawCountMtx.applymap(lambda x: 1 if x > 0 else 0).T
    cellExRatioMtxTr.columns = cellExRatioMtxTr.columns.astype("str")
    cellExRatioMtxTr = cellExRatioMtxTr.assign(
        cellType=lambda df: df.index.map(lambda x: markerRevDt.get(x, "Others"))
    )
    cellExRatioMtx = cellExRatioMtxTr.groupby("cellType").apply(lambda x: x.mean(0)).T
    cellExRatioMtxWithLabel = cellExRatioMtx.join(adata.obs[label])
    clusterExRatioMtx = cellExRatioMtxWithLabel.groupby(label).agg(method)

    finalMtx = (
        pd.concat([clusterTypeExpMtx.unstack(), clusterExRatioMtx.unstack()], axis=1)
        .reset_index()
        .rename({0: "express", 1: "ratio", "level_0": "cellType"}, axis=1)
    ).query("cellType != 'Others'")

    return finalMtx


def cellTypeAnnoByMarkerOld(
    adata, allMarkerUse, expressionMtx, zscoreby="cluster", method="mean"
):
    """
    通过marker基因表达量鉴定细胞类型

    adata:
        adata.obs中有louvain
        通过adata.raw.var来判断哪些基因表达

    allMarkerUse:
    {
     'Zhang et al.':
         {
            'Columella root cap':
                ['AT4G27400','AT3G18250', 'AT5G20045']
         }
    }
    expressionMtx:
         由adata.to_df()获得:
             没有log-transformed
             没有筛选基因
             经过normalize_sum

    zscoreby = cluster|cell

    method = mean|median
    """

    def _getMarkerExpressionGene(adata, allMarkerUse):
        """
        去除marker中不表达的基因
        adata 存在raw
        allMarkerUse
            {'Zhang et al.': {'Columella root cap': ['AT4G27400','AT3G18250', 'AT5G20045']}}
        """
        expressionGene = set(adata.raw.var.index)

        integrateMarkerGene = {}
        for x, y in allMarkerUse.items():
            singleMarkerGeneUse = {}
            for j, k in y.items():
                k = list(set(k) & expressionGene)
                singleMarkerGeneUse[j] = k
            integrateMarkerGene[x] = singleMarkerGeneUse
        return integrateMarkerGene

    expressionMtx = expressionMtx.copy(True)
    #     expressionMtx = np.log2(expressionMtx + 1)

    allMarkers = _getMarkerExpressionGene(adata, allMarkerUse)

    expressionMtx = expressionMtx.join(adata.obs["louvain"], how="inner")
    allLouvain = expressionMtx["louvain"].unique()
    expressionCounts = (
        expressionMtx.groupby("louvain")
        .apply(lambda x: x.drop("louvain", axis=1).pipe(lambda y: y.sum() / len(y)))
        .fillna(0)
    )
    expressionCounts = np.log2(expressionCounts + 1)
    expressionSizes = (
        expressionMtx.groupby("louvain")
        .apply(
            lambda x: x.drop("louvain", axis=1).pipe(lambda y: (y > 0).sum() / len(y))
        )
        .fillna(0)
    )
    if zscoreby == "cluster":
        expressionZscore = expressionCounts.apply(zscore)
    elif zscoreby == "cell":
        expressionMtx = np.log2(expressionMtx.drop("louvain", axis=1) + 1)
        expressionMtx = expressionMtx.apply(zscore)
        expressionZscore = (
            expressionMtx.join(adata.obs["louvain"], how="inner")
            .groupby("louvain")
            .apply(lambda x: x.drop("louvain", axis=1).pipe(lambda y: y.sum() / len(y)))
            .fillna(0)
        )
    #     expressionCounts = expressionMtx.groupby('louvain').apply(lambda x:x.drop('louvain', axis=1).pipe(lambda y: y.sum() / len(y))).fillna(0)
    #     expressionCounts = expressionMtx.groupby('louvain').apply(lambda x:x.drop('louvain', axis=1).pipe(lambda y: y.sum() / (y > 0).sum())).fillna(0)

    groupAllClustersExpressionCounts = []
    groupAllClustersExpressionZscore = []
    groupAllClustersExpressionSizes = []
    groupNames = []
    for stage, tissueGenes in allMarkers.items():
        for tissue, genes in tissueGenes.items():
            if method == "mean":
                groupGeneCountsDf = expressionCounts.loc[:, genes].mean(1)
                groupGeneZscoreDf = expressionZscore.loc[:, genes].mean(1)
                groupGeneSizesDf = expressionSizes.loc[:, genes].mean(1)
            elif method == "median":
                groupGeneCountsDf = expressionCounts.loc[:, genes].median(1)
                groupGeneZscoreDf = expressionZscore.loc[:, genes].median(1)
                groupGeneSizesDf = expressionSizes.loc[:, genes].median(1)
            groupGeneCountsDf.name = f"{stage} {tissue}"
            groupGeneZscoreDf.name = f"{stage} {tissue}"
            groupGeneSizesDf.name = f"{stage} {tissue}"
            groupAllClustersExpressionCounts.append(groupGeneCountsDf)
            groupAllClustersExpressionZscore.append(groupGeneZscoreDf)
            groupAllClustersExpressionSizes.append(groupGeneSizesDf)

    groupAllClustersExpressionCounts = pd.concat(groupAllClustersExpressionCounts, 1)
    groupAllClustersExpressionZscore = pd.concat(groupAllClustersExpressionZscore, 1)
    groupAllClustersExpressionSizes = pd.concat(groupAllClustersExpressionSizes, 1)
    groupAllClustersExpression = pd.concat(
        [
            groupAllClustersExpressionSizes.stack(),
            groupAllClustersExpressionZscore.stack(),
            groupAllClustersExpressionCounts.stack(),
        ],
        axis=1,
    )
    groupAllClustersExpression.reset_index(inplace=True)
    groupAllClustersExpression.columns = [
        "louvain",
        "tissues",
        "Percentage of expressed nuclei",
        "Z-score of Expression",
        "Average expression",
    ]
    #     groupAllClustersExpression = groupAllClustersExpression.reset_index()
    #     groupAllClustersExpression.columns = ['louvain','tissues','Percentage of expressed nuclei', 'Average expression']
    return groupAllClustersExpression


def cellTypeAnnoByClusterEnriched(
    arrayExpressDf_StageTissue,
    clusterEnrichedGeneDf,
    useCluster="all",
    useGeneCounts=10,
):
    """
    使用cluster enriched基因在bulk数据中的表达情况对cluster进行注释
    暂时只能用在胚乳数据上 待进一步优化

    arrayExpressDf_StageTissue: dataframe, 形如
                                             AT1G01010  AT1G01030  AT1G01040
    stage             correctedTissue
    bending cotyledon chalazal endosperm     3.018853   2.430005   8.284994
                      chalazal seed coat     2.385562   2.364294   8.674318
                      embryo proper          2.258559   2.249158   7.577717
                      general seed coat      2.000998   2.168115   7.721052
                      peripheral endosperm   2.503232   2.154924   8.002944

    clusterEnrichedGeneDf: getClusterEnrichedGene输出
    useCluster：'all'|['1', '2']
    useGeneCounts: 每个cluster使用的基因数
    """
    stageOrderLs = [
        "pre-globular",
        "globular",
        "heart",
        "linear-cotyledon",
        "bending cotyledon",
        "mature green",
    ]
    tissueOrderLs = [
        "chalazal endosperm",
        "micropylar endosperm",
        "peripheral endosperm",
        "chalazal seed coat",
        "general seed coat",
        "embryo proper",
        "suspensor",
    ]
    expressList = list(arrayExpressDf_StageTissue.columns)
    clusterEnrichedGene_FilteredDf = (
        clusterEnrichedGeneDf.sort_values(
            ["clusters", "logfoldchanges"], ascending=[True, False]
        )
        .groupby("clusters")
        .apply(lambda x: x.loc[x["names"].isin(expressList)].iloc[:useGeneCounts])
        .reset_index(drop=True)
    )

    clusterEnrichedGeneName_FilteredDf = clusterEnrichedGene_FilteredDf.groupby(
        "clusters"
    )["names"].agg(lambda x: list(x))

    clusterEnrichedGeneFc_FilteredDf = clusterEnrichedGene_FilteredDf.groupby(
        "clusters"
    )["logfoldchanges"].agg(lambda x: np.exp2(x).mean())

    print(clusterEnrichedGeneName_FilteredDf.map(lambda x: len(x)))

    print(clusterEnrichedGeneFc_FilteredDf)

    if useCluster == "all":
        useClusterLs = list(clusterEnrichedGeneName_FilteredDf.index)
    else:
        useClusterLs = useCluster

    #     return arrayExpressDf_StageTissue

    # clusterName = useClusterLs[0]
    #     import pdb;pdb.set_trace()
    for clusterName in useClusterLs:
        fig, ax = plt.subplots(figsize=[5, 3])
        clusterEnrichedGeneExpressPatternInBulkDf = (
            arrayExpressDf_StageTissue.loc[
                :, clusterEnrichedGeneName_FilteredDf[clusterName]
            ]
            .median(1)
            .unstack()
            .reindex(stageOrderLs)
            .reindex(tissueOrderLs, axis=1)
        )
        sns.heatmap(clusterEnrichedGeneExpressPatternInBulkDf, cmap="Reds", ax=ax)
        ax.set_title(f"Cluster {clusterName}")
        ax.set_xlabel("Tissue")
        ax.set_ylabel("Stage")

    EnrichedGeneExpressPatternInBulkDf = clusterEnrichedGeneName_FilteredDf.map(
        lambda x: arrayExpressDf_StageTissue.loc[:, x].median(1).idxmax()
    )
    return EnrichedGeneExpressPatternInBulkDf


def cellTypeAnnoByEnrichedScore(adata, label, markerGeneDt, threads=12, times=100):
    """
    通过enriched score对cluster进行注释

    adata:
        必须有raw且为log-transformed

    label:
        adata.obs中的列名

    markerGeneDt:
        形如:{
    "type1": [
        "AT5G42235",
        "AT4G00540",
        ],
    "type2": [
        "AT1G55650",
        "AT5G45980",
        ],
    }

    threads:
        使用核心数

    times:
        重排的次数
    """
    adata = adata.copy()
    adata = adata[:, ~adata.var.index.str.contains("_")]
    adataRaw = adata.raw.to_adata()
    adataRaw = adataRaw[:, ~adataRaw.var.index.str.contains("_")]
    adata.raw = adataRaw

    markerGeneLs = list(set([y for x in markerGeneDt.values() for y in x]))
    clusterEnrichedScoreDf = getEnrichedScore(
        adata, label, markerGeneLs, threads, times
    )
    clusterTypeAvgEnrichedScoreDf = calculateGeneAverageEx(
        clusterEnrichedScoreDf, markerGeneDt
    )
    return clusterTypeAvgEnrichedScoreDf


def cellTypeAnnoByCellScore(adata, markerGenesDt, clusterLabel):
    """
    利用cellscore计算每个细胞的type

    adata:
        anndata
    markerGenesDt:
        {type:[genes]}
    clusterLabel:
        cluster label

    return:
        cellScoreByGenesDf:
            每个细胞的cellScore
        clusterTypeRatio:
            每个cluster的type比例
    """
    adata = adata.copy()
    adata = adata[:, ~adata.var.index.str.contains("_")]
    adataRaw = adata.raw.to_adata()
    adataRaw = adataRaw[:, ~adataRaw.var.index.str.contains("_")]
    adata.raw = adataRaw
    
    for name, genes in markerGenesDt.items():
        sc.tl.score_genes(adata, genes, score_name=name, use_raw=True)

    cellScoreByGenesDf = adata.obs[markerGenesDt.keys()]
    cellScoreByGenesDf["maxType"], cellScoreByGenesDf["maxScore"] = (
        cellScoreByGenesDf.idxmax(1),
        cellScoreByGenesDf.max(1),
    )
    cellScoreByGenesDf["typeName"] = cellScoreByGenesDf["maxType"]
    cellScoreByGenesDf.loc[
        cellScoreByGenesDf.loc[:, "maxScore"] < 0, "typeName"
    ] = "Unknown"

    adata.obs["typeName"] = cellScoreByGenesDf["typeName"]

    clusterTypeRatio = (
        adata.obs.groupby(clusterLabel)["typeName"]
        .apply(lambda x: x.value_counts() / len(x))
        .unstack()
    )
    return cellScoreByGenesDf, clusterTypeRatio


#######
# 绘图
#######


def plotLabelPercentageInCluster(adata, groupby, label, labelColor):
    """
    根据label在adata.obs中groupby的占比绘图

    groupby:
        表明cluster。需要存在于adata.obs
    label:
        展示的占比。需要存在于adata.obs
    labelColor:
        label的颜色
    """
    groupbyWithLabelCountsDf = (
        adata.obs.groupby(groupby)[label].apply(lambda x: x.value_counts()).unstack()
    )
    groupbyWithLabelCounts_CumsumPercDf = groupbyWithLabelCountsDf.pipe(
        lambda x: x.cumsum(1).div(x.sum(1), 0) * 100
    )
    legendHandleLs = []
    legendLabelLs = []
    for singleLabel in groupbyWithLabelCounts_CumsumPercDf.columns[::-1]:
        ax = sns.barplot(
            x=groupbyWithLabelCounts_CumsumPercDf.index,
            y=groupbyWithLabelCounts_CumsumPercDf[singleLabel],
            color=labelColor[singleLabel],
        )
        legendHandleLs.append(
            plt.Rectangle((0, 0), 1, 1, fc=labelColor[singleLabel], edgecolor="none")
        )
        legendLabelLs.append(singleLabel)
    legendHandleLs, legendLabelLs = legendHandleLs[::-1], legendLabelLs[::-1]
    plt.legend(legendHandleLs, legendLabelLs, bbox_to_anchor=[1, 1], frameon=False)
    plt.xlabel(groupby.capitalize())
    plt.ylabel(f"Percentage")
    return ax


##########
##整合数据
##########
def combineAdataUseScanorama(
    adataLs, batchKey, batchCateLs, subSample=False, subSampleCounts=0
):
    """
    利用Scanorama整合不同adata
    adataLs:
        [adata1, adata2]
    batchKey:
        添加的label
    batchCateLs:
        每个batch的名称 需要和adataLs一致
    subSampleCounts:
        下采样的样本数。
    return:
        整合后的adata
    """
    import scanorama

    adataLs = [x.copy() for x in adataLs]
    if subSample:
        sampleSize = min([x.shape[0] for x in adataLs])
        if subSampleCounts:
            sampleSize = min(sampleSize, subSampleCounts)
        logger.info(f"sample size: {sampleSize}")
        [sc.pp.subsample(x, n_obs=sampleSize) for x in adataLs]

    for adata in adataLs:
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata, flavor="seurat", n_top_genes=2000, inplace=True
        )

    print(f"↓↓↓↓↓↓↓↓↓{batchCateLs}↓↓↓↓↓↓↓↓")
    combineScanoramaLs = scanorama.correct_scanpy(adataLs, return_dimred=True)
    print(f"↑↑↑↑↑↑↑↑↑{batchCateLs}↑↑↑↑↑↑↑↑")
    combineAdata = sc.concat(
        combineScanoramaLs, label=batchKey, index_unique="-", keys=batchCateLs
    )
    sc.pp.neighbors(combineAdata, n_pcs=50, use_rep="X_scanorama")
    sc.tl.umap(combineAdata)
    return combineAdata


def combineAdataUseScanoramaOld(
    adataLs, batchKey, batchCateLs, subSample=False, subSampleCounts=0
):
    """
    利用Scanorama整合不同adata
    adataLs:
        [adata1, adata2]
    batchKey:
        添加的label
    batchCateLs:
        每个batch的名称 需要和adataLs一致
    subSampleCounts:
        下采样的样本数。
    return:
        整合后的adata
    """
    import scanorama

    adataLs = [x.copy() for x in adataLs]
    if subSample:
        sampleSize = min([x.shape[0] for x in adataLs])
        if subSampleCounts:
            sampleSize = min(sampleSize, subSampleCounts)
        logger.info(f"sample size: {sampleSize}")
        [sc.pp.subsample(x, n_obs=sampleSize) for x in adataLs]

    combineAdata = adataLs[0].concatenate(
        adataLs[1:], batch_key=batchKey, batch_categories=batchCateLs
    )

    sc.pp.normalize_per_cell(combineAdata, counts_per_cell_after=1e4)
    sc.pp.log1p(combineAdata)

    sc.pp.highly_variable_genes(
        combineAdata, min_mean=0.0125, max_mean=3, min_disp=1.5, batch_key=batchKey
    )
    sc.pl.highly_variable_genes(combineAdata)

    varGenes = combineAdata.var.highly_variable

    varGenes = varGenes[varGenes].keys()

    varGenes = list(varGenes)

    alldata = {}

    for oneBatchName in combineAdata.obs[batchKey].unique():
        alldata[oneBatchName] = combineAdata[
            combineAdata.obs[batchKey] == oneBatchName, varGenes
        ]

    combineAdataLs = list(alldata.values())

    print(f"↓↓↓↓↓↓↓↓↓{batchCateLs}↓↓↓↓↓↓↓↓")
    combineScanoramaLs = scanorama.correct_scanpy(combineAdataLs, return_dimred=True)
    print(f"↑↑↑↑↑↑↑↑↑{batchCateLs}↑↑↑↑↑↑↑↑")
    combineAdata = sc.concat(combineScanoramaLs)
    #     import pdb; pdb.set_trace()
    #     combineScanoramaAr = np.concatenate(combineScanoramaLs)

    #     combineAdata.obsm["SC"] = combineScanoramaAr

    #     combineAdata.raw = combineAdata
    #     combineAdata = combineAdata[:, varGenes]
    #     sc.pp.scale(combineAdata, max_value=10)
    #     sc.tl.pca(combineAdata, svd_solver="arpack", n_comps=50)
    #     sc.pl.pca(combineAdata)
    sc.pp.neighbors(combineAdata, n_pcs=50, use_rep="X_scanorama")
    sc.tl.umap(combineAdata)
    return combineAdata


##########
##bustools
##########
def parseBustoolsIndex(t2gPath, t2cPath=False):
    """
    解析bustools的索引文件

    t2gPath:
        kbpython index t2g

    return:
        t2gDt: key为转录本名 value为基因名
        trLs: 所有的转录本 sorted
        geneLs: 所有的基因 sorted
    """
    t2cTrList = []
    if not t2cPath:
        pass
    else:
        with open(t2cPath) as t2cFh:
            for line in t2cFh:
                t2cTrList.append(line.strip())
    t2cTrSet = set(t2cTrList)

    t2gDt = {}
    trLs = []
    with open(t2gPath) as t2gFh:
        for line in t2gFh:
            lineLs = line.split()
            if (not t2cTrSet) | (lineLs[0] in t2cTrSet):
                t2gDt[lineLs[0]] = lineLs[1]
                trLs.append(lineLs[0])
            else:
                pass
    geneLs = list(set(t2gDt.values()))
    return t2gDt, trLs, geneLs


def parseMatEc(ecPath, t2gDt, trLs):
    """
    解析表达矩阵ec

    ecPath:
        kb-python 产生的bus文件
    t2gDt:
        parseBustoolsIndex 产生的字典 key为转录本名 value为基因名
    trLs:
        parseBustoolsIndex 产生的list 所有的转录本 sorted

    return:
        函数 输入ec 输出对应的geneLs
    """
    ecsDt = {}
    with open(ecPath) as ecFh:
        for line in ecFh:
            l = line.split()
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(",")]
            ecsDt[ec] = trs

    def __ec2g(ec):
        """
        函数 输入ec 输出对应的geneLs
        """
        if ec in ecsDt:
            return set(t2gDt[trLs[t]] for t in ecsDt[ec])
        else:
            return set([])

    return __ec2g


def getBustoolsMappingResult(t2gPath, ecPath, busPath, method="inner", filterUmi=False):
    """
    解析kbpython结果 并获得表达情况

    t2gPath:
        kbpython index t2g
    ecPath:
        kbpython matrix ec
    busPath:
        kbpython bus
    method:
        inner|outer  inner指与多个基因比对上取交集 outer指取并集
    filterUmi:
        是否只保留仅比对上一个基因的umi
    return:
        df:
            columns = ["barcode", "umi", "geneSet"]
    """
    logger.info("start parse index")
    t2gDt, trLs, geneLs = parseBustoolsIndex(t2gPath)

    logger.info("start parse mat")
    ec2GeneFun = parseMatEc(ecPath, t2gDt, trLs)

    busFh = StringIO()

    logger.info("start parse bus")
    sh.bustools.text("-p", busPath, _out=busFh)

    busFh = StringIO(busFh.getvalue())

    busDf = pd.read_csv(
        busFh, sep="\t", header=None, names=["barcode", "umi", "ec", "count"]
    )

    logger.info("start get mapped gene")
    busDf = busDf.assign(geneLs=lambda x: x["ec"].map(ec2GeneFun))

    def __getSetIntersect(*setList):
        if len(setList) == 1:
            return setList[0]
        else:
            return setList[0].intersection(*setList[1:])

    def __getSetOutersect(*setList):
        if len(setList) == 1:
            return setList[0]
        else:
            return setList[0].union(*setList[1:])

    setFc = {"outer": __getSetOutersect, "inner": __getSetIntersect}[method]

    logger.info("start get finnal results")

    busDf = (
        busDf.groupby(["barcode", "umi"])["geneLs"]
        .agg(lambda x: setFc(*x))
        .reset_index()
    ).assign(barcodeUmi=lambda df: df["barcode"] + "_" + df["umi"])
    logger.info(f"before filter {len(busDf)} Umis found")
    if filterUmi:
        logger.info("Umi filter: True, start filter")
        busDf = (
            busDf.assign(geneCounts=lambda df: df["geneLs"].map(len))
            .query("geneCounts == 1")
            .assign(geneLs=lambda df: df["geneLs"].map(lambda x: list(x)[0]))
        ).drop("geneCounts", axis=1)

        logger.info(f"after filter {len(busDf)} Umis found")
    else:
        logger.info("Umi filter: False")
    return busDf


def getAdataFromKbNucleiResult(
    t2gPath, ecPath, splicePath, unsplicePath, needUmiMappingInfo=False, adataPath=False
):
    """
    用于从kb的nuclei策略中获得adata
    t2gPath:
        kbpython index t2g
    ecPath:
        kbpython matrix ec
    splicePath:
        kbpython splice bus
    unsplicePath:
        kbpython unsplice bus
    needUmiMappingInfo:.
        need umi mapping info or not
    adataPath:
        adata store path

    return:
        adata
        (umiMappingDf)
    """
    logger.info("start parse splice bus")
    spliceBusDf = getBustoolsMappingResult(t2gPath, ecPath, splicePath, "inner", True)
    logger.info("start parse unsplice bus")
    unspliceBusDf = getBustoolsMappingResult(
        t2gPath, ecPath, unsplicePath, "inner", True
    )
    kbDf = pd.concat([spliceBusDf, unspliceBusDf])

    logger.info("start get overlap umi")
    kbUmiGeneCountsSr = kbDf.groupby("barcodeUmi")["geneLs"].agg("count")
    unoverlapUmiLs = list(kbUmiGeneCountsSr.pipe(lambda x: x[x == 1]).index)
    overlapUmiSt = set(kbUmiGeneCountsSr.pipe(lambda x: x[x != 1]).index)
    overlapUseSr = (
        kbDf.query("barcodeUmi in @overlapUmiSt")
        .groupby("barcodeUmi")["geneLs"]
        .apply(lambda df: True if df.iat[0] == df.iat[1] else False)
    )
    overlapUmiUseLs = list(overlapUseSr.pipe(lambda x: x[x]).index)
    useUmiLs = sorted([*unoverlapUmiLs, *overlapUmiUseLs])

    logger.info("start filter overlap umi and creat anndata")
    kbDf = kbDf.drop_duplicates("barcodeUmi").query("barcodeUmi in @useUmiLs")
    kbMtxDf = (
        kbDf.groupby("barcode")["geneLs"]
        .apply(lambda df: df.value_counts())
        .unstack()
        .fillna(0)
    )
    kbSplicedMtxDf = (
        spliceBusDf.query("barcodeUmi in @unoverlapUmiLs")
        .groupby("barcode")["geneLs"]
        .apply(lambda df: df.value_counts())
        .unstack()
        .reindex(kbMtxDf.index)
        .reindex(kbMtxDf.columns, axis=1)
        .fillna(0)
    )
    kbUnsplicedMtxDf = (
        unspliceBusDf.query("barcodeUmi in @unoverlapUmiLs")
        .groupby("barcode")["geneLs"]
        .apply(lambda df: df.value_counts())
        .unstack()
        .reindex(kbMtxDf.index)
        .reindex(kbMtxDf.columns, axis=1)
        .fillna(0)
    )

    kbAd = creatAnndataFromDf(
        kbMtxDf, spliced=kbSplicedMtxDf, unspliced=kbUnsplicedMtxDf
    )

    if adataPath:
        logger.info("start write anndata")
        kbAd.write(adataPath)

    if needUmiMappingInfo:
        return (kbAd, kbDf)
    else:
        return kbAd


######
# cellranger
#####
def extractReadCountsByUmiFromTenX(molInfoPath, kitVersion="v2", filtered=True):
    """
    用于从cellRanger文件中提取read counts

    molInfoPath:
        molecule_info.h5
    kitVersion:
        v2|v3。v2:umi 10bp, v3:umi 12bp
    filtered:
        只使用过滤后的

    return:
        dataframe: columns = ['barcodeUmi', 'featureName', 'readCount']

    """
    umiLength = {"v2": 10, "v3": 12}[kitVersion]

    def NumToSeq():
        nonlocal umiLength

        numToBase = {"00": "A", "01": "C", "10": "G", "11": "T"}

        def _numToSeq(num):
            num = int(num)
            numStr = f"{num:032b}"[-umiLength * 2 :]
            return "".join(
                [numToBase[numStr[2 * x : 2 * x + 2]] for x in range(umiLength)]
            )

        return _numToSeq

    numToSeq = NumToSeq()
    molInfo = h5py.File(molInfoPath, "r")
    allBarcode = molInfo["barcodes"][()].astype(str)
    allFeature = molInfo["features/id"][()].astype(str)

    allUmiCount = pd.DataFrame(
        np.c_[
            molInfo["umi"][()],
            allBarcode[molInfo["barcode_idx"][()]],
            molInfo["count"][()],
            allFeature[molInfo["feature_idx"][()]],
        ]
    )
    if filtered:
        allUmiCount = allUmiCount[
            allUmiCount[1].isin(
                allBarcode[molInfo["barcode_info/pass_filter"][()][:, 0]]
            )
        ]
    allUmiCount[0] = allUmiCount[0].map(numToSeq)
    allUmiCount["barcodeUmi"] = allUmiCount[1] + "_" + allUmiCount[0]
    allUmiCount = allUmiCount.reindex(["barcodeUmi", 3, 2], axis=1, copy=False)
    allUmiCount.columns = ["barcodeUmi", "featureName", "readCount"]
    return allUmiCount