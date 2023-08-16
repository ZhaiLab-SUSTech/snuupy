"""
Description: 
Author: Liuzj
Date: 2020-10-14 20:53:56
LastEditTime: 2020-11-24 13:23:42
LastEditors: Liuzj
"""
import pyranges as pr
import scanpy as sc
import pandas as pd
import numpy as np
import pysam
import portion
from collections import defaultdict
from loguru import logger
from collections import defaultdict
from tqdm import tqdm
import scipy.sparse as ss
import muon as mu


def parseReadApaInfo(apaClusterPath, inBamPath, geneTag, expressionInfo):
    """
    parse APA information, and classify each read into corresponding PAC
    """
    apaCluster = pr.read_bed(apaClusterPath, True)
    apaCluster["Chromosome"] = apaCluster["Chromosome"].astype(str)
    apaCluster["Name"] = apaCluster["Name"] + "_APA"
    apaCluster["geneName"] = apaCluster["Name"].str.split("_").str[:-2].str.join("_")
    apaCluster = apaCluster.reindex(["geneName", "Name", "Start", "End"], axis=1)
    apaClusterDict = defaultdict(lambda: {})
    for line in apaCluster.itertuples():
        apaClusterDict[line.geneName][line.Name] = portion.closedopen(
            line.Start, line.End
        )
    readsApaInfo = {}

    with pysam.AlignmentFile(inBamPath) as inBam:
        i = 0
        for read in inBam:
            i += 1
            if not read.has_tag(geneTag):
                continue
            readGene = read.get_tag(geneTag)
            geneApaInfo = apaClusterDict.get(readGene, None)
            if geneApaInfo is None:
                readApaName = f"{readGene}_N_APA"
            else:
                if read.is_reverse:
                    readEndPos = read.positions[0]
                else:
                    readEndPos = read.positions[-1]

                readApaName = f"{readGene}_N_APA"
                for apaName, apaSpan in geneApaInfo.items():
                    if readEndPos in apaSpan:
                        readApaName = apaName
                        break
            readsApaInfo[read.qname] = readApaName
            if i % 100000 == 0:
                logger.info(f"{i} reads processed")

    readsApaInfo = pd.Series(readsApaInfo)
    useUmi = list(set(readsApaInfo.index) & set(expressionInfo.index))
    expressionInfo = pd.concat([expressionInfo, readsApaInfo])
    expressionInfo = expressionInfo.loc[expressionInfo.index.isin(useUmi)]
    return expressionInfo


# def replaceNanoporeExpressionByIllumina(illuminaEx, nanoporeEx, nanoporeCorrectMtx):
#     if illuminaEx.endswith(".h5"):
#         illuAdata = sc.read_10x_h5(illuminaEx, genome=None, gex_only=True)
#     elif illuminaEx.endswith(".h5ad"):
#         illuAdata = sc.read_h5ad(illuminaEx)
#     else:
#         logger.error("Ambiguous Format !!!")
#         0 / 0
#     nanoAdata = sc.read_10x_mtx(nanoporeEx)
#     illuEx = illuAdata.to_df()
#     nanoEx = nanoAdata.to_df()
#     nanoEx = nanoEx.loc[:, nanoEx.columns.str.find("_") != -1]
#     nanoEx = nanoEx.join(illuEx, how="right")
#     nanoEx.fillna(0, inplace=True)
#     nanoEx.index = nanoEx.index.str.split("-").str[0]
#     transformExpressionMatrixTo10XMtx(nanoEx, nanoporeCorrectMtx)


def generateMtx(
    apaClusterPath,
    inBamPath,
    geneTag,
    inIrInfoPath,
    irMode,
    intronList,
    illuminaEx,
    onlyFullLength,
    h5muPath=None,
    outMtxDirPath=None,
    illuminaMtxDirPath=None,
):
    """
    generate expression mtx (format h5mu)
    """
    if intronList == 'False':  # compatible with click
        intronList = False

    if not h5muPath:
        h5muPath = illuminaMtxDirPath.rstrip("/") + ".h5mu"
        logger.warning(f"argument `illuminaMtxDirPath` is deprecated, output file will be redirected to `{h5muPath}`")

    mode = []
    if (apaClusterPath != "False") & (inBamPath != "False"):
        mode.append("apa")
        logger.warning("apa mode")
    if irMode:
        mode.append("ir")
        logger.warning("ir mode")
        if onlyFullLength:
            logger.warning("Only FullLength Mode")
    if not mode:
        logger.warning("expression only mode")

    readIrInfo = pd.read_csv(inIrInfoPath, sep="\t")

    # readIrInfo['fullySpliced'] = readIrInfo['IntronOverlapInfo'].isna().astype(str)
    # readIrInfo['fullySpliced'] = readIrInfo.pipe(lambda x:x['geneId'] + '_' + x['fullySpliced'] + '_fullySpliced')
    readIrInfo.set_index("Name", inplace=True)
    expressionGeneInfo = readIrInfo["geneId"].copy()
    if "ir" in mode:
        readIrInfo["ExonOverlapInfo"] = readIrInfo["ExonOverlapInfo"].fillna('') # some reads only mapped to intron, but not exon
        readIrInfo["exonOverlapCounts"] = (
            readIrInfo["ExonOverlapInfo"].str.split(",").map(lambda x: len(x))
        )
        # readIrInfo['IntronOverlapInfo'].fillna('nan', inplace=True)
        # readIrInfo = readIrInfo.query('exonOverlapCounts != 1 or IntronOverlapInfo != "nan"')
        # readIrInfo.loc[irAllReadInfo['IntronOverlapInfo'] == 'nan', 'IntronOverlapInfo'] = np.nan

        if intronList:
            useIntron = pd.read_table(intronList)
            useIntron["geneId"] = useIntron["intron_id"].str.split("\|").str[-1].str.split('_intron').str[0]
            useIntron["intronId"] = (
                useIntron["intron_id"].str.split("intron").str[-1].astype(int) - 1
            )
            useIntronDict = defaultdict(lambda: [])
            for line in useIntron.itertuples():
                useIntronDict[line.geneId].append(line.intronId)
            useIntronDict = dict(useIntronDict)
        else:
            useIntronDict = False

        def getIntronSpliceInfo(line):
            if onlyFullLength:
                if line.exonOverlapCounts != line.GeneExonCounts:
                    return "Ambiguous"
            if (line.exonOverlapCounts == 1) and pd.isna(line.IntronOverlapInfo):
                return "Ambiguous"
            if not useIntronDict:
                if pd.isna(
                    line.IntronOverlapInfo,
                ):
                    return "True"
                else:
                    return "False"
            else:
                if line.geneId not in useIntronDict.keys():
                    return "Ambiguous"
                else:
                    if pd.isna(
                        line.IntronOverlapInfo,
                    ):
                        return "True"
                    else:
                        retainedIntron = set(
                            np.fromstring(line.IntronOverlapInfo, int, sep=",")
                        )
                        if retainedIntron & set(useIntronDict[line.geneId]):
                            return "False"
                        else:
                            return "True"

        readIrInfo["intronSpliceInfo"] = readIrInfo.apply(getIntronSpliceInfo, axis=1)
        readIrInfo = readIrInfo.query("intronSpliceInfo != 0")
        expressionIrInfo = (
            readIrInfo["geneId"]
            + "_"
            + readIrInfo["intronSpliceInfo"]
            + "_fullySpliced"
        )

        expressionInfo = pd.concat([expressionGeneInfo, expressionIrInfo])
    else:
        expressionInfo = expressionGeneInfo

    if "apa" in mode:
        expressionInfo = parseReadApaInfo(
            apaClusterPath, inBamPath, geneTag, expressionInfo
        )

    expressionInfo = pd.DataFrame(expressionInfo)
    expressionInfo.reset_index(inplace=True)
    expressionInfo.columns = ["BcUmi", "expressionInfo"]
    expressionInfo["Bc"] = expressionInfo["BcUmi"].str.split("_").str[0]

    ls_bc = expressionInfo["Bc"].unique().tolist()
    ls_feature = expressionInfo["expressionInfo"].unique().tolist()

    dt_bc = {x: i for i, x in enumerate(ls_bc)}
    dt_feature = {x: i for i, x in enumerate(ls_feature)}

    ls_mtx = [[0 for x in ls_feature] for y in tqdm(ls_bc, "generate matrix content")]
    for line in tqdm(expressionInfo.itertuples(), total=len(expressionInfo)):
        ls_mtx[dt_bc[line.Bc]][dt_feature[line.expressionInfo]] += 1
    logger.info("transfer mtx format to sparse.matrix")
    ss_mtx = ss.csr_matrix(ls_mtx)
    ad_mtx = sc.AnnData(
        ss_mtx, obs=pd.DataFrame(index=ls_bc), var=pd.DataFrame(index=ls_feature)
    )
    ad_mtx.obs.index = ad_mtx.obs.index + "-1"
    ad_mtx = ad_mtx[sorted(ad_mtx.obs.index)]

    logger.info("generate mudata matrix")
    ad_illumina = sc.read_10x_h5(illuminaEx, genome=None, gex_only=True)
    ad_illumina.var_names_make_unique()
    ad_nanopore = ad_mtx[:, ~ad_mtx.var.index.str.contains(r"APA|fullySpliced")]
    ad_apa = ad_mtx[:, ad_mtx.var.index.str.contains(r"APA")]
    ad_splice = ad_mtx[:, ad_mtx.var.index.str.contains(r"fullySpliced")]
    dt_allAd = dict(
        illuminaAbu=ad_illumina, nanoporeAbu=ad_nanopore, apa=ad_apa, splice=ad_splice
    )
    if ad_apa.shape[1] == 0:
        del dt_allAd["ad_apa"]
    if ad_splice.shape[1] == 0:
        del dt_allAd["ad_splice"]
    md = mu.MuData(dt_allAd)
    md.write_h5mu(h5muPath)

    # expressionInfo = expressionInfo.groupby("Bc")["expressionInfo"].apply(
    #     pd.value_counts
    # )
    # expressionInfo = expressionInfo.unstack().fillna(0).astype(int)

    # transformExpressionMatrixTo10XMtx(expressionInfo, outMtxDirPath)

    # replaceNanoporeExpressionByIllumina(illuminaEx, outMtxDirPath, illuminaMtxDirPath)
