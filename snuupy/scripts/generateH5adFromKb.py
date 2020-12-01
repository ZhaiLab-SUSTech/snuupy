import os
import sh
import pyranges as pr
import pandas as pd
from io import StringIO
from collections import defaultdict
from more_itertools import chunked
from concurrent.futures import ThreadPoolExecutor
from loguru import logger
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import click
from . import kbParseTools


def getAdataFromKbNucleiResult(
    t2gPath, ecPath, splicePath, unsplicePath, adataPath, needUmiMappingInfo=False
):
    """
    get adata from kbpython result
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
    spliceBusDf = kbParseTools.getBustoolsMappingResult(
        t2gPath, ecPath, splicePath, "inner", True
    )
    logger.info("start parse unsplice bus")
    unspliceBusDf = kbParseTools.getBustoolsMappingResult(
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

    kbAd = kbParseTools.creatAnndataFromDf(
        kbMtxDf, spliced=kbSplicedMtxDf, unspliced=kbUnsplicedMtxDf
    )

    if adataPath:
        logger.info("start write anndata")
        kbAd.write(adataPath)

    if needUmiMappingInfo:
        return (kbAd, kbDf)
    else:
        return kbAd