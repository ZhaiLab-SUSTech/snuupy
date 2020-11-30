'''
Description: 
Author: Liuzj
Date: 2020-11-30 14:48:20
LastEditTime: 2020-11-30 15:54:21
LastEditors: Liuzj
'''
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


def __getSetOutersect(*setList):
    if len(setList) == 1:
        return setList[0]
    else:
        return setList[0].union(*setList[1:])


def writeWindowFasta(chromFastaDir, windowNum, fastaLs, i, totalCounts):
    fastaDir = f'{chromFastaDir}{windowNum}/'
    kbParseTools.mkdir(fastaDir)
    numLimit = 500
    for partialNum, partialFastaLs in enumerate((chunked(fastaLs, numLimit))):
        with open(f'{fastaDir}{partialNum}.fa', 'w') as fh:
            fh.write('\n'.join(partialFastaLs))
    if i % 1e4 == 0:
        logger.info(f"{i}/{totalCounts}  window were generated")

def generateIlluminaWindowFromKb(t2gPath, ecPath, splicePath, unsplicePath, gtfPath, illuminaWindowDir, windowSize):
    """
    generate illumina windows from kb_python results(workflow: nuclei)
    t2gPath: 
        index file
    ecPath: 
        matrix ec
    splicePath: 
        filtered spliced bus
    unsplicePath: 
        filtered spliced bus
    gtfPath: 
        gtf anno file, used to create kb ref
    illuminaWindowDir:
        dir stored illumina reads, end with '/'
    windowSize:
        windowSize
    """
    kbParseTools.mkdir(illuminaWindowDir)
    logger.info("start parse gff file")
    gtfDf = pr.read_gtf(gtfPath, as_df=True)
    gtfDf = gtfDf.query("Feature == 'exon'").reindex(['Chromosome', 'Start', 'End', 'gene_id'], axis=1)
    gtfDf = gtfDf.assign(Gene=lambda x: x["gene_id"]).groupby("Gene").agg(
        {"Chromosome": lambda x:x.iloc[0],
        "Start": 'min',
        "End": 'max'}
    )
    gtfDf = gtfDf.assign(
        StartWin=lambda x: x["Start"] // windowSize - 1,
        EndWin=lambda x: x["End"] // windowSize + 1,
    )
    gtfDf = gtfDf.to_dict('index')

    logger.info("start parse bus file")
    kbUmiSpliceMappingInfoDf = kbParseTools.getBustoolsMappingResult(t2gPath, ecPath, splicePath)
    kbUmiUnspliceMappingInfoDf = kbParseTools.getBustoolsMappingResult(t2gPath, ecPath, unsplicePath)
    kbUmiMappingInfoDf = pd.concat([kbUmiUnspliceMappingInfoDf,kbUmiSpliceMappingInfoDf])

    kbUmiMappingInfoDf = kbUmiMappingInfoDf.groupby('barcodeUmi').agg({'geneLs':lambda x:__getSetOutersect(*x)})
    kbUmiMappingInfoDf = kbUmiMappingInfoDf.assign(geneCounts = lambda df:df['geneLs'].map(len)).query("geneCounts >= 1")
    kbUmiMappingInfoDf = kbUmiMappingInfoDf.reset_index().assign(
        barcode=lambda df: df["barcodeUmi"].str.split("_").str[0],
        umi=lambda df: df["barcodeUmi"].str.split("_").str[1],
    ).assign(seq=lambda df: df["barcode"] + df["umi"])

    illuminaWindowContentDt = defaultdict(lambda: defaultdict(lambda :[]))
    for oneUmiNt in kbUmiMappingInfoDf.itertuples():
        for gene in oneUmiNt.geneLs:
            geneGtfDt = gtfDf[gene]
            geneChr = geneGtfDt['Chromosome']
            geneStartWin = geneGtfDt['StartWin']
            geneEndWin = geneGtfDt['EndWin']
            for singleWin in range(geneStartWin, geneEndWin+1):
                illuminaWindowContentDt[geneChr][singleWin].append(f'>{oneUmiNt.barcodeUmi}\n{oneUmiNt.seq}')


    i = 0
    totalCounts = sum([len(x) for x in illuminaWindowContentDt.values()])

    with ThreadPoolExecutor(24) as mtT:
        for chromNum, chromDt in illuminaWindowContentDt.items():
            chromFastaDir = f'{illuminaWindowDir}{chromNum}/'
            kbParseTools.mkdir(chromFastaDir)
            for windowNum, windowLs in chromDt.items():
                i += 1
                mtT.submit(writeWindowFasta, chromFastaDir, windowNum, windowLs, i, totalCounts)


