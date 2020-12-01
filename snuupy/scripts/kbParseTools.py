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
        obs=pd.DataFrame(index = df.index),
        var=pd.DataFrame(index = df.columns)
    )
    for layerName, layerMtx in layerInfoDt.items():
        transformedAd.layers[layerName] = layerMtx
        
    return transformedAd
    

def mkdir(dirPath):
    try:
        sh.mkdir(dirPath)
    except:
        logger.warning(f'{dirPath} existed!!')

        
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
            trs = [int(x) for x in l[1].split(',')]
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


def getBustoolsMappingResult(t2gPath, ecPath, busPath, method='inner', filterUmi=False):
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
    ec2GeneFun = parseMatEc(
        ecPath,
        t2gDt, trLs)

    busFh = StringIO()
    
    logger.info("start parse bus")
    sh.bustools.text('-p', busPath, _out = busFh)

    busFh = StringIO(busFh.getvalue())

    busDf = pd.read_csv(
        busFh, sep="\t", header=None, names=["barcode", "umi", "ec", "count"]
    )
    
    logger.info("start get mapped gene")
    busDf = busDf.assign(geneLs = lambda x:x['ec'].map(ec2GeneFun))

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
    
    setFc = {'outer':__getSetOutersect,
            'inner':__getSetIntersect}[method]
    
    logger.info("start get finnal results")

    busDf = (
        busDf.groupby(["barcode", "umi"])["geneLs"]
        .agg(lambda x: setFc(*x))
        .reset_index()
    ).assign(
        barcodeUmi = lambda df:df['barcode'] + '_' + df['umi']
    )
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