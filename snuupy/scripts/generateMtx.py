'''
Description: 
Author: Liuzj
Date: 2020-10-14 20:53:56
LastEditTime: 2020-11-24 13:23:42
LastEditors: Liuzj
'''
import pyranges as pr
import scanpy as sc
import pandas as pd
import numpy as np
import pysam
import portion
from collections import defaultdict
from loguru import logger
from .tools import transformExpressionMatrixTo10XMtx


def parseReadApaInfo(apaClusterPath, inBamPath, geneTag, expressionInfo):
    """
    解析apa信息 从wp pipeline中获得apacluster info 随后对每个read进行注释
    """
    apaCluster = pr.read_bed(apaClusterPath, True)
    apaCluster['Name'] = apaCluster['Name'] + '_APA'
    apaCluster['geneName'] = apaCluster['Name'].str.split('_').str[0]
    apaCluster = apaCluster.reindex(['geneName','Name','Start', 'End'], axis=1)
    apaClusterDict = defaultdict(lambda :{})
    for line in apaCluster.itertuples():
        apaClusterDict[line.geneName][line.Name] = portion.closedopen(line.Start, line.End)
    readsApaInfo = {}

    with pysam.AlignmentFile(inBamPath) as inBam:
        i = 0
        for read in inBam:
            i += 1
            readGene = read.get_tag(geneTag)
            geneApaInfo = apaClusterDict.get(readGene,'None')
            if geneApaInfo == 'None':
                readApaName = f'{readGene}_N_APA'
            else:
                if read.is_reverse:
                    readEndPos = read.positions[0]
                else:
                    readEndPos = read.positions[-1]

                readApaName = f'{readGene}_N_APA'
                for apaName, apaSpan in geneApaInfo.items():
                    if readEndPos in apaSpan:
                        readApaName = apaName
                        break
            readsApaInfo[read.qname] = readApaName
            if i % 100000 == 0:
                logger.info(f'{i} reads processed')

    readsApaInfo = pd.Series(readsApaInfo)
    useUmi = list(set(readsApaInfo.index) & set(expressionInfo.index))
    expressionInfo = pd.concat([expressionInfo, readsApaInfo])
    expressionInfo = expressionInfo.loc[expressionInfo.index.isin(useUmi)]
    return expressionInfo

def replaceNanoporeExpressionByIllumina(illuminaEx, nanoporeEx, nanoporeCorrectMtx):
    illuAdata = sc.read_10x_h5(illuminaEx, genome=None, gex_only=True)
    nanoAdata = sc.read_10x_mtx(nanoporeEx)
    illuEx = illuAdata.to_df()
    nanoEx = nanoAdata.to_df()
    nanoEx = nanoEx.loc[:,nanoEx.columns.str.find('_') != [-1]]
    nanoEx = nanoEx.join(illuEx, how='right')
    nanoEx.fillna(0, inplace=True)
    nanoEx.index = nanoEx.index.str.split('-').str[0]
    transformExpressionMatrixTo10XMtx(nanoEx, nanoporeCorrectMtx)




def generateMtx(apaClusterPath, inBamPath, geneTag, inIrInfoPath, outMtxDirPath, illuminaMtxDirPath, irMode, intronList, illuminaEx, onlyFullLength):
    """
    生成10X expression mtx 含有splice (APA) Info

    \b
    -a: polyAcluster bed ./single_cell_root.polya_cluster.bed
    --bam: 带有geneId的nanopore Bam
    --tag: geneId tag名 默认gi
    -i : ./step14_getIrInfo/irInfo.tsv
    -o : 文件夹 10X表达矩阵
    --intronList : 只使用该tsv内的intron
    --full-length
    """
    mode = []
    if (apaClusterPath != 'False') & (inBamPath  != 'False'):
        mode.append('apa')
        logger.warning('apa mode')
    if irMode:
        mode.append('ir')
        logger.warning('ir mode')
    if not mode:
        logger.warning('expression only mode')

    ## 读取ir情况 并且生成含有splice信息的表达矩阵
    readIrInfo = pd.read_csv(inIrInfoPath, sep='\t')

    # readIrInfo['fullySpliced'] = readIrInfo['IntronOverlapInfo'].isna().astype(str)
    # readIrInfo['fullySpliced'] = readIrInfo.pipe(lambda x:x['geneId'] + '_' + x['fullySpliced'] + '_fullySpliced')
    readIrInfo.set_index('Name', inplace=True)
    expressionGeneInfo = readIrInfo['geneId'].copy()
    if 'ir' in mode:


        readIrInfo['exonOverlapCounts'] = readIrInfo['ExonOverlapInfo'].str.split(',').map(lambda x:len(x))
        # readIrInfo['IntronOverlapInfo'].fillna('nan', inplace=True)
        # readIrInfo = readIrInfo.query('exonOverlapCounts != 1 or IntronOverlapInfo != "nan"')
        # readIrInfo.loc[irAllReadInfo['IntronOverlapInfo'] == 'nan', 'IntronOverlapInfo'] = np.nan

        if intronList:
            useIntron = pd.read_table(intronList)
            useIntron['geneId'] = useIntron['intron_id'].str.split('.').str[0]
            useIntron['intronId'] = useIntron['intron_id'].str.split('intron').str[-1].astype(int) - 1
            useIntronDict = defaultdict(lambda: [])
            for line in useIntron.itertuples():
                useIntronDict[line.geneId].append(line.intronId)
            useIntronDict = dict(useIntronDict)
        else:
            useIntronDict = False

        #0 去除 并且不计数
        #1 计数 不带intron
        #2 计数 带intron
        def getIntronSpliceInfo(line):
            if onlyFullLength:
                if line.exonOverlapCounts != line.GeneExonCounts:
                    return 'Ambiguous'
            if (line.exonOverlapCounts == 1) and pd.isna(line.IntronOverlapInfo):
                return 'Ambiguous'
            if not useIntronDict:
                if pd.isna(line.IntronOverlapInfo, ):
                    return 'True'
                else :
                    return 'False'
            else:
                if line.geneId not in useIntronDict.keys():
                    return 'Ambiguous'
                else:
                    if pd.isna(line.IntronOverlapInfo, ):
                        return 'True'
                    else:
                        retainedIntron = set(np.fromstring(line.IntronOverlapInfo, int, sep=','))
                        if retainedIntron & set(useIntronDict[line.geneId]):
                            return 'False'
                        else:
                            return 'True'
        readIrInfo['intronSpliceInfo'] = readIrInfo.apply(getIntronSpliceInfo, axis=1)
        readIrInfo = readIrInfo.query("intronSpliceInfo != 0")
        expressionIrInfo = readIrInfo['geneId'] + '_' + readIrInfo['intronSpliceInfo'] + '_fullySpliced'

        expressionInfo = pd.concat([expressionGeneInfo, expressionIrInfo])
    else:
        expressionInfo = expressionGeneInfo

    if 'apa' in mode:
        expressionInfo = parseReadApaInfo(apaClusterPath, inBamPath, geneTag, expressionInfo)


    expressionInfo = pd.DataFrame(expressionInfo)
    expressionInfo.reset_index(inplace=True)
    expressionInfo.columns = ['BcUmi', 'expressionInfo']
    expressionInfo['Bc'] = expressionInfo['BcUmi'].str.split('_').str[0]
    expressionInfo = expressionInfo.groupby('Bc')['expressionInfo'].apply(pd.value_counts)
    expressionInfo = expressionInfo.unstack().fillna(0).astype(int)

    transformExpressionMatrixTo10XMtx(expressionInfo, outMtxDirPath)

    replaceNanoporeExpressionByIllumina(illuminaEx, outMtxDirPath, illuminaMtxDirPath)
