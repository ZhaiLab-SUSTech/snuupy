"""
@Date: 2020-06-17 12:27:25
LastEditors: liuzj
LastEditTime: 2020-08-06 11:39:23
@Description: file content
@Author: liuzj
FilePath: /liuzj/projects/singleCell/01_pipeline/pipeline_planC/01_scripts/step08_polishRead.py
"""
import pandas as pd
import numpy as np
import time
import click
import pyfastx
import copy
import os
import subprocess
from loguru import logger
from concurrent.futures import ProcessPoolExecutor
from more_itertools import chunked
from itertools import repeat
from jpy_tools.ReadProcess import readFasta


def GetConsensusSeq():
    seq2numDict = {
    '-':4,
    'A':0,
    'T':1,
    'G':2,
    'C':3
}
    num2seqDict = {i:j for j,i in seq2numDict.items()}
    
    def _getMaxCountElement(Array):
        return np.argmax(np.bincount(Array))
    
    def _getConsensusSeq(alignedFasta):
        seqList = []
        fastaAligned = readFasta(alignedFasta)
        for x in fastaAligned:
            seqList.append([seq2numDict[x] for x in list(x.seq)])
        seqArray = np.array(seqList)
        consensusSeq = np.apply_along_axis(_getMaxCountElement, 0, seqArray)
        consensusSeq = consensusSeq[consensusSeq != 4]
        consensusSeq = ''.join([num2seqDict[x] for x in consensusSeq])
        return consensusSeq
    
    return _getConsensusSeq


def getConsensesFasta(readLs, tempPath, penaltyPath, finalPath=False):
    getConsensusSeq = GetConsensusSeq()
    
    with open(f'{tempPath}_all.fa', 'w') as fh:
        for read in readLs[:10]:
            fh.write(f'>{read[0]}\n{read[1]}\n')
    
    os.system(f'poa -read_fasta {tempPath}_all.fa -pir {tempPath}_all.aln {penaltyPath}')
    
    consensusSeq = getConsensusSeq(f'{tempPath}_all.aln')
    
    readName = tempPath.split('/')[-1]
    if not finalPath:
        finalPath = f'{tempPath}_round0.fa'

    with open(finalPath, 'w') as fh:
        fh.write(f'>{readName}_{len(readLs)}\n{consensusSeq}\n')
    
#     os.system(f'rm {tempPath}*')

def getPolishRead(tempPath, finalPath, times=1):
    ##效果不佳 暂且放弃
    def _getPolishCommand(commandExecuted, targetFasta, useFasta, polishedFasta):
        # commandExecuted += f'bwa index -a is {targetFasta} && bwa mem {targetFasta} {useFasta} > {targetFasta}.sam && racon {useFasta} {targetFasta}.sam {targetFasta} > {polishedFasta} &&'
        commandExecuted += f'minimap2 --secondary=no -ax map-ont {targetFasta} {useFasta} > {targetFasta}.sam && racon {useFasta} {targetFasta}.sam {targetFasta} > {polishedFasta} &&'
        return commandExecuted.strip()
    
    
    allPolishedPath = [tempPath + f'_round{x}.fa' for x in range(times + 1)]
    commandExecuted = '('
    for x in range(times):
        commandExecuted = _getPolishCommand(commandExecuted, allPolishedPath[x], f'{tempPath}_all.fa', allPolishedPath[x+1])
        
    commandExecuted += f'mv {allPolishedPath[x+1]} {finalPath} && rm {tempPath}* )'
    commandExecuted = commandExecuted.strip()
    return commandExecuted

def polishSeq(barcodeWithReadLs, nanoporeDict, tempDirPath, finalDirPath, penaltyPath):
    barcodeWithReadLs = copy.deepcopy(barcodeWithReadLs)
    barcode, readLs = barcodeWithReadLs
    finalPath = f'{finalDirPath}{barcode}.fa'
    returnList = [barcode, readLs]
    lengthList = []
    for read in readLs:
        readName, readStrand = read[0].split('_')
        if readStrand == '1':
            readSeq = nanoporeDict[readName].antisense
            seqLength = len(readSeq)
        else:
            readSeq = nanoporeDict[readName].seq
            seqLength = len(readSeq)
        lengthList.append(seqLength)
        read.extend([readSeq])

    if len(readLs) <= 2:
#         print(barcode,lengthList)
        lengthMinIndex = lengthList.index(max(lengthList))
        referenceRead = readLs.pop(lengthMinIndex)
        read = referenceRead
        with open(finalPath, 'w') as fh:
            fh.write(f'>{barcode}_{len(readLs) + 1}\n{read[1]}\n')
        return '(echo "no need")'

    else:
        tempReadPath = f'{tempDirPath}{barcode}'
        
        getConsensesFasta(readLs, tempReadPath, penaltyPath)
        return getPolishRead(tempReadPath, finalPath)

def chunkPolishSeq(chunkBarcodeWithReadIter, nanoporeReadPath, tempDirPath, finalDirPath, penaltyPath, i):
    nanoporeRead = pyfastx.Fasta(nanoporeReadPath)
    commandExecuted = list(map(polishSeq, chunkBarcodeWithReadIter, repeat(nanoporeRead), repeat(tempDirPath), repeat(finalDirPath), repeat(penaltyPath))) 
    commandExecuted = ' ;\\\n'.join(commandExecuted)
    os.system(commandExecuted)
    if i % 100 == 0 :
        logger.info(f'{i*100} reads processed')


@click.command()
@click.option("-i", "MISMATCH_RESULT", help="step07 output")
@click.option("-f", "NANOPORE_READ", help="original nanopore read")
@click.option("-T", "TEMP_DIR", help="temp dir; end with /")
@click.option("-F", "FINAL_DIR", help="polished read stored dir; end with /")
@click.option("-o", "POLISHED_READ", help="polished read")
@click.option("-p", "PENALTY_PATH", help="penalty matrix used by poa")
@click.option("-t", "THREADS", type=int, help="threads")
def main(
    MISMATCH_RESULT,
    NANOPORE_READ,
    TEMP_DIR,
    FINAL_DIR,
    POLISHED_READ,
    THREADS,
    PENALTY_PATH,
):
    os.mkdir(TEMP_DIR)
    os.mkdir(FINAL_DIR)
    logger.info('read mismatch results')
    mismatchResult = pd.read_feather(MISMATCH_RESULT)
    logger.info('prepare for polish')
    mismatchResult["readStrand"] = (
        mismatchResult["readStrand"] ^ mismatchResult["umiStrand"]
    )
    mismatchResult.drop("umiStrand", axis=1, inplace=True)
    mismatchResult["readStrand"] = mismatchResult["readStrand"].astype(str)
    mismatchResult["temp"] = mismatchResult["name"] + "_" + mismatchResult["readStrand"]
    sameUmiReadDt = mismatchResult.groupby("qseqid")["temp"].agg(lambda x: list(x))
    sameUmiReadDc = {i: [[k] for k in j] for i, j in sameUmiReadDt.items()}
    logger.info('start polish')
    umiReadDcIter = chunked(sameUmiReadDc.items(), 100)
    i = 0
    with ProcessPoolExecutor(THREADS) as multiP:
        for umiReadDtChunk in umiReadDcIter:
            i += 1
            multiP.submit(chunkPolishSeq, umiReadDtChunk, NANOPORE_READ, TEMP_DIR, FINAL_DIR, PENALTY_PATH, i)
    logger.info('merge all polished reads')
    time.sleep(10)
    os.system(
        f"""
    cat {FINAL_DIR}/* > {POLISHED_READ} && sleep 15 &&\
    rm -rf {FINAL_DIR} &&\
    rm -rf {TEMP_DIR}
    """
    )


main()
