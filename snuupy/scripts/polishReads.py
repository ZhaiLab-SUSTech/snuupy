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
from .tools import readFasta


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


def getConsensesFasta(readLs, tempPath, penaltyPath, poaPath, finalPath=False):
    getConsensusSeq = GetConsensusSeq()
    
    with open(f'{tempPath}_all.fa', 'w') as fh:
        for read in readLs[:10]:
            fh.write(f'>{read[0]}\n{read[1]}\n')
    
    os.system(f'{poaPath} -read_fasta {tempPath}_all.fa -pir {tempPath}_all.aln {penaltyPath}')
    
    
    consensusSeq = getConsensusSeq(f'{tempPath}_all.aln')
    
    readName = tempPath.split('/')[-1]
    if not finalPath:
        finalPath = f'{tempPath}_round0.fa'

    with open(finalPath, 'w') as fh:
        fh.write(f'>{readName}_{len(readLs)}\n{consensusSeq}\n')
    
#     os.system(f'rm {tempPath}*')

def getPolishRead(tempPath, finalPath, minimapPath, raconPath, times=1):
    def _getPolishCommand(commandExecuted, targetFasta, useFasta, polishedFasta):
        # commandExecuted += f'bwa index -a is {targetFasta} && bwa mem {targetFasta} {useFasta} > {targetFasta}.sam && racon {useFasta} {targetFasta}.sam {targetFasta} > {polishedFasta} &&'
        commandExecuted += f'{minimapPath} --secondary=no -ax map-ont {targetFasta} {useFasta} > {targetFasta}.sam && {raconPath} {useFasta} {targetFasta}.sam {targetFasta} > {polishedFasta} &&'
        return commandExecuted.strip()
    
    
    allPolishedPath = [tempPath + f'_round{x}.fa' for x in range(times + 1)]
    commandExecuted = '('
    for x in range(times):
        commandExecuted = _getPolishCommand(commandExecuted, allPolishedPath[x], f'{tempPath}_all.fa', allPolishedPath[x+1])
        
    commandExecuted += f'mv {allPolishedPath[x+1]} {finalPath} && rm {tempPath}* )'
    commandExecuted = commandExecuted.strip()
    return commandExecuted

def polishSeq(barcodeWithReadLs, nanoporeDict, tempDirPath, finalDirPath, penaltyPath, minimapPath, poaPath, raconPath):
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
        
        getConsensesFasta(readLs, tempReadPath, penaltyPath, poaPath)
        return getPolishRead(tempReadPath, finalPath,  minimapPath, raconPath)

def chunkPolishSeq(chunkBarcodeWithReadIter, nanoporeReadPath, tempDirPath, finalDirPath, penaltyPath, i, minimapPath, poaPath, raconPath):
    nanoporeRead = pyfastx.Fasta(nanoporeReadPath)
    commandExecuted = list(map(polishSeq, chunkBarcodeWithReadIter, repeat(nanoporeRead), repeat(tempDirPath), repeat(finalDirPath), repeat(penaltyPath), repeat(minimapPath), repeat(poaPath), repeat(raconPath)))
    commandExecuted = ' ;\\\n'.join(commandExecuted)
    os.system(commandExecuted)
    if i % 100 == 0 :
        logger.info(f'{i*100} reads processed')



def polishReads(
    MISMATCH_RESULT,
    NANOPORE_READ,
    TEMP_DIR,
    FINAL_DIR,
    POLISHED_READ,
    THREADS,
    PENALTY_PATH,
    minimapPath,
    poaPath,
    raconPath,
    seqkitPath
):  
    if os.path.exists(TEMP_DIR):
        logger.warning(f"{TEMP_DIR} existed!!")
    else:
        os.mkdir(TEMP_DIR)
    if os.path.exists(FINAL_DIR):
        logger.warning(f"{FINAL_DIR} existed!!")
    else:
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
    allResults = []
    with ProcessPoolExecutor(THREADS) as multiP:
        for umiReadDtChunk in umiReadDcIter:
            i += 1
            allResults.append(multiP.submit(chunkPolishSeq, umiReadDtChunk, NANOPORE_READ, TEMP_DIR, FINAL_DIR, PENALTY_PATH, i, minimapPath, poaPath, raconPath))
    [x.result() for x in allResults]
    logger.info('merge all polished reads')
    time.sleep(10)
    os.system(
        f"""
    cat {FINAL_DIR}* | {seqkitPath} seq -rp > {POLISHED_READ} && sleep 15 &&\
    rm -rf {FINAL_DIR} &&\
    rm -rf {TEMP_DIR}
    """
    )