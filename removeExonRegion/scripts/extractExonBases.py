import os
import sh
import re
from loguru import logger
import pickle
import glob
import lmdb
import pyfastx
import pandas as pd
import numpy as np
from collections import ChainMap, namedtuple
from more_itertools import chunked
from concurrent.futures import ProcessPoolExecutor as multiP
from concurrent.futures import ThreadPoolExecutor as multiT
from .tools import Jinterval as jI
from .tools import writeFastq, readFastq, getSubFastq


def splitBam(bamPath, splitedDir, splitedCounts, picardPath):
    sh.mkdir(f"{splitedDir}")
    os.system(
        f"java -jar {picardPath} SplitSamByNumberOfReads I={bamPath} OUTPUT={splitedDir} N_FILES={splitedCounts}"
    )

def getOverlapTsv(bam, bed, tsv):
    os.system(f"bedtools intersect -abam {bam} -b {bed} -split -bed -wo -s > {tsv}")

def getOverlapInfo(splitedDir, tsvPath, threads, bedAnnoPath, bufferSize):
    tempDir = "/".join(splitedDir.split("/")[:-2]) + "/getOverlapTemp/"
    os.mkdir(tempDir)

    allChunkedBam = [splitedDir + x for x in os.listdir(splitedDir)]
    allChunkedOverlapInfo = [
        ".".join(x.split(".")[:-1]) + "_overlap.tsv" for x in allChunkedBam
    ]

    with multiP(threads) as mP:
        for singleBam, singleTsv in zip(allChunkedBam, allChunkedOverlapInfo):
            mP.submit(getOverlapTsv, singleBam, bedAnnoPath, singleTsv)

    os.system(
        f"""
    LC_ALL=C cat {' '.join(allChunkedOverlapInfo)} |\
        LC_ALL=C sort -T {tempDir} -S {bufferSize} --parallel {threads} -k 4,4 -k 25nr,25 |\
            LC_ALL=C sort -T {tempDir} -S {bufferSize} --parallel {threads} -k 4,4 -u > {tsvPath} &&\
                rm -rf {tempDir}
    """
    )

def getGeneExon(line):
    geneBlockStarts = line.geneBlockStarts
    geneBlockSizes = line.geneBlockSizes
    if line.Strand == "+":
        geneBlockStarts = geneBlockStarts + line.geneStart
    #        geneBlockSizes = geneBlockSizes
    else:
        geneBlockStarts = geneBlockStarts[::-1] + line.geneStart
        geneBlockSizes = geneBlockSizes[::-1]
    for singleBlockStart, singleBlockSize in zip(geneBlockStarts, geneBlockSizes):
        yield jI(singleBlockStart, singleBlockStart + singleBlockSize, 0)

def getReadBlock(line):
    readBlockStarts = line.BlockStarts
    readBlockSizes = line.BlockSizes
    if line.Strand == "+":
        readBlockStarts = readBlockStarts + line.Start
    #         readBlockSizes = readBlockSizes
    else:
        readBlockStarts = readBlockStarts[::-1] + line.Start
        readBlockSizes = readBlockSizes[::-1]
    for singleBlockStart, singleBlockSize in zip(readBlockStarts, readBlockSizes):
        yield jI(singleBlockStart, singleBlockStart + singleBlockSize), singleBlockSize


def getGeneExonOverlap(line):

    strand = line.Strand
    overlapInfo = []
    exonGenerator = getGeneExon(line)
    blockGenerator = getReadBlock(line)

    currentBlockPos = 0
    currentBlock, currentBlockSize = next(blockGenerator)
    currentExon = next(exonGenerator)

    while True:
        currentOverlapInfo = currentExon & currentBlock
        if currentOverlapInfo:
            if strand == "+":
                overlapInfo.append(
                    [
                        x - currentBlock.lower + currentBlockPos
                        for x in currentOverlapInfo
                    ]
                )
            else:
                overlapInfo.append(
                    [
                        currentBlock.upper - x + currentBlockPos
                        for x in currentOverlapInfo[::-1]
                    ]
                )
        try:
            if strand == "+":
                if currentExon.lower >= currentBlock.upper:
                    currentBlockPos += currentBlockSize
                    currentBlock, currentBlockSize = next(blockGenerator)
                else:
                    currentExon = next(exonGenerator)
            elif strand == "-":
                if currentExon.upper > currentBlock.lower:
                    currentExon = next(exonGenerator)
                else:
                    currentBlockPos += currentBlockSize
                    currentBlock, currentBlockSize = next(blockGenerator)
        except StopIteration:
            break

    overlapInfo = np.array(overlapInfo).tobytes()

    return overlapInfo


def processOneChunk(bedFile):
    bedFile["BlockStarts"] = bedFile["BlockStarts"].map(
        lambda x: np.fromstring(x, sep=",", dtype=int)
    )
    bedFile["BlockSizes"] = bedFile["BlockSizes"].map(
        lambda x: np.fromstring(x, sep=",", dtype=int)
    )
    bedFile["geneBlockSizes"] = bedFile["geneBlockSizes"].map(
        lambda x: np.fromstring(x, sep=",", dtype=int)
    )
    bedFile["geneBlockStarts"] = bedFile["geneBlockStarts"].map(
        lambda x: np.fromstring(x, sep=",", dtype=int)
    )

    readUsefulRegion = {}
    for line in bedFile.itertuples():
        readUsefulRegion[line.Name] = getGeneExonOverlap(line)

    return readUsefulRegion


def processSeveralChunks(chunkGenerator, threads):
    readUsefulRegionSeveral = []
    severalChunk = next(chunkGenerator)
    with multiP(threads) as mP:
        for bedFile in severalChunk:
            readUsefulRegionSeveral.append(mP.submit(processOneChunk, bedFile))

    readUsefulRegionSeveral = [x.result() for x in readUsefulRegionSeveral]
    readUsefulRegionSeveral = ChainMap(*readUsefulRegionSeveral)
    return readUsefulRegionSeveral


def writeResult(mdbFile, severalChunkResult):
    if severalChunkResult:
        for key, value in severalChunkResult.items():
            mdbFile.put(key=key.encode(), value=value)
    else:
        pass

def getUsefulRegion(tsvPath, lmdbPath, threads):
    NAMES = [
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Score",
        "Strand",
        "ThickStart",
        "ThickEnd",
        "ItemRGB",
        "BlockCount",
        "BlockSizes",
        "BlockStarts",
        "geneChromosome",
        "geneStart",
        "geneEnd",
        "geneName",
        "geneScore",
        "geneStrand",
        "geneThickStart",
        "geneThickEnd",
        "geneItemRGB",
        "geneBlockCount",
        "geneBlockSizes",
        "geneBlockStarts",
        "cov",
    ]
    USECOLS = [
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Strand",
        "BlockSizes",
        "BlockStarts",
        "geneStart",
        "geneEnd",
        "geneName",
        "geneBlockCount",
        "geneBlockSizes",
        "geneBlockStarts",
    ]
    bedFileChunked = chunked(
        pd.read_table(tsvPath, names=NAMES, usecols=USECOLS, chunksize=48*1024), threads
    )

    severalChunkResult = []
    with lmdb.open(lmdbPath, map_size=1099511627776) as mdbDataBase:
        mdbFile = mdbDataBase.begin(write=True)
        while True:
            try:
                with multiT(2) as mT:
                    writeThread = mT.submit(writeResult, mdbFile, severalChunkResult)
                    severalChunkResult = mT.submit(
                        processSeveralChunks, bedFileChunked, threads
                    ).result()
                    writeThread.result()
            except StopIteration:
                break
        mdbFile.commit()

def processOneFastq(singleR1Path, singleR2Path, lmdbPath, outDir, cutoff):
    singleR1File, singleR2File = readFastq(singleR1Path), readFastq(singleR2Path)
    singleR1OutFile, singleR2OutFile = outDir + singleR1Path.split('/')[-1], outDir + singleR2Path.split('/')[-1]
    with lmdb.open(lmdbPath, map_size=1099511627776) as mdbDataBase, open(singleR1OutFile, 'w') as fh1, open(singleR2OutFile, 'w') as fh2:
        mdbFile = mdbDataBase.begin()
        for singleRead1, singleRead2 in zip(singleR1File, singleR2File):
            singleUsefulRegion = mdbFile.get(singleRead1.name.encode())
            if singleUsefulRegion:
                singleUsefulRegion = np.frombuffer(singleUsefulRegion, dtype=int).reshape(-1, 2)
                singleRead2Corrected=getSubFastq(singleRead2, singleUsefulRegion)
                if len(singleRead2Corrected.seq) >= cutoff :
                    writeFastq(singleRead1, fh1)
                    writeFastq(singleRead2Corrected, fh2)

def extractSeq(fastqDir, outDir, lmdbPath, threads, splitInput, cutoff):
    try:
        os.mkdir(outDir)
    except:
        logger.warning(f'{outDir} existed!!')
    if not splitInput:
        allR1Path = glob.glob(f'{fastqDir}*R1*')
        allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
    else:

        fastqTemp = outDir + 'tempSplited/'
        try:
            sh.mkdir(fastqTemp)
        except:
            logger.warning(f'{fastqTemp} existed!!')

        allR1Path = glob.glob(f'{fastqDir}*_R1*')
        allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
        allSplitedPath = [fastqTemp + re.search(r'(?<=/)[\w\W]+?(?=_R1)', x)[0] + '/' for x in allR1Path]

        if allR1Path[0].endswith('.gz'):
            formatGz = True
        else:
            formatGz = False

        splitedNum = threads // len(allSplitedPath)
        
        if splitedNum <= 1 :
            allR1Path = glob.glob(f'{fastqDir}*R1*')
            allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
            if allR1Path[0].endswith('.gz'):
                logger.error('format gz, please uncompress it.')
                1/0
        else:
            mPResults = []
            with multiP(threads//2) as mP:
                for singleR1Path, singleR2Path, singleSplitedPath in zip(allR1Path, allR2Path, allSplitedPath):
                    mPResults.append(mP.submit(sh.seqkit, "split2", "-f", "-1", singleR1Path, "-2", singleR2Path, p=splitedNum, O=singleSplitedPath, j=2))

            tempAllSplitedR1Path = glob.glob(f'{fastqTemp}*/*R1*')
            tempAllSplitedR2Path = [x.replace('R1', 'R2') for x in tempAllSplitedR1Path]
            sampleId = set([re.search(r'(?<=/)[\w\W]+?(?=_L)',x)[0] for x in tempAllSplitedR1Path])

            if len(sampleId) != 1:
                raise NameError("MORE THAN ONE INPUT SAMPLES")
            else:
                sampleId = sampleId.pop()

            i = 0
            formatGzUseThreadContents = []
            for tempSingleSplitedR1Path, tempSingleSplitedR2Path in zip(tempAllSplitedR1Path, tempAllSplitedR2Path):
                i += 1
                if formatGz:
                    sh.mv(tempSingleSplitedR1Path, f'{fastqTemp}{sampleId}_L{i:03}_R1_001.fastq.gz')
                    sh.mv(tempSingleSplitedR2Path, f'{fastqTemp}{sampleId}_L{i:03}_R2_001.fastq.gz')
                    formatGzUseThreadContents.append(sh.gzip('-d', f'{fastqTemp}{sampleId}_L{i:03}_R1_001.fastq.gz', _bg=True))
                    formatGzUseThreadContents.append(sh.gzip('-d', f'{fastqTemp}{sampleId}_L{i:03}_R2_001.fastq.gz', _bg=True))
                else:
                    sh.mv(tempSingleSplitedR1Path, f'{fastqTemp}{sampleId}_L{i:03}_R1_001.fastq')
                    sh.mv(tempSingleSplitedR2Path, f'{fastqTemp}{sampleId}_L{i:03}_R2_001.fastq')
            if formatGz:
                [x.wait() for x in formatGzUseThreadContents]

            for singleTempDir in glob.glob(f'{fastqTemp}*/'):
                sh.rmdir(singleTempDir)

            allR1Path = glob.glob(f'{fastqTemp}*R1*')
            allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
        
    
    allSubProcess = []
    with multiP(threads) as mP:
        for singleR1Path, singleR2Path in zip(allR1Path, allR2Path):
            allSubProcess.append(mP.submit(processOneFastq, singleR1Path, singleR2Path, lmdbPath, outDir, cutoff))
    [x.result() for x in allSubProcess]
    
    if not splitInput:
        pass
    else:
        sh.rm('-rf', fastqTemp)

def extractExonBases(bamPath, tempDir, threads, picardPath, bedAnnoPath, bufferSize, fastqDir, outDir, cutoff, splitInput=True):
    try:
        sh.mkdir(tempDir)
    except:
        logger.warning('{tempDir} existed')

    ### split bam
    logger.info('split bam start')
    splitedDir = f'{tempDir}splitedDir/'
    splitedCounts = threads
    splitBam(bamPath, splitedDir, splitedCounts, picardPath)
    logger.info('split bam end')

    ### get overlap region
    logger.info('get overlap region and generate lmdb database start')
    tsvPath = f'{tempDir}overlapRegion.tsv'
    getOverlapInfo(splitedDir, tsvPath, threads, bedAnnoPath, bufferSize)
    logger.info('get overlap region and generate lmdb database end')

    ### parse overlap region
    logger.info('parse temp file start')
    lmdbPath = f'{tempDir}overlapRegionLmdb/'
    getUsefulRegion(tsvPath, lmdbPath, threads)
    logger.info('parse temp file end')

    ### extract seq
    logger.info('extract exon base start')
    extractSeq(fastqDir, outDir, lmdbPath, threads, splitInput, cutoff)
    logger.info('extract exon base end')

