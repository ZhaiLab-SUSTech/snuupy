import pickle
import glob
import os
import lmdb
import pandas as pd
import numpy as np
from collections import ChainMap
from more_itertools import chunked
from concurrent.futures import ProcessPoolExecutor as multiP
from concurrent.futures import ThreadPoolExecutor as multiT
from jpy_tools.otherTools import Jinterval as jI
import click


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


@click.command()
@click.option("-i", "tsvPath")
@click.option("-o", "lmdbPath")
@click.option("-t", "threads", type=int)
def main(tsvPath, lmdbPath, threads):
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


main()
