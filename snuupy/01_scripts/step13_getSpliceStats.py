'''
@Date: 2020-07-27 16:47:47
LastEditors: Liuzj
LastEditTime: 2020-09-12 10:00:21
@Description: file content
@Author: liuzj
FilePath: /liuzj/scripts/pipeline/calcIrRatioNanopore/scripts/step13_getSpliceStats.py
'''
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor as multiP
from loguru import logger
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
    for singleBlockStart, singleBlockSize in zip(geneBlockStarts,
                                                 geneBlockSizes):
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
    for singleBlockStart, singleBlockSize in zip(readBlockStarts,
                                                 readBlockSizes):
        yield jI(singleBlockStart, singleBlockStart + singleBlockSize)


def getGeneIntron(strand):
    forwordExon = False
    backwordExon = False
    intron = jI(0, 0)
    while True:
        backwordExon = forwordExon
        forwordExon = yield intron
        if backwordExon:
            if strand == "+":
                forwordExonLower = forwordExon.lower
                backwordExonUpper = backwordExon.upper
                intron = jI(backwordExonUpper, forwordExonLower)
            elif strand == "-":
                forwordExonUpper = forwordExon.upper
                backwordExonLower = backwordExon.lower
                intron = jI(forwordExonUpper, backwordExonLower)
            else:
                raise ValueError("strand info error")


def getOverlapIntronAndExon(line, needOverlap):
    strand = line.Strand
    exonGenerator = getGeneExon(line)
    blockGenerator = getReadBlock(line)
    overlapExons = []
    overlapIntrons = []
    if needOverlap:
        overlapIntronsInfo = {}

    currentBlock = next(blockGenerator)
    currentExon = next(exonGenerator)
    exonNum = 0

    intronGenerator = getGeneIntron(strand)
    next(intronGenerator)
    currentIntron = intronGenerator.send(currentExon)

    while True:
        if currentExon & currentBlock:
            overlapExons.append(str(exonNum))
            if (needOverlap) & (exonNum >= 1):
                overlapIntronsInfo[str(exonNum - 1)] = str(
                    currentIntron.getOverlapRatio(currentBlock))
        if currentIntron & currentBlock:
            overlapIntrons.append(str(exonNum - 1))
            if needOverlap:
                overlapIntronsInfo[str(exonNum - 1)] = str(
                    currentIntron.getOverlapRatio(currentBlock))

        try:
            if strand == '+':
                if currentExon.lower >= currentBlock.upper:
                    currentBlock = next(blockGenerator)
                else:
                    currentExon = next(exonGenerator)
                    currentIntron = intronGenerator.send(currentExon)
                    exonNum += 1
            elif strand == '-':
                if currentExon.upper > currentBlock.lower:
                    currentExon = next(exonGenerator)
                    currentIntron = intronGenerator.send(currentExon)
                    exonNum += 1
                else:
                    currentBlock = next(blockGenerator)
        except StopIteration:
            if needOverlap:
                if not overlapIntronsInfo:
                    overlapIntronsInfo = ''
                else:
                    overlapIrIntronsInfo = {x:y for x, y in overlapIntronsInfo.items() if float(y) != 0}
                    overlapIntronsInfo = {x:y for x,y in overlapIntronsInfo.items() if x in overlapExons}
                    overlapIntronsInfo.update(overlapIrIntronsInfo)
                    overlapIntronsInfo = [f'{x}:{y}' for x, y in overlapIntronsInfo.items()]
            break
    if needOverlap:
        return overlapExons, overlapIntrons, overlapIntronsInfo
    else:
        return overlapExons, overlapIntrons


@click.command()
@click.option('-i', 'FILE_PATH')
@click.option('-o', 'OUT_PATH')
@click.option('--ratio',
              'NEED_RATIO',
              is_flag=True,
              help='need retention intron overlap ratio or not')
def main(FILE_PATH, OUT_PATH, NEED_RATIO):
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
    logger.info('read bedtools result')
    bedFile = pd.read_table(FILE_PATH,
                            header=None,
                            names=NAMES,
                            usecols=USECOLS)
    logger.info('read bedtools result over; start transform format')
    bedFile['BlockStarts'] = bedFile['BlockStarts'].map(
        lambda x: np.fromstring(x, sep=',', dtype=int))
    bedFile['BlockSizes'] = bedFile['BlockSizes'].map(
        lambda x: np.fromstring(x, sep=',', dtype=int))
    bedFile['geneBlockSizes'] = bedFile['geneBlockSizes'].map(
        lambda x: np.fromstring(x, sep=',', dtype=int))
    bedFile['geneBlockStarts'] = bedFile['geneBlockStarts'].map(
        lambda x: np.fromstring(x, sep=',', dtype=int))
    with open(OUT_PATH, 'w') as fh:
        if NEED_RATIO:
            header = f'Name\tGeneId\tStrand\tGeneExonCounts\tExonOverlapInfo\tIntronOverlapInfo\tintronOverlapRatioInfo\n'
            fh.write(header)
            i = 0

            for line in bedFile.itertuples():
                if i % 100000 == 0:
                    logger.info(f'{i} lines were processed, waiting')
                i += 1

                lineExonInfo, lineIntronInfo, lineIntronOverlapInfo = getOverlapIntronAndExon(
                    line, NEED_RATIO)
                lineExonInfo = ','.join(lineExonInfo)
                lineIntronInfo = ','.join(lineIntronInfo)
                lineIntronOverlapInfo = ','.join(lineIntronOverlapInfo)
                lineStrand = line.Strand
                lineName = line.Name
                lineGene = line.geneName
                lineGeneExonCounts = line.geneBlockCount
                lineContent = f'{lineName}\t{lineGene}\t{lineStrand}\t{lineGeneExonCounts}\t{lineExonInfo}\t{lineIntronInfo}\t{lineIntronOverlapInfo}\n'
                fh.write(lineContent)

        else:
            header = f'Name\tGeneId\tStrand\tGeneExonCounts\tExonOverlapInfo\tIntronOverlapInfo\n'
            fh.write(header)
            i = 0

            for line in bedFile.itertuples():
                if i % 100000 == 0:
                    logger.info(f'{i} lines were processed, waiting')
                i += 1

                lineExonInfo, lineIntronInfo = getOverlapIntronAndExon(
                    line, NEED_RATIO)
                lineExonInfo = ','.join(lineExonInfo)
                lineIntronInfo = ','.join(lineIntronInfo)
                lineStrand = line.Strand
                lineName = line.Name
                lineGene = line.geneName
                lineGeneExonCounts = line.geneBlockCount
                lineContent = f'{lineName}\t{lineGene}\t{lineStrand}\t{lineGeneExonCounts}\t{lineExonInfo}\t{lineIntronInfo}\n'
                fh.write(lineContent)


main()