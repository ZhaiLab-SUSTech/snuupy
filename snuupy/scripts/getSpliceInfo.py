import pandas as pd
import numpy as np
import sh
from concurrent.futures import ProcessPoolExecutor as multiP
from loguru import logger
from .tools import Jinterval as jI
from .tools import bedtoolsGetIntersect
from io import StringIO
import click
import pickle


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
                    overlapIrIntronsInfo = {
                        x: y
                        for x, y in overlapIntronsInfo.items() if float(y) != 0
                    }
                    overlapIntronsInfo = {
                        x: y
                        for x, y in overlapIntronsInfo.items()
                        if x in overlapExons
                    }
                    overlapIntronsInfo.update(overlapIrIntronsInfo)
                    overlapIntronsInfo = [
                        f'{x}:{y}' for x, y in overlapIntronsInfo.items()
                    ]
            break
    if needOverlap:
        return overlapExons, overlapIntrons, overlapIntronsInfo
    else:
        return overlapExons, overlapIntrons


def filterResultsBasedOnGeneName(INTRON_INFO, GENE_NAME_INFO, OUT_PATH):
    with open(GENE_NAME_INFO, 'rb') as fh:
        geneNameInfo = pickle.load(fh)
    intronInfo = pd.read_table(INTRON_INFO)
    intronInfo["geneIdNonTrans"] = intronInfo["GeneId"].str.split(".").str[0]
    intronInfo["baselineGeneId"] = intronInfo["Name"].map(
        lambda x: geneNameInfo.get(x, {"gene_id": 0})["gene_id"])
    intronInfo.query("geneIdNonTrans == baselineGeneId", inplace=True)
    intronInfo.drop(["baselineGeneId", "GeneId"], axis=1, inplace=True)
    intronInfo.rename({"geneIdNonTrans": "geneId"}, axis=1, inplace=True)
    intronInfo.to_csv(OUT_PATH, sep="\t", index=False)



def getSpliceInfo(INBAM_PATH, BED_REPRE_ANNO, GENE_NAME_INFO, OUT_PATH,
                  NEED_RATIO, bedtoolsPath):
    intersectResults = bedtoolsGetIntersect(INBAM_PATH, BED_REPRE_ANNO,
                                            bedtoolsPath)

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
    bedFile = pd.read_table(intersectResults,
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
    fh = StringIO()
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
    logger.info(f'{i} lines were processed, waiting')
    
    fh.seek(0)
    filterResultsBasedOnGeneName(fh, GENE_NAME_INFO, OUT_PATH)
