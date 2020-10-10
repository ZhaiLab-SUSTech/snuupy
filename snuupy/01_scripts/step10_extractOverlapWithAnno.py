'''
Description: 
Author: Liuzj
Date: 2020-09-25 10:36:52
LastEditTime: 2020-09-25 10:37:53
LastEditors: Liuzj
'''
import numpy as np
import pandas as pd
import portion
import click
import sh
from collections import defaultdict

BED_NAMES=[
    'Chromosome', 'Start', 'End', 'geneName', 'Score', 'Strand', 
    'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'geneSizes', 'geneStarts', 
    ]

@click.command()
@click.option('-b', 'bedFilePath')
@click.option('-i','inBam')
@click.option('-o', 'outBed')
def main(bedFilePath, inBam, outBed):
    bedFile = pd.read_table(bedFilePath, names=BED_NAMES)

    bedFile['geneSizes'] = bedFile['geneSizes'].map(lambda x: np.fromstring(x, sep=',',dtype=int))
    bedFile['geneStarts'] = bedFile['geneStarts'].map(lambda x: np.fromstring(x, sep=',',dtype=int))
    bedFile['geneStarts']  = bedFile['geneStarts'] + bedFile['Start']
    bedFile['geneEnds'] = bedFile['geneStarts'] + bedFile['geneSizes']
    bedFile['isoformName'] = bedFile['geneName'].str.split('\|\|').str[0]
    bedFile['geneName'] = bedFile['geneName'].str.split('.').str[0]
    allLineExonsInterval = defaultdict(lambda :[])
    i = 0
    for singleLine in bedFile.itertuples():
        i += 1
        linePos = portion.closedopen(0,0)
        for singleExonStart, singleExonEnd in zip(singleLine.geneStarts, singleLine.geneEnds):
            linePos = linePos | portion.closedopen(singleExonStart, singleExonEnd)
        allLineExonsInterval[singleLine.geneName].append(linePos)
        if i % 1000 == 0:
            print(i)

    consensusExonsPos = {}
    i = 0
    for geneName, allIsoformExonPos in allLineExonsInterval.items():
        i += 1
        geneConsensusPos = ~portion.closedopen(0,0)
        for singleIsoformExonPos in allIsoformExonPos:
            geneConsensusPos = geneConsensusPos & singleIsoformExonPos
        consensusExonsPos[geneName] = geneConsensusPos
        if i % 1000 == 0:
            print(i)

    specificExonPos = defaultdict(lambda :[])
    i = 0
    for geneName, allIsoformExonPos in allLineExonsInterval.items():
        i += 1
        singleGeneConsensusPos = consensusExonsPos[geneName]
        for singleIsoformExonPos in allIsoformExonPos:
            singleSpecificExonPos = singleIsoformExonPos - singleGeneConsensusPos
            specificExonPos[geneName].append(singleSpecificExonPos)
        if i % 1000 == 0:
            print(i)

    specificExonStart, specificExonEnd = [],[]
    i = 0
    for singleGeneName, singleGeneSpecExonsPos in specificExonPos.items():
        i += 1
        for singleIsoSpecExonsPos in singleGeneSpecExonsPos:
            if singleIsoSpecExonsPos.empty:
                specificExonStart.append([])
                specificExonEnd.append([])
            else:
                specificExonStart.append([x.lower for x in singleIsoSpecExonsPos])
                specificExonEnd.append([x.upper for x in singleIsoSpecExonsPos])
        if i % 1000 == 0:
            print(i)

    bedFileSpec = bedFile.copy(deep=True)
    specificExonStart = [np.array(x) for x in specificExonStart]
    specificExonEnd = [np.array(x) for x in specificExonEnd]
    bedFileSpec['geneStarts'] = specificExonStart
    bedFileSpec['geneEnds'] = specificExonEnd
    bedFileSpec = bedFileSpec.loc[bedFileSpec['geneStarts'].map(lambda x: False if not x.any() else True)]
    bedFileSpec['geneSizes'] = bedFileSpec.pipe(lambda x: x.geneEnds - x.geneStarts)
    bedFileSpec['geneName'] = bedFileSpec.pipe(lambda x: x.isoformName + '_specific')
    bedFileSpec['Start'] = bedFileSpec['geneStarts'].map(lambda x: x[0])
    bedFileSpec['End'] = bedFileSpec['geneEnds'].map(lambda x: x[-1])
    bedFileSpec['geneStarts'] = bedFileSpec['geneStarts'] - bedFileSpec['Start']
    bedFileSpec['BlockCount'] = bedFileSpec['geneStarts'].map(lambda x: len(x))
    bedFileSpec['geneSizes'] = bedFileSpec['geneSizes'].map(lambda x:','.join(list(x.astype(str)))) + ','
    bedFileSpec['geneStarts'] = bedFileSpec['geneStarts'].map(lambda x:','.join(list(x.astype(str))))+ ','
    bedFileSpec.drop(['geneEnds','isoformName'], axis = 1, inplace=True)
    bedFile = pd.read_table(bedFilePath, names=BED_NAMES)
    bedFileSpec = pd.concat([bedFile, bedFileSpec])
    bedFileSpec.sort_values(['Chromosome', 'Start'], inplace=True)
    bedFileSpecStr = []
    for singleLine in bedFileSpec.iterrows():
        bedFileSpecStr.append('\t'.join(list(singleLine[1].astype('str'))))
    bedFileSpecStr = '\n'.join(bedFileSpecStr)

    sh.bedtools.intersect('-abam', inBam, '-wo', '-s', '-split', '-bed' , b='stdin', _in=bedFileSpecStr, _out=outBed)


main()