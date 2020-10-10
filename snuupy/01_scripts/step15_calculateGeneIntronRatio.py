'''
@Date: 2020-07-27 20:54:32
@LastEditors: liuzj
@LastEditTime: 2020-07-28 19:03:09
@Description: file content
@Author: liuzj
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/02_pipeline0710/01_scripts/step15_calculateGeneIntronRatio.py
'''
import pickle
import click
import pandas as pd
import numpy as np
from functools import reduce

def parseOneRead(line):
    lineIntrons = line.ExonOverlapInfo.split(',')
    lineIntrons = np.array(lineIntrons).astype(int)
    
    lineGeneIntronCounts = line.GeneExonCounts - 1
    intronIrDenominator = np.zeros(lineGeneIntronCounts)

    intronIrDenominator[min(lineIntrons):max(lineIntrons)] = 1
    
    intronIrNumerator = np.zeros(lineGeneIntronCounts)
    
    if pd.isna(line.IntronOverlapInfo):
        return np.array([intronIrNumerator, intronIrDenominator])
    else:
        lineIrIntrons = line.IntronOverlapInfo.split(',')
        for singleLineIrIntron in lineIrIntrons:
            singleLineIrIntron = int(singleLineIrIntron)
            intronIrNumerator[singleLineIrIntron] += 1
        intronIrDenominator = (intronIrDenominator + intronIrNumerator).astype(bool).astype(int)
        return np.array([intronIrNumerator, intronIrDenominator])


def parseOneGene(dtframe):
    intronIrFraction = [parseOneRead(x) for x in dtframe.itertuples()]
    intronIrFraction = reduce(lambda a, b: a + b, intronIrFraction)
    intronIrFraction = np.concatenate(
        [intronIrFraction, (intronIrFraction[0] / intronIrFraction[1]).reshape(1, -1)]
    )
    intronIrFraction = "\t".join([",".join(x) for x in intronIrFraction.astype(str)])
    
    geneReadCounts = len(dtframe)
    geneIr = len(dtframe.dropna())/geneReadCounts
    geneName = dtframe.iloc[0,-1]
    geneIntronCounts = dtframe.iloc[0,2] - 1
    return f'{geneName}\t{geneIntronCounts}\t{geneReadCounts}\t{intronIrFraction}\t{geneIr}\n'


def getAllIrResults(dtframe, outPath):
    with open(outPath, "w") as fh:
        header = [
            "Name",
            "intronCounts",
            "readCoverage",
            "intronCoverage",
            "intronRetentionReadCounts",
            "intronRetentionRatio",
            "readIrRatio",
        ]
        header = "\t".join(header) + "\n"
        fh.write(header)
        dtframeGroupby = iter(dtframe.groupby("geneId"))
        for dtframeChunk in dtframeGroupby:
            chunkContent = parseOneGene(dtframeChunk[1])
            fh.write(chunkContent)

@click.command()
@click.option('-i', 'IN_PATH')
@click.option('-o', 'OUT_PATH')
def main(IN_PATH, OUT_PATH):
    intronInfo = pd.read_table(IN_PATH)
    getAllIrResults(intronInfo, OUT_PATH)

main()