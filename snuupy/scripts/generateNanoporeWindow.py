'''
Description: 
Author: Liuzj
Date: 2020-10-13 11:56:51
LastEditTime: 2021-02-19 10:04:30
LastEditors: liuzj
'''
import os
import pysam
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def getGenomeUpper(genomeIndexPath, window, i):
    """
    get the length up chromosome
    i: only use <i> chromosome(s)
    """
    genomeIndex = pd.read_table(genomeIndexPath, header=None)
    genomeIndex = genomeIndex.iloc[:i, 0:2]
    genomeIndex.set_index([0], inplace=True)
    genomeIndex = {x:y[1]//window  for x,y in genomeIndex.to_dict('index').items()}
    return genomeIndex


def parseOneReadToWindow(read, window, upperLimit, outputPath):
    name = read.qname
    windowStart = read.reference_start//window
    windowEnd = read.reference_end//window
    
    if windowStart != 0:
        windowStart -= 1
    if windowEnd != upperLimit:
        windowEnd != 1

    unmappedBaseE, unmappedBaseLengthE = read.get_tag('ES'), read.get_tag('EL')
    unmappedBaseF, unmappedBaseLengthF = read.get_tag('FS'), read.get_tag('FL')
    nameE = f'{name}_e_{unmappedBaseLengthE}'
    nameF = f'{name}_f_{unmappedBaseLengthF}'
    contentE = f'>{nameE}\n{unmappedBaseE}\n'
    contentF = f'>{nameF}\n{unmappedBaseF}\n'
    content = contentF + contentE
    for window in range(windowStart, windowEnd+1):
        windowPath = f'{outputPath}{window}.fa'
        with open(windowPath,'a') as fh:
            fh.write(content)


def parseOneChr(bamChrFetch, window, upperLimit, outputPath):
    os.mkdir(outputPath)   
    for read in bamChrFetch:
        parseOneReadToWindow(read, window, upperLimit, outputPath)


def generateNanoporeWindow(GENOME_INDEX, BAM_PATH, OUT_PATH, WINDOW, useColumn):
    """
    output nanopore reads based on mapping info

    GENOME_INDEX: the fai format file.
    useColumn: chromosome counts, should same as <parseIllumina>.
    BAM_PATH: bam added unmapped base tag ; format bam.
    OUT_PAT: output dir; end with /
    WINDOW: window size, same as <parseIllumina>
    """
    os.system(f'samtools index {BAM_PATH}')
    genomeUpper = getGenomeUpper(GENOME_INDEX, WINDOW, useColumn)
    with ThreadPoolExecutor(5) as multiT:
        os.mkdir(OUT_PATH)
        for chr_, upperLimit in genomeUpper.items():
            bamChrFetch = pysam.AlignmentFile(BAM_PATH).fetch(chr_)
            outputPath = f'{OUT_PATH}{chr_}/'
            multiT.submit(parseOneChr, bamChrFetch, WINDOW, upperLimit, outputPath)