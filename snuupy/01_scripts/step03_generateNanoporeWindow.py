'''
@Date: 2020-06-08 11:13:53
@LastEditors: liuzj
@LastEditTime: 2020-06-10 18:45:48
@Description: 用于获得每个window里的nanopore read
@Author: liuzj
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/00_pipeline/finalVersion/step03_generateNanoporeWindow.py
'''
import os
import pysam
import click
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


def getGenomeUpper(genomeIndexPath, window, i=2):
    """
    获得染色体的长度
    i 为忽略最后i条
    """
    genomeIndex = pd.read_table(genomeIndexPath, header=None)
    genomeIndex = genomeIndex.iloc[:-i, 0:2]
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

@click.command()
@click.option('-g', 'GENOME_INDEX', help = 'genome index; format fai')
@click.option('-b', 'BAM_PATH', help ='bam added unmapped base tag ; format bam')
@click.option('-o', 'OUT_PATH', help = 'output dir; should end with / ')
@click.option('-w', 'WINDOW', type=int, help = 'window size')
def main(GENOME_INDEX, BAM_PATH, OUT_PATH, WINDOW):
    os.system(f'samtools index {BAM_PATH}')
    genomeUpper = getGenomeUpper(GENOME_INDEX, WINDOW)
    with ThreadPoolExecutor(5) as multiT:
        os.mkdir(OUT_PATH)
        for chr_, upperLimit in genomeUpper.items():
            bamChrFetch = pysam.AlignmentFile(BAM_PATH).fetch(chr_)
            outputPath = f'{OUT_PATH}{chr_}/'
            multiT.submit(parseOneChr, bamChrFetch, WINDOW, upperLimit, outputPath)
main()