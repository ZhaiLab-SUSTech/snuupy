
import pandas as pd
import numpy as np
import pysam
import os
import h5py

def parseIllumina(secondBam, secondIndex, genomeIndex, windowSize, parsedIndex, useColumn):
    genomeIndex = pd.read_csv(genomeIndex, sep='\t',
                              header=None).iloc[:useColumn]  #只保留染色体部分
    secondIndex = pd.read_csv(secondIndex, header=None)
    secondIndex[0] = secondIndex[0].str.split('-', expand=True)[0]
    secondIndex = set(secondIndex[0])

    #初始化index容器
    indexContentDict = {}
    for chrName, chrLength in zip(genomeIndex[0], genomeIndex[1]):
        indexContentDict[chrName] = {}
        i = 0
        while i * windowSize < chrLength:
            indexContentDict[chrName][i] = set()
            i += 1
    #解析cellranger结果
    secondBam = pysam.AlignmentFile(secondBam)
    for read in secondBam:
        if read.reference_name in indexContentDict.keys():
            if read.has_tag('CB') & read.has_tag('UB'):
                readBc = read.get_tag('CB').split('-')[0]
                if readBc in secondIndex:
                    readUmi = read.get_tag('UB')
                    subIndexDict = indexContentDict[read.reference_name]
                    readPos = int(read.reference_start / windowSize)
                    subIndexDict[readPos].add(readBc + readUmi)

    #储存结果
    indexH5 = h5py.File(parsedIndex, 'w')
    for chrNum in indexContentDict.keys():
        for chrNumCount in indexContentDict[chrNum].keys():
            indexH5[f'/{chrNum}/{chrNumCount}'] = [
                x.encode() for x in indexContentDict[chrNum][chrNumCount]
            ]

    #输出结果
    indexH5.close()