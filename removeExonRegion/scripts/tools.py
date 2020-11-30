'''
Description: 
Author: Liuzj
Date: 2020-10-12 21:43:11
LastEditTime: 2020-10-12 22:30:37
LastEditors: Liuzj
'''
import os
import sh
import glob
import pysam
import h5py
import numpy as np
import pandas as pd
from collections import ChainMap, namedtuple
from loguru import logger

class Jinterval:
    def __init__(self, lower, upper, overlapLimit=0.5):
        self.lower, self.upper = lower, upper
        self.interval = [lower, upper]
        self.overlapLimit = overlapLimit

    def __repr__(self):
        return f"Jinterval{self.interval}"

    def __str__(self):
        return f"Jinterval{self.interval}"

    def __and__(self, otherInterval):
        minn = max(self.lower, otherInterval.lower)
        maxn = min(self.upper, otherInterval.upper)
        if (maxn - minn) / (self.upper - self.lower) > self.overlapLimit:
            return [minn, maxn]
        else:
            return False

    def getOverlapRatio(self, otherInterval):
        minn = max(self.lower, otherInterval.lower)
        maxn = min(self.upper, otherInterval.upper)
        return max((maxn - minn) / (self.upper - self.lower), 0)


def writeFastq(read, fh):
    '''
    @description: 用于将pyfastx的read输出为fastq
    @param:
        read: pyfastx fastq
        fh: file fh mode w
    @return: None
    '''
    readContent = f'@{read.name}\n{read.seq}\n{read.desc}\n{read.qual}\n'
    fh.write(readContent)


def readFastq(path, length=False):
    '''
    @description: 读fastq
    @param {type} fastq路径, 读取长度从3'算
    @return: 一个迭代器
    '''
    FastqRead = namedtuple('FastqRead', ['name', 'seq', 'desc', 'qual'])

    def _readFastq(path):
        with open(path, 'r') as fh:
            i = 0
            readContent = []
            while True:
                lineContent = fh.readline()
                if lineContent == '':
                    break
                i += 1
                readContent.append(lineContent.strip())
                if i % 4 == 0:
                    if not length:
                        read = FastqRead(name=readContent[0][1:].split(' ')[0],
                                         seq=readContent[1],
                                         desc=readContent[2],
                                         qual=readContent[3])
                    else:
                        read = FastqRead(name=readContent[0][1:].split(' ')[0],
                                         seq=readContent[1][:length],
                                         desc=readContent[2],
                                         qual=readContent[3][:length])
                    yield read
                    readContent = []

    return _readFastq(path)


def getSubFastq(fastqRead, subRegion):
    """
    subRegion: np.array shape n*2
    """
    FastqRead = namedtuple("FastqRead", ["name", "seq", "desc", "qual"])
    name = fastqRead.name
    desc = fastqRead.desc
    seq = ""
    qual = ""
    for singleSubRegion in subRegion:
        seq += fastqRead.seq[singleSubRegion[0]:singleSubRegion[1]]
        qual += fastqRead.qual[singleSubRegion[0]:singleSubRegion[1]]
    read = FastqRead(name=name, desc=desc, seq=seq, qual=qual)
    return read

