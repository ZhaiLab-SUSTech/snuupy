import pysam
import click
import pyfastx
import itertools
import numpy as np
import pandas as pd
import more_itertools as mlt
from concurrent.futures import ThreadPoolExecutor
from .tools import isOne, getBlock

def isExceedExtend(read, introns):
    """
    判断最后一个intron是否过长
    100 没有intron
    00 intron正常
    10 intron 5'异常
    01 intron 3'异常
    11 均异常
    """
    if len(introns) == 0:
        return 100
    else:
        exons = np.array(getBlock(read, introns))
        introns = np.array(introns)
        exonLength = exons[:, 1] - exons[:, 0]
        intronLength = introns[:, 1] - introns[:, 0]
        result = 0
        if exonLength[-1] / intronLength[-1] <= 0.01:
            result += 1
        if exonLength[0] / intronLength[0] <= 0.01:
            result += 10
        return result


def getClipLength(cigar, exceedExtend, pos):
    """
    用于取出read未比对上基因组的区域并额外取30nt。
    params:
        cigar: pysam.cigar
        exceedExtend: isExceedExtend result
        pos: 0代表3' 1代表5'
    """
    if pos == 0:
        cigar = cigar[::-1]
    if isOne(exceedExtend, pos):
        Length = 0
        for singleCigar in cigar:
            if singleCigar[0] == 3:
                break
            Length += singleCigar[1]
        Length += 30
    else:
        if (cigar[0][0] == 4) | (cigar[0][0] == 5):
            Length = cigar[0][1] + 30
        else:
            Length = 30
    return Length


def getFasta(seq, length):
    return [seq[:length[0]], seq[-length[-1]:]]


def singleReadProcess(read, allFasta):
    """
    对一条read进行处理，并获得加上TAG的read
    params:
        read:pysam.read
        allFasta:diction
    """
    name = read.reference_name
    if (name != 'chrC') | (name != 'chrM'):
        introns = list(bamFile.find_introns([read]))
        exceedExtend = isExceedExtend(read, introns)
        cigar = read.cigar
        fiveLength = getClipLength(cigar, exceedExtend, 1)
        threeLength = getClipLength(cigar, exceedExtend, 0)

        if (fiveLength > 150) or (threeLength > 150):
            return False

        length = [fiveLength, threeLength]
        seq = allFasta[read.qname].antisense if read.is_reverse else allFasta[
            read.qname].seq
        seq = getFasta(seq, length)
        read.set_tag('JI', exceedExtend)
        read.set_tag('FL', fiveLength)
        read.set_tag('EL', threeLength)
        read.set_tag('FS', seq[0])
        read.set_tag('ES', seq[1])
        return read


def singleThread(chunkReads, allFasta):
    '''
    @description: 用于处理二级线程的read
    @param {type} :
        chunkReads: 迭代器 迭代出每一条read
        allFasta: 字典 key为fastaName value为pyfastx对象
    @return: 
        该线程添加tag后的read
    '''
    readProcessList = []
    for read in chunkReads:
        singleReadProcessResult = singleReadProcess(read, allFasta)
        if singleReadProcessResult:
            readProcessList.append(singleReadProcessResult)
    return readProcessList


def bamProcess(readGenerate, fastaDict):
    '''
    @description: 
        用于生成处理bam文件的二级线程。
    @param {type} 
        readGenerate: 迭代器 迭代出每一个二级线程处理的迭代器
        fastaDict: 字典 key为fastaName value为pyfastx对象
    @return: 
        这个chunk的结果
    '''
    chunkProcessList = []
    with ThreadPoolExecutor(max_workers=5) as multiT:
        for _, chunkReads in enumerate(readGenerate):
            chunkProcessList.append(
                multiT.submit(singleThread, chunkReads, fastaDict))
    chunkProcessList = [singleProcessRead for singleChunk in chunkProcessList \
                       for singleProcessRead in singleChunk.result()]
    return chunkProcessList


def outputProcessedRead(bamFileOut, processedReadList):
    '''
    @description: 
        用于输出bam
    @param {type} 
        bamFileOut: 输出文件句柄
        processReadList: 处理后的read列表
    '''
    for read in processedReadList:
        bamFileOut.write(read)



def addUnmappedBaseTag(BAM_PATH, NANOPORE_FASTA, BAM_PATH_OUT):
    global bamFile
    bamFile = pysam.AlignmentFile(BAM_PATH, 'rb')
    bamFileOut = pysam.AlignmentFile(BAM_PATH_OUT, 'wbu', template=bamFile)

    allFasta = pyfastx.Fasta(NANOPORE_FASTA)
    allFastaDict = {}
    for i, x in enumerate(allFasta):
        allFastaDict[x.name] = x

    allReadGenerate = mlt.chunked(bamFile, 40000)
    splicedReadGenerate = mlt.ichunked(allReadGenerate, 8)

    splicedResult = []
    for singleRawChunk in splicedReadGenerate:
        with ThreadPoolExecutor(max_workers=2) as multiT:
            multiT.submit(outputProcessedRead, bamFileOut, splicedResult)
            splicedResult = multiT.submit(bamProcess, singleRawChunk,
                                          allFastaDict).result()
    outputProcessedRead(bamFileOut, splicedResult)
    bamFileOut.close()
