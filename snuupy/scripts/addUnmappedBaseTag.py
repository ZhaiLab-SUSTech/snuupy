import pysam
import click
import itertools
import numpy as np
import pandas as pd
import more_itertools as mlt
from concurrent.futures import ThreadPoolExecutor
from .tools import isOne, getBlock, FastaContent


def isExceedExtend(read, introns):
    """
    Determine if the last intron is too long
    100 intron not detected
    00 intron normal
    10 intron 5' abnormal
    01 intron 3' abnormal
    11 both abnormal
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
    obtaining the region of the unmapped reads and take an additional 30nt
    params:
        cigar: pysam.cigar
        exceedExtend: isExceedExtend result
        pos: 0 represents 3' 1 represents 5'
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
    process one read, and add clipped region to the tag part
    params:
        read:pysam.read
        allFasta:diction
    """
    name = read.reference_name
    if (name != "chrC") | (name != "chrM"):
        introns = list(bamFile.find_introns([read]))
        exceedExtend = isExceedExtend(read, introns)
        cigar = read.cigar
        fiveLength = getClipLength(cigar, exceedExtend, 1)
        threeLength = getClipLength(cigar, exceedExtend, 0)

        if (fiveLength > 180) or (threeLength > 180):  # 150 + 30
            return False

        length = [fiveLength, threeLength]
        seq = (
            allFasta[read.qname].getAnti().seq
            if read.is_reverse
            else allFasta[read.qname].seq
        )
        seq = getFasta(seq, length)
        read.set_tag("JI", exceedExtend)
        read.set_tag("FL", fiveLength)
        read.set_tag("EL", threeLength)
        read.set_tag("FS", seq[0])
        read.set_tag("ES", seq[1])
        return read


def singleThread(chunkReads, allFasta):
    """
    @description: process reads
    @param {type} :
        chunkReads: read generator.
        allFasta: dict. key: fastaName value: pyfastx fasta
    @return:
        list, content is read with clipped region tag
    """
    readProcessList = []
    for read in chunkReads:
        singleReadProcessResult = singleReadProcess(read, allFasta)
        if singleReadProcessResult:
            readProcessList.append(singleReadProcessResult)
    return readProcessList


def bamProcess(readGenerate, fastaDict):
    """
    @description:
        process reads
    @return:
        list, content is read with clipped region tag
    """
    chunkProcessList = []
    with ThreadPoolExecutor(max_workers=5) as multiT:
        for _, chunkReads in enumerate(readGenerate):
            chunkProcessList.append(multiT.submit(singleThread, chunkReads, fastaDict))
    chunkProcessList = [
        singleProcessRead
        for singleChunk in chunkProcessList
        for singleProcessRead in singleChunk.result()
    ]
    return chunkProcessList


def outputProcessedRead(bamFileOut, processedReadList):
    """
    @description:
        write processed bam
    @param {type}
        bamFileOut: file handle to output
        processReadList:  list, content is read with clipped region tag
    """
    for read in processedReadList:
        bamFileOut.write(read)


def addUnmappedBaseTag(BAM_PATH, NANOPORE_FASTA, BAM_PATH_OUT):
    global bamFile
    bamFile = pysam.AlignmentFile(BAM_PATH, "rb")
    bamFileOut = pysam.AlignmentFile(BAM_PATH_OUT, "wbu", template=bamFile)

    allFasta = FastaContent(NANOPORE_FASTA)
    allFastaDict = {}
    for i, x in enumerate(allFasta.iter()):
        allFastaDict[x.name] = x

    allReadGenerate = mlt.chunked(bamFile, 40000)
    splicedReadGenerate = mlt.ichunked(allReadGenerate, 8)

    splicedResult = []
    for singleRawChunk in splicedReadGenerate:
        with ThreadPoolExecutor(max_workers=2) as multiT:
            multiT.submit(outputProcessedRead, bamFileOut, splicedResult)
            splicedResult = multiT.submit(
                bamProcess, singleRawChunk, allFastaDict
            ).result()
    outputProcessedRead(bamFileOut, splicedResult)
    bamFileOut.close()
    pysam.index(BAM_PATH_OUT)
