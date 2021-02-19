"""
Author: liuzj
Date: 2021-02-19 09:47:43
LastEditors: liuzj
LastEditTime: 2021-02-19 09:48:31
Description: file content
FilePath: /snuupy/snuupy/scripts/parseIllumina.py
"""

import pandas as pd
import numpy as np
import pysam
import h5py


def parseIllumina(
    secondBam, secondIndex, genomeIndex, windowSize, parsedIndex, useColumn
):
    """
    parse Illumina bam file and generate Illumina index.

    secondBam: cellranger result.
    secondIndex: the filtered barcode list which is output by cellranger; format tsv NOT tsv.gz!!!
    genomeIndex: the fai format file of genome
    windowSize: window size; must same with Nanopore
    parsedIndex: output file; end with .index
    useColumn: chromosome counts; except mitochondria and chloroplast
    """
    genomeIndex = pd.read_csv(genomeIndex, sep="\t", header=None).iloc[
        :useColumn
    ]  # remove mito and chlo
    secondIndex = pd.read_csv(secondIndex, header=None)
    secondIndex[0] = secondIndex[0].str.split("-", expand=True)[0]
    secondIndex = set(secondIndex[0])

    # initialize index content
    indexContentDict = {}
    for chrName, chrLength in zip(genomeIndex[0], genomeIndex[1]):
        indexContentDict[chrName] = {}
        i = 0
        while i * windowSize < chrLength:
            indexContentDict[chrName][i] = set()
            i += 1
    # parse cellranger bam
    secondBam = pysam.AlignmentFile(secondBam)
    for read in secondBam:
        if read.reference_name in indexContentDict.keys():
            if read.has_tag("CB") & read.has_tag("UB"):
                readBc = read.get_tag("CB").split("-")[0]
                if readBc in secondIndex:
                    readUmi = read.get_tag("UB")
                    subIndexDict = indexContentDict[read.reference_name]
                    readPos = int(read.reference_start / windowSize)
                    subIndexDict[readPos].add(readBc + readUmi)

    # store parsed result
    indexH5 = h5py.File(parsedIndex, "w")
    for chrNum in indexContentDict.keys():
        for chrNumCount in indexContentDict[chrNum].keys():
            indexH5[f"/{chrNum}/{chrNumCount}"] = [
                x.encode() for x in indexContentDict[chrNum][chrNumCount]
            ]
    indexH5.close()