'''
Description: 
Author: Liuzj
Date: 2020-10-14 15:44:39
LastEditTime: 2020-10-14 19:36:17
LastEditors: Liuzj
'''
import pysam
import sh
import pandas as pd
import pysam
import click
from loguru import logger
from .polyACallerDir.buildIndex import build_index
from .polyACallerDir.polyACaller import polyACaller


def addPolyATag(infile, genome, threads, f5dir, f5summary, bed, tempDir, fp, ep,
           featherPath, bamFilePath, addPolyAFilePath, polyATag, minimapPath):
    """
    polyACaller and add polyA tag
    """
    try:
        sh.mkdir(tempDir)
    except:
        logger.warning(f'{tempDir} existed!')
        
    outfile = tempDir + 'polyACaller.index'
    build_index(infile, genome, threads, f5dir, f5summary, bed, outfile,
                minimapPath)

    polyACallerResults = tempDir + 'polyACaller.h5'
    polyACaller(outfile, polyACallerResults, threads, fp, ep)

    umiReadMapDt = pd.read_feather(featherPath)
    mapNameToId = umiReadMapDt.loc[:, ['qseqid', 'name']].set_index(
        'name').to_dict()['qseqid']
    umiReadMapDt = umiReadMapDt.groupby('qseqid')['name'].agg(
        lambda x: list(x)).to_dict()
    umiLabel = set(umiReadMapDt.keys())

    addPolyAResult = pd.read_hdf(polyACallerResults)
    addPolyAResult = addPolyAResult.query(
        'readType not in  ["invalid", "non-polyA/T"]')
    addPolyAResult['umiBarcode'] = addPolyAResult.index.map(mapNameToId)
    addPolyAResult = addPolyAResult.loc[:, ['umiBarcode', 'tailLength']]
    addPolyAResult = addPolyAResult.groupby('umiBarcode')['tailLength'].agg(
        'mean')
    addPolyAResult = addPolyAResult.reindex(umiLabel)
    addPolyAResult.fillna(0, inplace=True)
    addPolyAResult = addPolyAResult.to_dict()

    bamFile = pysam.AlignmentFile(bamFilePath)
    outBamFile = pysam.AlignmentFile(addPolyAFilePath, 'wb', template=bamFile)

    for read in bamFile:
        readUmiBarcode = read.qname[:27]
        polyALength = addPolyAResult[readUmiBarcode]
        read.set_tag(polyATag, polyALength, 'f')
        outBamFile.write(read)
    outBamFile.close()

    pysam.index(addPolyAFilePath)

    