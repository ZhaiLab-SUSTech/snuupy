"""
Description: 
Author: Liuzj
Date: 2020-10-14 15:44:39
LastEditTime: 2020-10-14 19:36:17
LastEditors: Liuzj
"""
import pysam
import sh
import os
import pandas as pd
import pysam
import click
from loguru import logger
from .polyACallerDir.adapterFinder import main as adapterFinder
from .polyACallerDir.PolyACaller import main as polyACaller


def addPolyATag(
    infile,
    genome,
    threads,
    f5dir,
    f5summary,
    bed,
    tempDir,
    fp,
    ep,
    featherPath,
    bamFilePath,
    addPolyAFilePath,
    polyATag,
    minimapPath,
    samtoolsPath
):
    """
    polyACaller and add polyA tag
    """
    try:
        sh.mkdir(tempDir)
    except:
        logger.warning(f"{tempDir} existed!")

    outBam = tempDir + "outputBam.bam"
    outAdapter = tempDir + "adapter.tsv"
    polyACallerResults = tempDir + "polyACaller.tsv"

    if os.path.exists(polyACallerResults):
        pass
    else:
        os.system(f"{minimapPath} -t {threads} -ax splice --secondary=no -G 12000 {genome} {infile} |\
            {samtoolsPath} sort -@ {threads} -o {outBam} - &&\
                {samtoolsPath} index -@ {threads} {outBam}")

        adapterFinder(outBam, infile, outAdapter, threads)
        polyACaller(outAdapter, f5summary, f5dir, polyACallerResults, threads)

    addPolyAResult = (
        pd.read_table(polyACallerResults)
        .reindex(
            [
                "read_core_id",
                "polya_start_base",
                "polya_end_base",
                "polya_length",
                "polya_type",
            ],
            axis=1,
        )
        .rename(
            {
                "read_core_id": "readId",
                "polya_start_base": "tailBaseStart",
                "polya_end_base": "tailBaseEnd",
                "polya_length": "tailLength",
                "polya_type": "readType",
            },
            axis=1,
        ).assign(readId = lambda df:df['readId'].str.split(',').str[0])
        .set_index("readId", drop=False)
        .rename_axis("qname")
    )

    umiReadMapDt = pd.read_feather(featherPath)
    mapNameToId = (
        umiReadMapDt.loc[:, ["qseqid", "name"]].set_index("name").to_dict()["qseqid"]
    )
    umiReadMapDt = (
        umiReadMapDt.groupby("qseqid")["name"].agg(lambda x: list(x)).to_dict()
    )
    umiLabel = set(umiReadMapDt.keys())

    addPolyAResult = addPolyAResult.query('readType not in  ["invalid", "non-polyA/T"]')
    addPolyAResult["umiBarcode"] = addPolyAResult.index.map(mapNameToId)
    addPolyAResult = addPolyAResult.loc[:, ["umiBarcode", "tailLength"]]
    addPolyAResult = addPolyAResult.groupby("umiBarcode")["tailLength"].agg("mean")
    addPolyAResult = addPolyAResult.reindex(umiLabel)
    addPolyAResult.fillna(0, inplace=True)
    addPolyAResult = addPolyAResult.to_dict()

    bamFile = pysam.AlignmentFile(bamFilePath)
    outBamFile = pysam.AlignmentFile(addPolyAFilePath, "wb", template=bamFile)

    for read in bamFile:
        readUmiBarcode = '_'.join(read.qname.split('_')[:2])
        polyALength = addPolyAResult[readUmiBarcode]
        read.set_tag(polyATag, polyALength, "f")
        outBamFile.write(read)
    outBamFile.close()

    pysam.index(addPolyAFilePath)
