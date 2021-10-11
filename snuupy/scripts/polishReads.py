import pandas as pd
import numpy as np
import os
from loguru import logger
from itertools import repeat
import typing
from more_itertools import chunked
from concurrent.futures import ProcessPoolExecutor
import pyabpoa as poa
import mappy as mm
import time
from .tools import FastaContent, writeFasta


def consensusByPoa(readLs, umi, outputPath, readCounts) -> str:
    seqLs = [x.seq for x in readLs]
    poaAligner = poa.msa_aligner(match=5, extra_b=-1)
    poaRes = poaAligner.msa(seqLs, out_cons=True, out_msa=False)
    poaConsSeq = poaRes.cons_seq[0]
    with open(outputPath, "w") as fh:
        print(f">{umi}_{readCounts}\n{poaConsSeq}", file=fh)
    return poaConsSeq


def mappingByMinimap2(consSeq, readLs, umi, outputPath, readCounts) -> None:
    mmAlign = mm.Aligner(seq=consSeq, preset="map-ont")
    with open(outputPath, "w") as fh:
        for read in readLs:
            for hit in mmAlign.map(read.seq):
                print(
                    f"{read.name}\t{len(read.seq)}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\t{umi}_{readCounts}\t{hit.ctg_len}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.mapq}",
                    file=fh,
                )


def processOneUmi(
    umiWithReadId: typing.Tuple[str, typing.List[str]],
    nanoporeReadContent: FastaContent,
    tempDirPath: str,
    finalDirPath: str,
    raconPath: str,
) -> str:
    umi, readIdWithStrandLs = umiWithReadId
    readCounts = len(readIdWithStrandLs)

    readIdLs = ["_".join(x.split("_")[:-1]) for x in readIdWithStrandLs]
    readStrandLs = [x.split("_")[-1] for x in readIdWithStrandLs]
    readLs = []
    for readId, readStrand in zip(readIdLs, readStrandLs):
        if readStrand == "1":
            readLs.append(nanoporeReadContent[readId].getAnti())
        elif readStrand == "0":
            readLs.append(nanoporeReadContent[readId])
        else:
            logger.error(f"can't identify strand information of {readId}")

    readFinalOutputPath = f"{finalDirPath}{umi}_{readCounts}.final.fa"
    umiAllReadPath = f"{tempDirPath}{umi}_all.fa"

    if readCounts == 1:
        read = readLs[0]
        read.name = f"{umi}_{readCounts}"
        with open(readFinalOutputPath, "w") as fh:
            writeFasta(readLs[0], fh)
        cmdStr = 'true'
    else:
        readLs = readLs[:10]
        with open(umiAllReadPath, "w") as fh:
            for read in readLs:
                writeFasta(read, fh)

        poaConsusReadPath = f"{tempDirPath}{umi}_{readCounts}_poa.fa"
        poaConsSeq = consensusByPoa(readLs, umi, poaConsusReadPath, readCounts)

        minimapPafPath = f"{tempDirPath}{umi}_{readCounts}_minimap2.paf"
        mappingByMinimap2(poaConsSeq, readLs, umi, minimapPafPath, readCounts)
        raconReadPath = f"{tempDirPath}{umi}_{readCounts}_racon.fa"

        cmdStr = f"{raconPath} {umiAllReadPath} {minimapPafPath} {poaConsusReadPath} > {raconReadPath} 2>/dev/null &&\
            mv {raconReadPath} {readFinalOutputPath} &&\
            rm {umiAllReadPath} {minimapPafPath} {poaConsusReadPath}"

    return cmdStr


def processOneChunk(
    umiWithReadIdLs: typing.Sequence[typing.Tuple[str, typing.List[str]]],
    nanoporeFaPath: str,
    tempDirPath: str,
    finalDirPath: str,
    raconPath: str,
    i: int,
    totalUmiCounts:int
) -> None:
    nanoporeReadFastaContent = FastaContent(nanoporeFaPath, True)
    commandExecutedCmdLs = []
    for umiWithReadId in umiWithReadIdLs:
        try:
            commandExecutedCmdLs.append(processOneUmi(umiWithReadId, nanoporeReadFastaContent, tempDirPath, finalDirPath, raconPath))
        except:
            logger.warning(f"detect error when process {umiWithReadId[0]}")
    commandExecutedCmd = " ;\\\n".join(commandExecutedCmdLs)
    os.system(commandExecutedCmd)
    logger.info(f"{i*64} UMIs processed; total {totalUmiCounts}")


def polishReads(
    barcodeAssignedPath,
    raconPath,
    seqkitPath,
    nanoporeFaPath,
    tempResultsDir,
    finalResultsDir,
    threads,
    polishedRead,
):
    if not os.path.exists(tempResultsDir):
        os.mkdir(tempResultsDir)
    if not os.path.exists(finalResultsDir):
        os.mkdir(finalResultsDir)

    logger.info('start build fasta lmda database')
    fastaContent = FastaContent(nanoporeFaPath, True) # ensure lmdb database of fasta is built

    logger.info('start parse barcodeAssigned result')
    barcodeAssignedDf = pd.read_feather(barcodeAssignedPath)

    barcodeAssignedDf["readStrand"] = (
        barcodeAssignedDf["readStrand"] ^ barcodeAssignedDf["umiStrand"]
    ).astype(str)
    barcodeAssignedDf.drop("umiStrand", axis=1, inplace=True)
    barcodeAssignedDf["temp"] = (
        barcodeAssignedDf["name"] + "_" + barcodeAssignedDf["readStrand"]
    )
    sameUmiReadDt = (
        barcodeAssignedDf.groupby("qseqid")["temp"].agg(lambda x: list(x)).to_dict()
    )

    logger.info('start polish')
    totalUmiCounts = len(sameUmiReadDt)
    umiWithReadIdLsIter = chunked(sameUmiReadDt.items(), 64)
    allResults = []
    with ProcessPoolExecutor(threads) as multiP:
        for i, umiWithReadIdLs in enumerate(umiWithReadIdLsIter):
            allResults.append(
                multiP.submit(
                    processOneChunk,
                    umiWithReadIdLs,
                    nanoporeFaPath,
                    tempResultsDir,
                    finalResultsDir,
                    raconPath,
                    i,
                    totalUmiCounts
                )
            )
    [x.result() for x in allResults]

    logger.info("merge all polished reads")
    time.sleep(10)
    
    os.system(
        f"""
    cat {finalResultsDir}* | {seqkitPath} seq -rp > {polishedRead} && sleep 15 &&\
    mkdir /tmp/empty &&\
    rsync --delete-before -av /tmp/empty/ {finalResultsDir}/ &&\
    rsync --delete-before -av /tmp/empty/ {tempResultsDir}/ 
    """
    )
