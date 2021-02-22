"""
Description: 
Author: Liuzj
Date: 2020-10-13 16:20:00
LastEditTime: 2020-10-13 18:12:29
LastEditors: Liuzj
"""


import os
import sh
import click
import pandas as pd
from joblib import Parallel, delayed
from .tools import creatUnexistedDir


def scanRefFasta(refPath):
    listPath = []
    for singleChr in os.listdir(refPath):
        currentChrPath = refPath + singleChr + "/"
        for singleWindow in os.listdir(currentChrPath):
            singleWindowFastas = os.listdir(f"{currentChrPath}/{singleWindow}")
            for singleWindowFasta in singleWindowFastas:
                listPath.append([refPath, singleChr, singleWindow, singleWindowFasta])
    return listPath


def blastMapping(line, blastPath, blastMinScore):
    illuminaFilePath = f"{line.illuminaDir}{line.chr}/{line.window}/{line.fasta}"
    nanoporeFilePath = f"{line.nanoporeDir}{line.chr}/{line.window}.fa"
    name = line.fasta.split(".")[0]
    blastWindowPath = f"{line.blastDir}{line.chr}/{line.window}"
    blastFilePath = f"{line.blastDir}{line.chr}/{line.window}/{name}.result"
    if os.path.exists(nanoporeFilePath):
        creatUnexistedDir(blastWindowPath)
        os.system(
            f"""{blastPath}makeblastdb -dbtype nucl -in {illuminaFilePath} -hash_index -parse_seqids &&\
        {blastPath}blastn -query {nanoporeFilePath} -num_threads 2 -word_size 7 -gapopen 0 -gapextend 2 -penalty -1 -reward 1  -db {illuminaFilePath} -outfmt "6 sseqid qseqid length qstart qend mismatch gaps bitscore score" |\
        sort -k 9nr |\
        awk '$9>={blastMinScore} {{print $0}}' > {blastFilePath} &&\
        rm {illuminaFilePath}.*
        """
        )


def windowBlast(ILLUMINA_DIR, NANOPORE_DIR, BLAST_DIR, THREADS, BLAST_PATH, KIT):
    kitMinScore = {"v2": 26 - (2 * 6), "v3": 28 - (2 * 6)}
    blastMinScore = kitMinScore[KIT]

    ILLUMINA_DIR = str(sh.realpath(ILLUMINA_DIR)).strip() + "/"
    NANOPORE_DIR = str(sh.realpath(NANOPORE_DIR)).strip() + "/"
    BLAST_DIR = str(sh.realpath(BLAST_DIR)).strip() + "/"
    BLAST_PATH = str(sh.realpath(BLAST_PATH)).strip() + "/"
    print(ILLUMINA_DIR, NANOPORE_DIR, BLAST_DIR, THREADS, BLAST_PATH)

    blastPath = BLAST_PATH
    creatUnexistedDir(BLAST_DIR)
    illuminaFastaList = pd.DataFrame(scanRefFasta(ILLUMINA_DIR))
    illuminaFastaList[4] = NANOPORE_DIR
    illuminaFastaList[5] = BLAST_DIR
    illuminaFastaList.columns = [
        "illuminaDir",
        "chr",
        "window",
        "fasta",
        "nanoporeDir",
        "blastDir",
    ]
    illuminaFastaList = illuminaFastaList.sample(frac=1)

    for chr_ in illuminaFastaList["chr"].unique():
        creatUnexistedDir(f"{BLAST_DIR}{chr_}")

    Parallel(n_jobs=THREADS, batch_size=100)(
        delayed(blastMapping)(line, blastPath, blastMinScore)
        for line in illuminaFastaList.itertuples()
    )

    os.system(f"cat {BLAST_DIR}*/*/*.result > {BLAST_DIR}allResult.result")