#!/usr/bin/env python
'''
@Author: liuzj
@Date: 2020-05-27 12:48:48
@LastEditors: liuzj
@LastEditTime: 2020-06-09 20:45:32
@Description: 用于处理第三、四步的输出 获得blast结果
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/00_pipeline/finalVersion/step05_windowBlast.py
'''

import os
import click
import pandas as pd
from joblib import Parallel, delayed
from jpy_tools.otherTools import creatUnexistedDir


def scanRefFasta(refPath):
    listPath = []
    for singleChr in os.listdir(refPath):
        currentChrPath = refPath + singleChr + "/"
        for singleWindow in os.listdir(currentChrPath):
            singleWindowFastas = os.listdir(f'{currentChrPath}/{singleWindow}')
            for singleWindowFasta in singleWindowFastas:
                listPath.append([refPath, singleChr ,singleWindow, singleWindowFasta])
    return listPath


def blastMapping(line):
    illuminaFilePath = f'{line.illuminaDir}{line.chr}/{line.window}/{line.fasta}'
    nanoporeFilePath = f'{line.nanoporeDir}{line.chr}/{line.window}.fa'
    name = line.fasta.split('.')[0]
    blastWindowPath = f'{line.blastDir}{line.chr}/{line.window}'
    blastFilePath = f'{line.blastDir}{line.chr}/{line.window}/{name}.result'
    
    creatUnexistedDir(blastWindowPath)

    os.system(f"""makeblastdb -dbtype nucl -in {illuminaFilePath} -hash_index -parse_seqids &&\
    blastn -query {nanoporeFilePath} -num_threads 2 -word_size 7 -gapopen 0 -gapextend 2 -penalty -1 -reward 1  -db {illuminaFilePath} -outfmt "6 sseqid qseqid length qstart qend mismatch gaps bitscore score" |\
    sort -k 9nr |\
    awk '$9>=14 {{print $0}}' > {blastFilePath} &&\
    rm {illuminaFilePath}.*
    """)


@click.command()
@click.option('-i','ILLUMINA_DIR', help = 'illumina dir; should end with "/"')
@click.option('-n','NANOPORE_DIR', help = 'nanopore dir; should end with "/"')
@click.option('-o','BLAST_DIR', help = 'result dir; should end with "/"')
@click.option('-t','THREADS', type = int, help = 'threads')
def main(ILLUMINA_DIR, NANOPORE_DIR, BLAST_DIR, THREADS):
    os.mkdir(BLAST_DIR)
    illuminaFastaList = pd.DataFrame(scanRefFasta(ILLUMINA_DIR))
    illuminaFastaList[4] = NANOPORE_DIR
    illuminaFastaList[5] = BLAST_DIR
    illuminaFastaList.columns = ['illuminaDir', 'chr', 'window', 'fasta', 'nanoporeDir', 'blastDir']
    illuminaFastaList = illuminaFastaList.sample(frac=1)

    for chr_ in illuminaFastaList['chr'].unique():
        os.mkdir(f'{BLAST_DIR}{chr_}')
        
    Parallel(n_jobs=THREADS, batch_size = 10)(\
    delayed(blastMapping)(\
    line) for line in illuminaFastaList.itertuples())


main()
