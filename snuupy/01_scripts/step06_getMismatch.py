'''
@Date: 2020-06-10 12:27:26
@LastEditors: liuzj
@LastEditTime: 2020-06-16 12:49:03
@Description: 利用bio.aligen解析blast结果 获得mismatch
@Author: liuzj
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/00_pipeline/finalVersion/step06_getMismatch.py
'''
import click
import pysam
import pandas as pd
import jpy_tools.ReadProcess as jrp
from Bio import Align
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from more_itertools import chunked

def getAlignScore(line):
    barcodeUmi = line[0]
    unmappedSeq = line[1]
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.gap_score = -2
    aligner.query_gap_score = -1
    aligner.target_end_gap_score  = 0

    barcodeUmi = barcodeUmi.split('_')
    barcode = barcodeUmi[0]
    umi = barcodeUmi[1]
    barcodeUmi = ''.join(barcodeUmi)

    mappingResult = [aligner.align(barcodeUmi, x) for x in unmappedSeq]
    mappingScore = [x.score for x in mappingResult]
    barcodeUmiScore = max(mappingScore)
    bestScoreIndex = mappingScore.index(barcodeUmiScore)
    if bestScoreIndex % 2 == 0:
        mappingStrand = 0
    else:
        mappingStrand = 1

    bestAlign = mappingResult[bestScoreIndex ][0]
    seqAlignedSeq = bestAlign.query[bestAlign.aligned[1][0][0]:bestAlign.aligned[1][-1][-1]]

    umiScore = aligner.align(umi, seqAlignedSeq).score
    barcodeScore = aligner.align(barcode, seqAlignedSeq).score
    return [str(x) for x in [barcodeUmiScore, barcodeScore, umiScore, mappingStrand]]


@click.command()
@click.option('-i', 'MAPPING_RESULT', help = 'mergerd blast result')
@click.option('-b', 'ADD_SEQ_BAM', help = 'bam added unmapped seq tag')
@click.option('-o', 'OUT_FEATHER', help = 'output feather')
@click.option('-t', 'THREADS', type=int, help = 'threads')
def main(MAPPING_RESULT, ADD_SEQ_BAM, OUT_FEATHER, THREADS):
    blastResult = pd.read_csv(MAPPING_RESULT, sep='\t', header=None,\
                         names =  'qseqid sseqid length qstart qend mismatch gaps bitscore score'.split(' '))
    addSeqBam = pysam.AlignmentFile(ADD_SEQ_BAM,'r')

    seqTransformer = jrp.sequence()
    bamDict = defaultdict(lambda:[])

    for x in addSeqBam:
        readESSeq = x.get_tag('ES')
        readFSSeq = x.get_tag('FS')
        if x.is_reverse:
            readStrand = 1
        else:
            readStrand = 0
        bamDict[x.qname].extend([readESSeq,seqTransformer.reverseComplement(readESSeq),\
        readFSSeq, seqTransformer.reverseComplement(readFSSeq),readStrand])

    blastResult.drop_duplicates(['qseqid','sseqid'], inplace = True)
    blastResult['name'] = blastResult['sseqid'].str.split('_',expand=True)[0]
    blastResult['unmappedSeq'] = blastResult['name'].map(bamDict)
    blastResult['unmappedSeq'], blastResult['readStrand'] = blastResult['unmappedSeq'].str[:4], blastResult['unmappedSeq'].str[4]

    iterBlastResult = chunked(zip(blastResult['qseqid'].values,blastResult['unmappedSeq'].values), THREADS*5000*2)
    alignResult = []
    for chunkBlastResult in iterBlastResult:
        with ProcessPoolExecutor(THREADS) as multiP:
            subAlignResult = multiP.map(getAlignScore, chunkBlastResult, chunksize=5000)
        alignResult.extend(list(subAlignResult))

    blastResult['alignResult'] = alignResult
    blastResult['barcodeUmiScore'] = blastResult['alignResult'].str[0].astype(float)
    blastResult['barcodeScore'] = blastResult['alignResult'].str[1].astype(float)
    blastResult['umiScore'] = blastResult['alignResult'].str[2].astype(float)
    blastResult['umiStrand'] = blastResult['alignResult'].str[3].astype(int)

    blastResult['barcodeUmiMismatch'] = (26 - blastResult['barcodeUmiScore'])/2
    blastResult['barcodeMismatch'] = (16 - blastResult['barcodeScore'])/2
    blastResult['umiMismatch'] = (10 - blastResult['umiScore'])/2

    blastResult = blastResult.reindex(['name','qseqid','barcodeUmiMismatch','barcodeMismatch','umiMismatch', 'readStrand', 'umiStrand'], axis=1, copy=False)
    blastResult.reset_index(drop=True, inplace=True)

    blastResult.to_feather(OUT_FEATHER)


main()
