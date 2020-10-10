'''
Date: 2020-08-12 11:28:09
LastEditors: Liuzj
LastEditTime: 2020-10-08 12:50:57
Description: file content
Author: liuzj
FilePath: /liuzj/scripts/pipeline/extractUsefulBaseForCellranger/scripts/step05_extractSeq.py
'''
import os
import sh
import re
import pyfastx
import glob
import lmdb
import click
import numpy as np
from loguru import logger
from concurrent.futures import ProcessPoolExecutor as multiP
from jpy_tools.ReadProcess import writeFastq, readFastq, getSubFastq

def processOneFastq(singleR1Path, singleR2Path, lmdbPath, outDir, cutoff):
    singleR1File, singleR2File = readFastq(singleR1Path), readFastq(singleR2Path)
    singleR1OutFile, singleR2OutFile = outDir + singleR1Path.split('/')[-1], outDir + singleR2Path.split('/')[-1]
    with lmdb.open(lmdbPath, map_size=1099511627776) as mdbDataBase, open(singleR1OutFile, 'w') as fh1, open(singleR2OutFile, 'w') as fh2:
        mdbFile = mdbDataBase.begin()
        for singleRead1, singleRead2 in zip(singleR1File, singleR2File):
            singleUsefulRegion = mdbFile.get(singleRead1.name.encode())
            if singleUsefulRegion:
                singleUsefulRegion = np.frombuffer(singleUsefulRegion, dtype=int).reshape(-1, 2)
                singleRead2Corrected=getSubFastq(singleRead2, singleUsefulRegion)
                if len(singleRead2Corrected.seq) >= cutoff :
                    writeFastq(singleRead1, fh1)
                    writeFastq(singleRead2Corrected, fh2)


@click.command()
@click.option('-i', 'fastqDir')
@click.option('-o', 'outDir')
@click.option('-l', 'lmdbPath', help = 'lmdbPath')
@click.option('-t', 'threads', type=int)
@click.option('-s', 'splitInput', is_flag=True)
@click.option('-c', 'cutoff', type=int, default=75)
def main(fastqDir, outDir, lmdbPath, threads, splitInput, cutoff):
    try:
        os.mkdir(outDir)
    except:
        logger.warning(f'{outDir} existed!!')
    if not splitInput:
        allR1Path = glob.glob(f'{fastqDir}*R1*')
        allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
    else:

        fastqTemp = outDir + 'tempSplited/'
        try:
            sh.mkdir(fastqTemp)
        except:
            logger.warning(f'{fastqTemp} existed!!')

        allR1Path = glob.glob(f'{fastqDir}*_R1*')
        allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
        allSplitedPath = [fastqTemp + re.search(r'(?<=/)\w+?(?=_R1)', x)[0] + '/' for x in allR1Path]

        if allR1Path[0].endswith('.gz'):
            formatGz = True
        else:
            formatGz = False

        splitedNum = threads // len(allSplitedPath)
        
        if splitedNum <= 1 :
            allR1Path = glob.glob(f'{fastqDir}*R1*')
            allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
            if allR1Path[0].endswith('.gz'):
                logger.error('format gz, please uncompress it.')
                1/0
        else:
            mPResults = []
            with multiP(threads//2) as mP:
                for singleR1Path, singleR2Path, singleSplitedPath in zip(allR1Path, allR2Path, allSplitedPath):
                    mPResults.append(mP.submit(sh.seqkit, "split2", "-f", "-1", singleR1Path, "-2", singleR2Path, p=splitedNum, O=singleSplitedPath, j=2))

            tempAllSplitedR1Path = glob.glob(f'{fastqTemp}*/*R1*')
            tempAllSplitedR2Path = [x.replace('R1', 'R2') for x in tempAllSplitedR1Path]
            sampleId = set([re.search(r'(?<=/)\w+?(?=_L)',x)[0] for x in tempAllSplitedR1Path])

            if len(sampleId) != 1:
                raise NameError("MORE THAN ONE INPUT SAMPLES")
            else:
                sampleId = sampleId.pop()

            i = 0
            formatGzUseThreadContents = []
            for tempSingleSplitedR1Path, tempSingleSplitedR2Path in zip(tempAllSplitedR1Path, tempAllSplitedR2Path):
                i += 1
                if formatGz:
                    sh.mv(tempSingleSplitedR1Path, f'{fastqTemp}{sampleId}_L{i:03}_R1_001.fastq.gz')
                    sh.mv(tempSingleSplitedR2Path, f'{fastqTemp}{sampleId}_L{i:03}_R2_001.fastq.gz')
                    formatGzUseThreadContents.append(sh.gzip('-d', f'{fastqTemp}{sampleId}_L{i:03}_R1_001.fastq.gz', _bg=True))
                    formatGzUseThreadContents.append(sh.gzip('-d', f'{fastqTemp}{sampleId}_L{i:03}_R2_001.fastq.gz', _bg=True))
                else:
                    sh.mv(tempSingleSplitedR1Path, f'{fastqTemp}{sampleId}_L{i:03}_R1_001.fastq')
                    sh.mv(tempSingleSplitedR2Path, f'{fastqTemp}{sampleId}_L{i:03}_R2_001.fastq')
            if formatGz:
                [x.wait() for x in formatGzUseThreadContents]

            for singleTempDir in glob.glob(f'{fastqTemp}*/'):
                sh.rmdir(singleTempDir)

            allR1Path = glob.glob(f'{fastqTemp}*R1*')
            allR2Path = [x.replace('R1', 'R2') for x in allR1Path]
        
    
    allSubProcess = []
    with multiP(threads) as mP:
        for singleR1Path, singleR2Path in zip(allR1Path, allR2Path):
            allSubProcess.append(mP.submit(processOneFastq, singleR1Path, singleR2Path, lmdbPath, outDir, cutoff))
    [x.result() for x in allSubProcess]
    
    if not splitInput:
        pass
    else:
        sh.rm('-rf', fastqTemp)

        
if __name__ == '__main__':
    main()
