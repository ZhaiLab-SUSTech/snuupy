"""
Date: 2020-08-11 14:18:58
LastEditors: liuzj
LastEditTime: 2020-08-11 14:27:47
Description: file content
Author: liuzj
FilePath: /liuzj/scripts/pipeline/extractUsefulBaseForCellranger/scripts/step02_getOverlapInfo.py
"""
import os
import click
from concurrent.futures import ProcessPoolExecutor as multiP


def getOverlapTsv(bam, bed, tsv):
    os.system(f"bedtools intersect -abam {bam} -b {bed} -split -bed -wo -s > {tsv}")


@click.command()
@click.option("-i", "splitedDir")
@click.option("-o", "tsvPath")
@click.option("-t", "threads", type=int)
@click.option(
    "--bed", "bedAnnoPath", default="~/data/Araport11/gene.bed", show_default=True
)
@click.option("--bS", "bufferSize", default="30G", show_default=True)
def main(splitedDir, tsvPath, threads, bedAnnoPath, bufferSize):
    tempDir = "/".join(splitedDir.split("/")[:-2]) + "/getOverlapTemp/"
    os.mkdir(tempDir)

    allChunkedBam = [splitedDir + x for x in os.listdir(splitedDir)]
    allChunkedOverlapInfo = [
        ".".join(x.split(".")[:-1]) + "_overlap.tsv" for x in allChunkedBam
    ]

    with multiP(threads) as mP:
        for singleBam, singleTsv in zip(allChunkedBam, allChunkedOverlapInfo):
            mP.submit(getOverlapTsv, singleBam, bedAnnoPath, singleTsv)

    os.system(
        f"""
    LC_ALL=C cat {' '.join(allChunkedOverlapInfo)} |\
        LC_ALL=C sort -T {tempDir} -S {bufferSize} --parallel {threads} -k 4,4 -k 25nr,25 |\
            LC_ALL=C sort -T {tempDir} -S {bufferSize} --parallel {threads} -k 4,4 -u > {tsvPath} &&\
                rm -rf {tempDir}
    """
    )


main()
