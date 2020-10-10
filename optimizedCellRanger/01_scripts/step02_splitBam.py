"""
Date: 2020-08-11 14:09:25
LastEditors: liuzj
LastEditTime: 2020-08-11 14:11:06
Description: file content
Author: liuzj
FilePath: /liuzj/scripts/pipeline/extractUsefulBaseForCellranger/scripts/step01_splitBam.py
"""
import os
import click


@click.command()
@click.option("-i", "bamPath")
@click.option("-o", "splitedDir")
@click.option("-t", "splitedCounts")
@click.option(
    "-p", "picardPath", default="~/softwares/picard/picard.jar", show_default=True
)
def main(bamPath, splitedDir, splitedCounts, picardPath):
    """
    分割bam文件
    """
    os.mkdir(f"{splitedDir}")
    os.system(
        f"java -jar {picardPath} SplitSamByNumberOfReads I={bamPath} OUTPUT={splitedDir} N_FILES={splitedCounts}"
    )


main()
