'''
Description: 
Author: Liuzj
Date: 2020-10-12 10:23:14
LastEditTime: 2020-10-12 22:18:55
LastEditors: Liuzj
'''
import click
from scripts.runCellRanger import runCellRanger
from scripts.extractExonBases import extractExonBases

@click.group()
def main():
    pass

@main.command('runCellRanger')
@click.option('-p', 'parameterPath', help='cellranger count parameters; \nneed two additional parameters --cellRangerPath outDir; \nyaml format')
def _runCellRanger(parameterPath):
    """
    \b
    run cellRanger count
    """
    runCellRanger(parameterPath)

@main.command('extractExonBases')
@click.option('--bam', 'bamPath', help='cellRanger bam')
@click.option('--temp', 'tempDir', help='Dir stored temp files; end with "\\"')
@click.option('-t', 'threads', type=int, help='threads')
@click.option('--picard', 'picardPath', help='picard path')
@click.option('--bed', 'bedAnnoPath', help='gene annotation file; bed12 format')
@click.option('--buffer', 'bufferSize', default='30G', show_default=True, help='split bam buffer size')
@click.option('--fastq', 'fastqDir', help='Dir stored 10x illumina reads; end with "\\"')
@click.option('--out', 'outDir', help='output dir; end with "\\"')
@click.option('--cutoff', 'cutoff', type=int, default=75, help='reads chopped shorter than <cutoff> will be discarded')
def _extractExonBases(bamPath, tempDir, threads, picardPath, bedAnnoPath, bufferSize, fastqDir, outDir, cutoff):
    """
    extract reads exon mapping region bases
    """
    extractExonBases(bamPath, tempDir, threads, picardPath, bedAnnoPath, bufferSize, fastqDir, outDir, cutoff)

main()