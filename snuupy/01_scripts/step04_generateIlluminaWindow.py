'''
@Date: 2020-06-08 12:43:38
@LastEditors: liuzj
@LastEditTime: 2020-06-08 14:18:15
@Description: 用于获得每个window里的illumina read
@Author: liuzj
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/00_pipeline/finalVersion/step04_generateIlluminaWindow.py
'''
import os
import h5py
import click

def generateOneWindow(chrPath, barcodes):
    os.mkdir(chrPath)
    windowBcCounts = len(barcodes)
    subWindowUpper = 500
    subWindowCounts = windowBcCounts // subWindowUpper + 1
    for subWindowNum in range(subWindowCounts):
        subWindowPath = f'{chrPath}{subWindowNum}.fa'
        subWindowBarcodes = barcodes[subWindowNum * subWindowUpper: (subWindowNum+1) * subWindowUpper]
        generateOneSubWindow(subWindowPath, subWindowBarcodes)


def generateOneSubWindow(subWindowPath, barcodes):
    with open(subWindowPath, 'w') as fh:
        for singleBc in barcodes:
            singleName = f'>{singleBc[:16]}_{singleBc[16:]}'
            singleContent = singleName + '\n' + singleBc + '\n'
            fh.write(singleContent)


@click.command()
@click.option('-i', 'ILLUMINA_INDEX', help='step01 output')
@click.option('-o', 'OUT_DIR', help='illumina window dir; must end with /')
def main(ILLUMINA_INDEX, OUT_DIR):
    os.system(f'mkdir {OUT_DIR}')
    parsedIndex = h5py.File(ILLUMINA_INDEX,'r')
    for chr_ in parsedIndex.keys():
        currentChr = OUT_DIR + chr_ + '/'
        os.system(f'mkdir {currentChr}')
        for singleWindow in parsedIndex[chr_].keys():
            currentWindowPath = currentChr + singleWindow + '/'
            windowBarcodes = parsedIndex[f'{chr_}/{singleWindow}'].value.astype(str)
            generateOneWindow(currentWindowPath, windowBarcodes)
main()