"""
Description: 
Author: Liuzj
Date: 2020-10-13 16:01:15
LastEditTime: 2020-10-13 16:06:45
LastEditors: Liuzj
"""
import os
import h5py


def generateOneWindow(chrPath, barcodes):
    os.mkdir(chrPath)
    windowBcCounts = len(barcodes)
    subWindowUpper = 500
    subWindowCounts = windowBcCounts // subWindowUpper + 1
    for subWindowNum in range(subWindowCounts):
        subWindowPath = f"{chrPath}{subWindowNum}.fa"
        subWindowBarcodes = barcodes[
            subWindowNum * subWindowUpper : (subWindowNum + 1) * subWindowUpper
        ]
        if len(subWindowBarcodes) > 0:
            generateOneSubWindow(subWindowPath, subWindowBarcodes)


def generateOneSubWindow(subWindowPath, barcodes):
    with open(subWindowPath, "w") as fh:
        for singleBc in barcodes:
            singleName = f">{singleBc[:16]}_{singleBc[16:]}"
            singleContent = singleName + "\n" + singleBc + "\n"
            fh.write(singleContent)


def generateIlluminaWindow(ILLUMINA_INDEX, OUT_DIR):
    os.system(f"mkdir {OUT_DIR}")
    parsedIndex = h5py.File(ILLUMINA_INDEX, "r")
    for chr_ in parsedIndex.keys():
        currentChr = OUT_DIR + chr_ + "/"
        os.system(f"mkdir {currentChr}")
        for singleWindow in parsedIndex[chr_].keys():
            currentWindowPath = currentChr + singleWindow + "/"
            windowBarcodes = parsedIndex[f"{chr_}/{singleWindow}"][()].astype(str)
            if len(windowBarcodes) > 0:
                generateOneWindow(currentWindowPath, windowBarcodes)