import pandas as pd
import os
import sh
from tqdm import tqdm
import h5py
from .tools import FastaContent
from more_itertools import chunked
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import math


def getMaxSeedLength(x):
    x = math.log(x // 8, 4)
    x = int(x)
    return x


def deleteIndex(ls_nanoporeFa):
    ls_cmd = [
        f"rm {x}.al1 {x}.bck {x}.bwt {x}.des {x}.lcp {x}.llv {x}.ois {x}.prj {x}.sds {x}.skp {x}.ssp {x}.sti1 {x}.suf {x}.tis"
        for x in ls_nanoporeFa
    ]
    cmd = " && ".join(ls_cmd)
    if os.system(cmd) != 0:
        assert False, "triger error during delete index"


def alignByVmatch(
    ls_nanoporeFa,
    dir_illumina,
    dir_tempOutput,
    dir_vmatch,
    lengthBcUmi,
    totalEd,
    seedLength,
):
    ls_cmd = []
    for path_nanoporeFa in ls_nanoporeFa:
        position = "/".join(path_nanoporeFa.split("/")[-2:]).split(".fa")[0]
        dir_illuminaPos = dir_illumina + "/" + position
        if not os.path.exists(dir_illuminaPos):
            continue
        for num_illumina in os.listdir(dir_illuminaPos):
            path_illumina = dir_illuminaPos + "/" + num_illumina
            path_output = (
                dir_tempOutput
                + "/"
                + position.replace("/", "_")
                + "_"
                + num_illumina
                + ".tsv"
            )
            ls_cmd.append(
                f"vmatch -q {path_illumina}  -l {lengthBcUmi} -showdesc 0 -noevalue -noscore -noidentity -e {totalEd} -seedlength {seedLength} {path_nanoporeFa} | sed '1d' > {path_output}"
            )
    cmd = " && ".join(ls_cmd)
    if os.system(cmd) != 0:
        assert False, "triger error during alignment"


def buildIndex(ls_nanoporeFa, seedLength, dir_vmatch):
    ls_cmd = []
    for path_nanoporeFa in ls_nanoporeFa:
        fa_nanopore = FastaContent(path_nanoporeFa)
        fa_counts = sum([len(x.seq) for x in fa_nanopore.iter()])
        _seedLength = int(min(seedLength, getMaxSeedLength(fa_counts)))
        ls_cmd.append(
            f"{dir_vmatch}/mkvtree -db {path_nanoporeFa} -dna -pl {_seedLength} -allout -indexname {path_nanoporeFa}"
        )

    cmd = " && ".join(ls_cmd)
    if os.system(cmd) != 0:
        assert False, "triger error during build index"


def main(
    dir_illumina,
    dir_nanopore,
    dir_output,
    dir_vmatch,
    threads,
    kit,
    bcEd,
    umiEd,
    seedLength,
):
    totalEd = bcEd + umiEd
    lengthBcUmi = {"v2": 26, "v3": 28}[kit]
    path_outFile = f"{dir_output}/all.tsv"

    ls_allNanoporeFa = []
    for chrName in tqdm(os.listdir(dir_nanopore), "get all nanopore fa"):
        ls_chrPath = [dir_nanopore, chrName]
        for name in os.listdir(f"{dir_nanopore}/{chrName}"):
            ls_allNanoporeFa.append("/".join([*ls_chrPath, name]))

    iterLs_NanoporeFa = chunked(ls_allNanoporeFa, 128)
    iterCounts = (len(ls_allNanoporeFa) // 128) + 1

    with ProcessPoolExecutor(min(threads, 4)) as mtp:
        ls_results = list(
            tqdm(
                mtp.map(
                    buildIndex,
                    iterLs_NanoporeFa,
                    repeat(seedLength),
                    repeat(dir_vmatch),
                ),
                desc="build index",
                total=iterCounts,
            )
        )

    dir_tempOutput = f"{dir_output}/temp/"
    sh.mkdir(dir_output, p=True)
    sh.mkdir(dir_tempOutput, p=True)

    iterLs_NanoporeFa = chunked(ls_allNanoporeFa, 128)
    iterCounts = (len(ls_allNanoporeFa) // 128) + 1

    with ProcessPoolExecutor(threads) as mtp:
        ls_results = list(
            tqdm(
                mtp.map(
                    alignByVmatch,
                    iterLs_NanoporeFa,
                    repeat(dir_illumina),
                    repeat(dir_tempOutput),
                    repeat(dir_vmatch),
                    repeat(lengthBcUmi),
                    repeat(totalEd),
                    repeat(seedLength),
                ),
                desc="alignment",
                total=iterCounts,
            )
        )

    iterLs_NanoporeFa = chunked(ls_allNanoporeFa, 64)
    iterCounts = (len(ls_allNanoporeFa) // 64) + 1

    for _ls_nanoporeFa in tqdm(
        iterLs_NanoporeFa, total=iterCounts, desc="delete index"
    ):
        deleteIndex(_ls_nanoporeFa)

    if os.system(f"cat {dir_tempOutput}* > {path_outFile}") != 0:
        assert False, "triger error during cat all results"
