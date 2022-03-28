import sh
import pysam
import tqdm
import os
from loguru import logger
from joblib import Parallel, delayed
import scanpy as sc
import muon as mu
from functools import reduce
from .polyAClusterDetected import polyAClusterDetected
from .generateMtx import generateMtx


def _getAllDir(sql_dir_war):
    dbtype_list = os.listdir(sql_dir_war)
    for dbtype in dbtype_list:
        if os.path.isfile(os.path.join(sql_dir_war, dbtype)):
            dbtype_list.remove(dbtype)
    return dbtype_list


def _concatMd(dtMd, label, indexUnique, mod="inner"):
    dt_sampleMod = {x: set(list(y.mod.keys())) for x, y in dtMd.items()}
    if mod == "inner":
        ls_useMod = list(
            reduce(lambda x, y: list(set(x) & set(y)), list(dt_sampleMod.values()))
        )
    elif mod == "outer":
        ls_useMod = list(
            reduce(lambda x, y: list(set(x) | set(y)), list(dt_sampleMod.values()))
        )
    else:
        assert False, "Unknown mod parameter"
    dtAd_mod = {}
    for mod in ls_useMod:
        dtAd_singleMod = {x: y[mod] for x, y in dtMd.items() if mod in y.mod}
        ad_mod = sc.concat(
            dtAd_singleMod, join="outer", label=label, index_unique=indexUnique
        )
        dtAd_mod[mod] = ad_mod
    md = mu.MuData(dtAd_mod)
    return md


def main(
    dir_allResults,
    dir_output,
    path_genome,
    path_geneBed,
    threads=16,
    ls_sample=None,
    isBed12=False,
    geneTag="gi",
    irMode=True,
    path_usedIntron=None,
    onlyFullLength=True,
    path_samtools="samtools",
    path_bedtools="bedtools",
):

    BAM_RELATIVE_PATH = "./results/step14_addPolyATag/addGNPABam.bam"
    IR_INFO_RELATIVE_PATH = "./results/step13_getSpliceInfo/splicingInfo.tsv"
    ILLUMINA_H5_RELATIVE_PATH = (
        "./results/step1_runCellRanger/test/outs/filtered_feature_bc_matrix.h5"
    )

    sh.mkdir(dir_output, p=True)
    if ls_sample is None:
        ls_sample = [
            x
            for x in _getAllDir(dir_allResults)
            if os.path.exists(f"{dir_allResults}/{x}/results/")
        ]
        logger.info(f"detected sample: {ls_sample}")

    ## merge tagged bam files
    path_mergedBam = dir_output + "/" + "merged.bam"
    path_mergedSortedBam = path_mergedBam + ".sorted.bam"
    lsPath_allTagedBam = [
        f"{dir_allResults}/{x}/{BAM_RELATIVE_PATH}" for x in ls_sample
    ]
    bam_merged = pysam.AlignmentFile(
        path_mergedBam, "w", template=pysam.AlignmentFile(lsPath_allTagedBam[0], "r")
    )
    for sample, path_bam in zip(ls_sample, lsPath_allTagedBam):
        bam_sample = pysam.AlignmentFile(path_bam, "r")
        bamCounts = bam_sample.count()
        bam_sample.reset()
        for read in tqdm.tqdm(bam_sample, sample, total=bamCounts):
            read.qname = sample + "-" + read.qname
            bam_merged.write(read)
        bam_sample.close()
    bam_merged.close()
    sh.Command(path_samtools).sort(
        path_mergedBam, O="BAM", o=path_mergedSortedBam, **{"@": threads}
    )
    sh.Command(path_samtools).index(path_mergedSortedBam, **{"@": threads})
    sh.rm(path_mergedBam)

    ## get pac
    dir_pacResults = dir_output + "/pac/"
    polyAClusterDetected(
        path_genome,
        path_mergedSortedBam,
        path_geneBed,
        dir_pacResults,
        threads,
        isBed12,
        path_bedtools,
    )

    ## generate matrix
    path_pacBed = dir_pacResults + "/polya_cluster.filtered.bed"
    Parallel(n_jobs=4)(
        delayed(generateMtx)(
            path_pacBed,
            f"{dir_allResults}/{sample}/{BAM_RELATIVE_PATH}",
            geneTag,
            f"{dir_allResults}/{sample}/{IR_INFO_RELATIVE_PATH}",
            True,
            None,
            f"{dir_allResults}/{sample}/{ILLUMINA_H5_RELATIVE_PATH}",
            False,
            dir_output + f"/{sample}.h5mu",
        )
        for sample in ls_sample
    )

    ## merge matrix
    dtMd = {
        sample: mu.read_h5mu(dir_output + f"/{sample}.h5mu") for sample in ls_sample
    }
    md = _concatMd(dtMd, label="sample", indexUnique="-")
    md.write(dir_output + "/merged.h5mu")