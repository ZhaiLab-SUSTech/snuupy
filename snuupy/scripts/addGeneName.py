from collections import defaultdict
import pickle
import numpy as np
import pandas as pd
import pysam
from loguru import logger
from .tools import bedtoolsGetIntersect


def addGeneTag(readWithGeneInformation, inBamPath, outBamPath, geneIdTag):
    readWithGeneInformation = {
        x: y["gene_id"] for x, y in readWithGeneInformation.items()
    }
    with pysam.AlignmentFile(inBamPath) as inBam:
        with pysam.AlignmentFile(outBamPath, "wb", template=inBam) as outBam:
            for read in inBam:
                readGeneId = readWithGeneInformation.get(read.qname, "None")
                read.set_tag(geneIdTag, readGeneId)
                outBam.write(read)
    pysam.index(outBamPath)


NAMES = [
    "Chromosome",
    "Start",
    "End",
    "Name",
    "Score",
    "Strand",
    "ThickStart",
    "ThickEnd",
    "ItemRGB",
    "BlockCount",
    "BlockSizes",
    "BlockStarts",
    "geneChromosome",
    "geneStart",
    "geneEnd",
    "geneName",
    "geneScore",
    "geneStrand",
    "geneThickStart",
    "geneThickEnd",
    "geneItemRGB",
    "geneBlockCount",
    "geneBlockSizes",
    "geneBlockStarts",
    "cov",
]
USECOLS = [
    "Chromosome",
    "Start",
    "End",
    "Name",
    "Strand",
    "geneStart",
    "geneEnd",
    "geneName",
    "geneStrand",
    "geneBlockSizes",
    "geneBlockStarts",
    "cov",
]


def addGeneName(inBamPath, bedAnno, outfile, bedtoolsPath, outBamPath, geneIdTag):
    logger.info("Start get intersect between bam and bed12")
    intersectBuff = bedtoolsGetIntersect(inBamPath, bedAnno, bedtoolsPath)

    logger.info("Start read intersect result")
    df = pd.read_csv(intersectBuff, sep="\t", names=NAMES, usecols=USECOLS, header=None)
    logger.info("Read csv Done!")

    logger.info("Start find Splice Sites")
    df["geneBlockSizes"] = df["geneBlockSizes"].map(lambda x: np.fromstring(x, sep=","))
    df["geneBlockStarts"] = df["geneBlockStarts"].map(
        lambda x: np.fromstring(x, sep=",")
    )

    differentSr = df["geneBlockStarts"].map(len) < df["geneBlockSizes"].map(len)
    sizeDiffStartUmi = '\t'.join(list(df.loc[differentSr]['Name']))
    if sizeDiffStartUmi:
        logger.warning(f"block size counts different with block start counts, ignored: {sizeDiffStartUmi}")
        df = df.loc[~differentSr]
        # df['geneBlockStartCounts'] = df["geneBlockStarts"].map(len)
        # newGeneBlockSizesLs = []
        # for lineTp in df.itertuples():
        #     newGeneBlockSizesLs.append(lineTp.geneBlockStarts[:lineTp.geneBlockStartCounts])
        # df['geneBlockSizes'] = newGeneBlockSizesLs
        # del(df['geneBlockStartCounts'])

    df["five_ss"] = (
        df["geneStart"] + df["geneBlockSizes"] + df["geneBlockStarts"]
    ).map(lambda x: x[:-1])
    df["three_ss"] = (df["geneStart"] + df["geneBlockStarts"]).map(lambda x: x[1:])
    logger.info("Find Splice Sites Done!")

    logger.info("Main function")

    df["geneName"] = df["geneName"].str.split("\|").str[-1]
    specificLine = df["geneName"].str.endswith("_specific")
    specificDf = df.loc[specificLine]
    isoformDf = df.loc[~specificLine]

    results = defaultdict(
        lambda: {
            "cov": 0,
            "is_splicing_intermediates": False,
            "gene_len": 0,
            "gene_id": None,
            "isoform_mapping_ratio": {},
        }
    )
    for item in isoformDf.itertuples():
        # Determine if it is a splice intermediate
        if item.Strand == "+":
            is_splicing_intermediates = (abs(item.End - item.five_ss) <= 10).any()
        else:
            is_splicing_intermediates = (abs(item.Start - item.three_ss) <= 10).any()
        results[item.Name]["is_splicing_intermediates"] = (
            results[item.Name]["is_splicing_intermediates"] or is_splicing_intermediates
        )
        results[item.Name]["isoform_mapping_ratio"][item.geneName] = item.cov / int(
            item.geneBlockSizes.sum()
        )

        # Take the one with the greatest exon coverage as the gene annotation of this read
        if results[item.Name]["cov"] < item.cov:
            results[item.Name]["cov"] = item.cov
            results[item.Name]["gene_id"] = item.geneName
            results[item.Name]["gene_len"] = item.geneEnd - item.geneStart

    logger.info("Gene Assign Done!")

    # specificDf["geneNameTrue"] = specificDf["geneName"].str.split(".").str[0]
    # specificDf["geneAssignGene"] = specificDf.pipe(
    #     lambda x: x["Name"].map(lambda y: results[y]["gene_id"])
    # )
    # specificDf = specificDf.query("geneAssignGene == geneNameTrue")
    # specificDf["isoformName"] = specificDf["geneName"].str.split("_").str[0]
    # specificDf = specificDf.loc[:, ["Name", "cov", "isoformName"]]
    # specificDfGroup = specificDf.groupby("Name").apply(
    #     lambda z: {x: y for x, y in zip(z["isoformName"], z["cov"])}
    # )
    # specificDfGroup = specificDfGroup.to_dict()

    # for readName, geneBedInfo in results.items():
    #     isoformMappingRatio = geneBedInfo["isoform_mapping_ratio"]
    #     if len(isoformMappingRatio) == 1:
    #         isoformName = list(isoformMappingRatio.keys())[0]
    #     else:
    #         readSpecificDfGroup = specificDfGroup.get(readName, {})
    #         readSpecificDfGroup = {
    #             x: readSpecificDfGroup.get(x, 0) for x in isoformMappingRatio.keys()
    #         }
    #         maxIsoformCoverageLength = max(readSpecificDfGroup.values())
    #         cutLength = maxIsoformCoverageLength - 15
    #         putativeIsoforms = [
    #             x for x, y in readSpecificDfGroup.items() if y > cutLength
    #         ]
    #         if len(putativeIsoforms) == 1:
    #             isoformName = putativeIsoforms[0]
    #         else:
    #             putativeIsoforms.sort(key=lambda x: isoformMappingRatio[x])
    #             putativeIsoformMappingRatio = np.array(
    #                 [isoformMappingRatio[x] for x in putativeIsoforms]
    #             )
    #             if (
    #                 putativeIsoformMappingRatio[-1] - putativeIsoformMappingRatio[-2]
    #                 > 0.1
    #             ):
    #                 isoformName = putativeIsoforms[-1]
    #             else:
    #                 isoformName = geneBedInfo["gene_id"] + ".N"
    #     results[readName]["isoform_id"] = isoformName

    logger.info("Main function Done!")

    with open(outfile, "wb") as o:
        pickle.dump(dict(results), o)

    addGeneTag(dict(results), inBamPath, outBamPath, geneIdTag)
