import pandas as pd
import numpy as np
import sh
from concurrent.futures import ProcessPoolExecutor as multiP
from loguru import logger
from io import StringIO
import pyranges as pr
from tempfile import NamedTemporaryFile
import click
import pickle
import portion as P
from collections import defaultdict
from tqdm import tqdm


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
    "BlockSizes",
    "BlockStarts",
    "geneStart",
    "geneEnd",
    "geneName",
    "geneBlockCount",
    "geneBlockSizes",
    "geneBlockStarts",
    "cov",
]


def bedtoolsGetIntersect(inBam, bedAnno, bedtoolsPath, bed12=True):
    intersectBuff = StringIO()
    if bed12:
        sh.Command(bedtoolsPath).intersect(
            "-abam",
            inBam,
            "-b",
            bedAnno,
            "-wo",
            "-s",
            "-split",
            "-bed",
            _out=intersectBuff,
        )
        intersectBuff.seek(0)
        bedFile = pd.read_table(
            intersectBuff, header=None, names=NAMES, usecols=USECOLS
        )
        bedFile["BlockStarts"] = (
            bedFile["BlockStarts"]
            .astype(str)
            .map(lambda x: np.fromstring(x, sep=",", dtype=int))
        )
        bedFile["BlockSizes"] = (
            bedFile["BlockSizes"]
            .astype(str)
            .map(lambda x: np.fromstring(x, sep=",", dtype=int))
        )
        bedFile["geneBlockSizes"] = (
            bedFile["geneBlockSizes"]
            .astype(str)
            .map(lambda x: np.fromstring(x, sep=",", dtype=int))
        )
        bedFile["geneBlockStarts"] = (
            bedFile["geneBlockStarts"]
            .astype(str)
            .map(lambda x: np.fromstring(x, sep=",", dtype=int))
        )
    else:
        sh.Command(bedtoolsPath).intersect(
            "-abam", inBam, "-b", bedAnno, "-wo", "-s", "-bed", "-split", _out=intersectBuff
        )
        intersectBuff.seek(0)
        bedFile = pd.read_table(intersectBuff, header=None)
        bedFile = bedFile[[3, 13, 14, 15, 18]].rename(
            columns={3: "Name", 13: "geneStart", 14: "geneEnd", 15: "geneName", 18: "cov"}
        )
    return bedFile

def getLongestIsoform(path_bed, path_tempBed):
    df_bed = pr.read_bed(path_bed, as_df=True)
    df_bed["IsoformLength"] = df_bed["BlockSizes"].map(
        lambda z: sum([int(x) for x in z.split(",")[:-1]])
    ) - df_bed.eval("`ThickStart` - Start + End - `ThickEnd`")
    df_bed["Gene"] = df_bed["Name"].str.split("\|").str[-1]
    df_bed = (
        df_bed.sort_values("IsoformLength", ascending=False)
        .drop_duplicates("Gene")
        .sort_values(["Chromosome", "Start"])
        .drop(columns=["IsoformLength", "Gene"])
    )
    df_bed.to_csv(path_tempBed, sep="\t", header=None, index=None)
    return path_tempBed

def generateGeneBed(path_bed, path_tempBed):
    df_bed = pr.read_bed(path_bed, as_df=True)
    df_bed = df_bed.pipe(
        lambda df: df.assign(
            ItemRGB=0.0,
            BlockCount=1,
            BlockSizes=(df["End"] - df["Start"]).astype(str) + ",",
            BlockStarts="0,",
        )
    )
    df_bed.to_csv(path_tempBed, sep="\t", header=None, index=None)
    return path_tempBed

def _getIntrons(line, bed12):
    if int(line.BlockCount) <= 1:
        return None

    ls_tuple = []
    for start, length in zip(
        line.BlockStarts.split(",")[:-1], line.BlockSizes.split(",")[:-1]
    ):
        start = int(start)
        length = int(length)
        ls_tuple.append(P.closedopen(start, start + length))
    iv_exon = P.Interval(*ls_tuple)
    iv_gene = P.closedopen(0, int(line.End) - int(line.Start))
    iv_intron = iv_gene - iv_exon
    ls_intron = list(iv_intron)
    if not bed12:
        if line.Strand == "-":
            ls_intron = ls_intron[::-1]

        ls_intronFeature = []
        for intronNum, iv_singleIntron in zip(range(1, 1 + len(ls_intron)), ls_intron):
            ls_intronFeature.append(
                [
                    f"{line.Name}_intron{intronNum}",
                    line.Chromosome,
                    line.Start + iv_singleIntron.lower,
                    line.Start + iv_singleIntron.upper,
                    line.Strand,
                ]
            )
        return ls_intronFeature
    else:
        Start = line.Start + ls_intron[0].lower
        BlockStarts = (
            ",".join([str(x.lower - ls_intron[0].lower) for x in ls_intron]) + ","
        )
        BlockSizes = ",".join([str(x.upper - x.lower) for x in ls_intron]) + ","
        BlockCount = len(ls_intron)
        End = line.Start + ls_intron[-1].upper
        sr_intron = pd.Series(
            dict(
                Chromosome=line.Chromosome,
                Start=Start,
                End=End,
                Name=line.Name,
                Score=line.Score,
                Strand=line.Strand,
                ThickStart=Start,
                ThickEnd=End,
                ItemRGB=line.ItemRGB,
                BlockCounts=BlockCount,
                BlockSizes=BlockSizes,
                BlockStarts=BlockStarts,
            )
        )
        return sr_intron


def generateIntrons(path_bed, path_tempBed, bed12=True) -> pd.DataFrame:
    df_bed = pr.read_bed(path_bed, as_df=True)
    ls_introns = [
        _getIntrons(x, bed12)
        for x in tqdm(
            df_bed.query("BlockCount > 1").itertuples(),
            total=len(df_bed.query("BlockCount > 1")),
        )
    ]
    if not bed12:
        df_introns = pd.DataFrame(
            [y for x in ls_introns for y in x],
            columns=["Name", "Chromosome", "Start", "End", "Strand"],
        ).assign(Score = 0)
    
        df_introns = df_introns[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]
        df_introns["Chromosome"] = df_introns["Chromosome"].astype(str)
        df_introns = df_introns.sort_values(["Chromosome", "Start"])
    else:
        df_introns = pd.DataFrame(
            ls_introns,
        )
        df_introns["Chromosome"] = df_introns["Chromosome"].astype(str)
        df_introns = df_introns.sort_values(["Chromosome", "Start"])
    df_introns.to_csv(path_tempBed, sep="\t", header=None, index=None)
    return path_tempBed


def _getExons(line):
    ls_tuple = []
    for start, length in zip(
        line.BlockStarts.split(",")[:-1], line.BlockSizes.split(",")[:-1]
    ):
        start = int(start)
        length = int(length)
        ls_tuple.append(P.closedopen(start, start + length))
    iv_exon = P.Interval(*ls_tuple)
    ls_exon = list(iv_exon)
    if line.Strand == "-":
        ls_exon = ls_exon[::-1]

    ls_exonFeature = []
    for exonNum, iv_singleExon in zip(range(1, 1 + len(ls_exon)), ls_exon):
        ls_exonFeature.append(
            [
                f"{line.Name}_exon{exonNum}",
                line.Chromosome,
                line.Start + iv_singleExon.lower,
                line.Start + iv_singleExon.upper,
                line.Strand,
            ]
        )
    return ls_exonFeature


def generateExons(path_bed, path_tempBed) -> pd.DataFrame:
    df_bed = pr.read_bed(path_bed, as_df=True)
    ls_exons = [
        _getExons(x)
        for x in tqdm(
            df_bed.itertuples(),
            total=len(df_bed),
        )
    ]
    df_exons = pd.DataFrame(
        [y for x in ls_exons for y in x],
        columns=["Name", "Chromosome", "Start", "End", "Strand"],
    ).assign(Score = 0)
    
    df_exons = df_exons[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]
    
    df_exons['Chromosome'] = df_exons['Chromosome'].astype(str)
    df_exons = df_exons.sort_values(['Chromosome', 'Start'])
    df_exons.to_csv(path_tempBed, sep="\t", header=None, index=None)
    return path_tempBed 


def getSpliceInfo(
    INBAM_PATH, BED_REPRE_ANNO, GENE_NAME_INFO, OUT_PATH, NEED_RATIO, bedtoolsPath, intronCutoff=0.1, INCLUDE_INTRON=True
):
    if GENE_NAME_INFO:
        logger.warning(f"parameter GENE_NAME_INFO is deprecated")
    if NEED_RATIO:
        logger.warning(f"parameter NEED_RATIO is deprecated")

    path_tempBed = NamedTemporaryFile("w", suffix=".bed")
    BED_REPRE_ANNO = getLongestIsoform(BED_REPRE_ANNO, path_tempBed.name)
    df_bed = pr.read_bed(BED_REPRE_ANNO, True).set_index('Name')
    dt_geneExonCounts = df_bed['BlockCount'].to_dict()
    dt_geneStrand = df_bed['Strand'].to_dict()

    if INCLUDE_INTRON:
        path_tempBedGene = NamedTemporaryFile("w", suffix=".bed")
        BED_GENE = generateGeneBed(BED_REPRE_ANNO, path_tempBedGene.name)
    else:
        BED_GENE = BED_REPRE_ANNO

    df_geneIntersect = bedtoolsGetIntersect(INBAM_PATH, BED_GENE, bedtoolsPath)

    dt_barcodeWithGene = (
        df_geneIntersect.sort_values(["Name", "cov"], ascending=False)
        .drop_duplicates(["Name"])
        .set_index("Name")["geneName"]
        .to_dict()
    )

    path_tempBedExon = NamedTemporaryFile("w", suffix=".bed")
    BED_EXON = generateExons(BED_REPRE_ANNO, path_tempBedExon.name)
    df_exonIntersect = bedtoolsGetIntersect(INBAM_PATH, BED_EXON, bedtoolsPath, bed12=False)

    df_exonIntersect = df_exonIntersect.assign(
        trueName=lambda df: df["Name"].map(dt_barcodeWithGene),
        exonName=lambda df: df["geneName"],
    ).assign(geneName=lambda df: df['geneName'].str.split('_exon').str[0])

    df_exonIntersect = df_exonIntersect.query("geneName == trueName")
    df_exonIntersect = df_exonIntersect.assign(exonNum = lambda df:df['exonName'].str.split('_exon').str[-1].astype(int))
    df_exonCov = df_exonIntersect.groupby('Name')['exonNum'].agg(['min', 'max'])
    dt_exonCov = df_exonCov.to_dict('index')

    path_tempBedIntron = NamedTemporaryFile("w", suffix=".bed")
    BED_INTRON = generateIntrons(BED_REPRE_ANNO, path_tempBedIntron.name, bed12=False)
    df_intronIntersect = bedtoolsGetIntersect(INBAM_PATH, BED_INTRON, bedtoolsPath, bed12=False)

    df_intronIntersect = df_intronIntersect.assign(
        trueName=lambda df: df["Name"].map(dt_barcodeWithGene),
        intronName=lambda df: df["geneName"],
    ).assign(geneName=lambda df: df['geneName'].str.split('_intron').str[0])

    df_intronIntersect = df_intronIntersect.query("geneName == trueName")
    df_intronIntersect = df_intronIntersect.assign(
        intronNum=lambda df: df["intronName"].str.split("_intron").str[-1].astype(int),
        intronRatio = lambda df:df.eval("cov / (geneEnd - geneStart)")
    )

    dt_intronCovInf = defaultdict(lambda :{})
    for line in tqdm(df_intronIntersect.itertuples(), total = len(df_intronIntersect)):
        dt_UmiIntron = dt_intronCovInf[line.Name]
        dt_UmiIntron[line.intronNum] = line.intronRatio
    dt_intronCovInf = dict(dt_intronCovInf)

    ls_results = []
    intronCutoff = 0.1
    for umi, gene in tqdm(dt_barcodeWithGene.items()):
        exonCounts = dt_geneExonCounts[gene]
        strand = dt_geneStrand[gene]
        # record exon info
        if umi in dt_exonCov:
            exonOverlap = dt_exonCov[umi]
            exonMin = exonOverlap["min"] - 1  # 1-base to 0-base
            exonMax = exonOverlap["max"]
            exonOverlapInfo = ",".join([str(x) for x in range(exonMin, exonMax)])
        else:
            exonOverlapInfo = ""

        # record intron info
        if umi in dt_intronCovInf:
            dt_intronOverlapRatioInfo = {x - 1:y for x,y in dt_intronCovInf[umi].items()} # 1-base to 0-base
        else:
            dt_intronOverlapRatioInfo = {}
        if exonOverlapInfo:
            intronMin = exonMin
            intronMax = exonMax - 1
            _dt = {x: 0.0 for x in range(intronMin, intronMax)}
            _dt.update(dt_intronOverlapRatioInfo)
            dt_intronOverlapRatioInfo = _dt
        intronOverlapRatio = ",".join(
            [f"{x}:{y}" for x, y in dt_intronOverlapRatioInfo.items()]
        )
        intronRetained = ",".join(
            [str(x) for x, y in dt_intronOverlapRatioInfo.items() if y > intronCutoff]
        )
    #     sr_result = pd.Series(
    #         dict(
    #             Name=umi,
    #             Strand=strand,
    #             GeneExonCounts=exonCounts,
    #             ExonOverlapInfo=exonOverlapInfo,
    #             IntronOverlapInfo=intronRetained,
    #             intronOverlapRatioInfo=intronOverlapRatio,
    #             geneId=gene.split("|")[-1],
    #         )
    #     ) # too slow
        ls_results.append([umi, strand, exonCounts, exonOverlapInfo, intronRetained, intronOverlapRatio, gene.split("|")[-1]])

    df_results = pd.DataFrame(
        ls_results,
        columns=[
            "Name",
            "Strand",
            "GeneExonCounts",
            "ExonOverlapInfo",
            "IntronOverlapInfo",
            "intronOverlapRatioInfo",
            "geneId",
        ],
    ).sort_values(['geneId', 'Name'])
    df_results.to_csv(OUT_PATH, sep="\t", index=False)