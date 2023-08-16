"""
@Author       : windz
@Date         : 2020-05-24 19:07:01
LastEditTime: 2020-10-14 20:22:53
@Description  :
"""


from collections import Counter
import os
from loguru import logger
import numpy as np
import pysam
import click
import pyranges as pr
import pyfastx
from concurrent.futures import ProcessPoolExecutor
import sh
from tempfile import TemporaryDirectory
import pandas as pd


def getSeqFromBedFile(df_bed, path_genome, path_bedtools):
    dir_temp = TemporaryDirectory()
    df_bed.to_csv(dir_temp.name + "t.bed", sep="\t", header=None, index=None)
    sh.Command(path_bedtools).getfasta(
        fi=path_genome,
        bed=dir_temp.name + "t.bed",
        bedOut=True,
        s=True,
        _long_sep=" ",
        _long_prefix="-",
        _out=dir_temp.name + "t.withFa.bed",
    )
    df_bedSeq = pd.read_table(
        dir_temp.name + "t.withFa.bed", names=[*df_bed.columns, "Seq"]
    )
    return df_bedSeq


def get_entropy(site_counter, total_count):
    """
    calculate polyA site entropy of this PAC
    """
    entropy = 0
    for site in site_counter:
        p = site_counter[site] / total_count
        entropy += -p * np.log2(p)
    return entropy


STRAND_TO_BOOL = {"-": True, "+": False}


def get_three_end(infile, gene_id, gene_model):
    """
    obtain position of the polyadenylated read 3’end
    """
    chrom, start, end, _, strand = gene_model.loc[gene_id, :].values
    strand = STRAND_TO_BOOL[strand]
    read_count = 0
    polya_sites = []
    with pysam.AlignmentFile(infile, "rb") as inbam:
        for read in inbam.fetch(chrom, start, end):
            # keep the same strand reads
            if strand is read.is_reverse:
                polya_len = read.get_tag("pa")
                read_gene_id = read.get_tag("gi")
                if polya_len >= 15 and read_gene_id in {"None", gene_id}:
                    read_count += 1
                    if not read.is_reverse:
                        polya_sites.append(read.reference_end)
                    else:
                        polya_sites.append(read.reference_start * -1)

    # choose these gene with high abundance
    if read_count < 10:
        return None

    polya_sites.sort()

    # merge site within 24nt
    total_site_count = len(polya_sites)
    pac_list = None
    for polya_site in polya_sites:
        if pac_list is None:
            pac_list = [[polya_site]]
        else:
            if polya_site - pac_list[-1][-1] <= 24:
                pac_list[-1].append(polya_site)
            else:
                pac_list.append([polya_site])

    major_cluster_site_count = 0
    summit = []

    polya_cluster = []
    polya_cluster_summit = []

    polya_cluster_major = ""
    polya_cluster_summit_major = ""

    polya_cluster_last = ""
    polya_cluster_summit_last = ""

    *_, strand = gene_model.loc[gene_id, :].values
    n = 0
    for pac in pac_list:
        # remove PAC with a read counts lower than 3
        # remove PAC with a read counts lower than 1% of total reads
        if len(pac) < 3 or len(pac) / total_site_count < 0.1:
            continue

        site_counter = Counter(pac)
        site_counter_most = site_counter.most_common(3)

        if site_counter_most[0][1] >= 3:
            n += 1
            start = min(abs(pac[0]), abs(pac[-1]))
            end = max(abs(pac[0]), abs(pac[-1]))
            entropy = get_entropy(site_counter, len(pac))
            polya_cluster.append(
                f"{chrom}\t{start-1}\t{end}\t{gene_id}_{n}\t{len(pac)}\t{strand}\t{entropy:.3f}\n"
            )

            # calculate the summit of cluster
            summit = abs(site_counter_most[0][0])
            # last column is summit/pac, and higher value indicates a significant summit
            polya_cluster_summit.append(
                f"{chrom}\t{summit-1}\t{summit}\t{gene_id}_{n}\t{site_counter_most[0][1]}\t{strand}\t{site_counter_most[0][1]/len(pac):.3f}\n"
            )
            if major_cluster_site_count < len(pac):
                major_cluster_site_count = len(pac)
                polya_cluster_major = polya_cluster[-1]
                polya_cluster_summit_major = polya_cluster_summit[-1]

    if len(polya_cluster) == 0:
        return None

    polya_cluster_last = polya_cluster[-1]
    polya_cluster_summit_last = polya_cluster_summit[-1]

    return (
        polya_cluster,
        polya_cluster_summit,
        polya_cluster_major,
        polya_cluster_summit_major,
        polya_cluster_last,
        polya_cluster_summit_last,
    )

def _fc(line, length=20):
    dt_line = line.to_dict()
    if dt_line['Strand'] == '+':
        dt_line['Start'] = dt_line['End'] - length
    else:
        dt_line['Start'] =  dt_line['Start'] + 1 
        dt_line['End'] = dt_line['Start'] + length
    sr_result = pd.Series(dt_line)
    return sr_result

def filterPAC(fastaPath, bedPath, bedSummitPath, fillterPolyASitePath, path_bedtools):
    df_pac = pr.read_bed(bedSummitPath, as_df=True)
    df_pacTerminal = df_pac.apply(_fc, axis=1, length=20)
    df_pacTerminalSeq = getSeqFromBedFile(df_pacTerminal, fastaPath, path_bedtools)
    df_pacTerminalSeqFiltered = df_pacTerminalSeq.loc[~(df_pacTerminalSeq['Seq'].str[-3:] == 'AAA')]
    # df_pacTerminalSeqFiltered = df_pacTerminalSeqFiltered.assign(
    #     Gene=lambda df: df["Name"].str.split("_").str[:-1].str.join("_")
    # )
    # ls_usedGene = df_pacTerminalSeqFiltered.value_counts('Gene').pipe(lambda sr:sr[sr > 1]).index.tolist()
    # df_pacTerminalSeqFiltered = df_pacTerminalSeqFiltered.query("Gene in @ls_usedGene")
    ls_usePac = df_pacTerminalSeqFiltered['Name'].to_list()

    # old version
    # genomeFa = pyfastx.Fasta(fastaPath)
    # polyAClusterBed = pr.read_bed(bedSummitPath, True)
    # polyAClusterBed["Chromosome"] = polyAClusterBed["Chromosome"].astype(str)
    # polyAClusterBed["seq"] = polyAClusterBed.apply(
    #     lambda x: genomeFa[x["Chromosome"]][x["Start"] - 10 : x["End"] + 10].seq, axis=1
    # )
    # polyAClusterBed["seqLength"] = polyAClusterBed["seq"].map(len)
    # polyAClusterBed["Ratio"] = (
    #     polyAClusterBed.apply(
    #         lambda x: x["seq"].count("A")
    #         if x["Strand"] == "+"
    #         else x["seq"].count("T"),
    #         axis=1,
    #     )
    #     / polyAClusterBed["seqLength"]
    # )
    # usePolyASite = polyAClusterBed.query("Ratio <= 0.5")["Name"]
    polyAClusterRawRangeBed = pr.read_bed(bedPath, True)
    polyAClusterRawRangeBed["Chromosome"] = polyAClusterRawRangeBed[
        "Chromosome"
    ].astype(str)
    polyAClusterPassedRangeBed = polyAClusterRawRangeBed.query("Name in @ls_usePac")
    polyAClusterPassedRangeBed.to_csv(
        fillterPolyASitePath, sep="\t", header=None, index=None
    )


def polyAClusterDetected(
    fastaPath, infile, gene_bed, out_suffix, threads, is_bed12=False, path_bedtools='bedtools'
):
    try:
        os.mkdir(out_suffix)
    except:
        logger.warning(f"{out_suffix} existed!")
    gene_model = pr.read_bed(gene_bed, as_df=True)
    gene_model["Chromosome"] = gene_model["Chromosome"].astype(str)
    if is_bed12:
        # transfer bed12 to bed
        gene_model = gene_model.assign(
            Gene=lambda df: df["Name"].str.split("\|").str[-1],
            GeneBiotype=lambda df: df["Name"].str.split("\|").str[-2],
        )
        gene_model = gene_model.groupby("Gene").agg(
            {
                "Chromosome": lambda x: x.iat[0],
                "Start": "min",
                "End": "max",
                "Strand": lambda x: x.iat[0],
                "GeneBiotype": lambda x: x.iat[0],
            }
        )
        gene_model = gene_model.reset_index().pipe(
            lambda df: df.assign(
                Name=df["Gene"] + "|" + df["GeneBiotype"] + "|" + df["Gene"],
                Score=0,
                ThickStart=df["Start"],
                ThickEnd=df["End"],
                ItemRGB=0.0,
                BlockCount=1,
                BlockSizes=(df["End"] - df["Start"]).astype(str) + ",",
                BlockStarts="0,",
            )
        )
        gene_model = (
            gene_model.reindex(
                columns=[
                    "Chromosome",
                    "Start",
                    "End",
                    "Gene",
                    "Score",
                    "Strand",
                ]
            )
            .sort_values(["Chromosome", "Start"])
            .rename(columns={"Gene": "Name"})
        )

    gene_model = gene_model.set_index(["Name"])

    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        for gene_id in gene_model.index:
            if gene_model.at[gene_id, "Chromosome"] not in {"Mt", "Pt"}:
                results.append(e.submit(get_three_end, infile, gene_id, gene_model))

    o_polya_cluster = open(f"{out_suffix}polya_cluster.bed", "w")
    o_polya_cluster_summit = open(f"{out_suffix}polya_cluster.summit.bed", "w")
    o_major_polya_cluster = open(f"{out_suffix}major_polya_cluster.bed", "w")
    o_major_polya_cluster_summit = open(
        f"{out_suffix}major_polya_cluster_summit.bed", "w"
    )
    o_last_polya_cluster = open(f"{out_suffix}last_polya_cluster.bed", "w")
    o_last_polya_cluster_summit = open(
        f"{out_suffix}last_polya_cluster_summit.bed", "w"
    )

    for res in results:
        result = res.result()
        if result is not None:
            (
                polya_cluster,
                polya_cluster_summit,
                polya_cluster_major,
                polya_cluster_summit_major,
                polya_cluster_last,
                polya_cluster_summit_last,
            ) = result
            for item in polya_cluster:
                o_polya_cluster.write(item)
            for item in polya_cluster_summit:
                o_polya_cluster_summit.write(item)

            o_major_polya_cluster.write(polya_cluster_major)
            o_major_polya_cluster_summit.write(polya_cluster_summit_major)
            o_last_polya_cluster.write(polya_cluster_last)
            o_last_polya_cluster_summit.write(polya_cluster_summit_last)

    o_polya_cluster.close()
    o_polya_cluster_summit.close()
    o_major_polya_cluster.close()
    o_major_polya_cluster_summit.close()
    o_last_polya_cluster.close()
    o_last_polya_cluster_summit.close()

    filterPAC(
        fastaPath,
        f"{out_suffix}polya_cluster.bed",
        f"{out_suffix}polya_cluster.summit.bed",
        f"{out_suffix}polya_cluster.filtered.bed",
        path_bedtools
    )
