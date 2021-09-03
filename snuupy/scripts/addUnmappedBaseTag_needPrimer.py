import pysam
from Bio import Align
from tqdm import tqdm
from typing import List, Optional
from concurrent.futures import ProcessPoolExecutor
from more_itertools import chunked
from .tools import getAntisense, writeFasta


def get_putative(
    read: pysam.AlignedSegment, seq_primer5="CTACACGACGCTCTTCCGATCT"
) -> Optional[pysam.AlignedSegment]:
    if (read.query_alignment_start > 150) or (read.query_length - read.query_alignment_end > 150):
        return None

    seq = read.query_sequence
    seq_leftClip = seq[: read.query_alignment_start]
    seq_rightClip = seq[read.query_alignment_end:]
    if not seq_leftClip:
        seq_leftClip = "A"
    if not seq_rightClip:
        seq_rightClip = "A"

    aligner = Align.PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.gap_score = -2
    aligner.query_gap_score = -1
    aligner.target_end_gap_score = 0
    seq_primer5 = seq_primer5
    seq_primer5_rev = getAntisense(seq_primer5)

    primer5FwdMapping = aligner.align(seq_primer5, seq_leftClip)
    primer5RevMapping = aligner.align(seq_primer5_rev, seq_rightClip)

    primer5FwdScore = primer5FwdMapping.score
    primer5RevScore = primer5RevMapping.score
    ls_score = [primer5FwdScore, primer5RevScore]
    ls_unmapBase = [
        (len(seq_primer5) - primer5FwdScore) / 2,
        (len(seq_primer5_rev) - primer5RevScore) / 2,
    ]

    if (max(ls_unmapBase) < 5) | (min(ls_unmapBase) > 5):
        return None

    _index = ls_unmapBase.index(min(ls_unmapBase))
    mapping = [primer5FwdMapping, primer5RevMapping][_index]
    seq_clip = [seq_leftClip, seq_rightClip][_index]

    for x in mapping:  # get align pos
        break
    if _index == 1:
        seq_putative = seq_clip[: x.aligned[1][0][0]]
        cat = 'mapping_R5'
    else:
        seq_putative = seq_clip[x.aligned[1][-1][-1] :]
        cat = 'mapping_F5'
    if not seq_putative:
        seq_putative = 'A'
    read.set_tag("PS", seq_putative)
    read.set_tag("PC", cat)

    return read


def processChunkBam(bamFile: List[pysam.AlignedSegment]) -> List[pysam.AlignedSegment]:
    ls_reads = []
    for i, read in enumerate(bamFile):
        if (read.query_alignment_start > 150) or (
            read.query_length - read.query_alignment_end > 150
        ):
            continue
        read_afterProcess = get_putative(read)
        if read_afterProcess:
            ls_reads.append(read)
    return ls_reads


def main(path_bamFile, path_outFile):
    bamFile = pysam.AlignmentFile(path_bamFile, "rb")
    bamFileOut = pysam.AlignmentFile(path_outFile, "wbu", template=bamFile)
    counts_bam = bamFile.mapped + bamFile.unmapped


    it_bam = chunked(bamFile, 1000)
    for _ls_bam in tqdm(it_bam, 'extract unmaped seq', total=(counts_bam // (1000) + 1)):
        ls_processedReads = map(get_putative, _ls_bam)
        ls_processedReads = [x for x in ls_processedReads if x]
        for read in ls_processedReads:
            bamFileOut.write(read)

    # it_bam = chunked(bamFile, threads * 10000)
    # for _ls_bam in tqdm(it_bam, 'extract unmaped seq', total=(counts_bam // (threads * 10000) + 1)):
    #     with ProcessPoolExecutor(threads) as mtp:
    #         ls_processedReads = mtp.map(get_putative, _ls_bam)
    #     ls_processedReads = [x for x in ls_processedReads if x]
    #     for read in ls_processedReads:
    #         bamFileOut.write(read)
    bamFileOut.close()
    pysam.index(path_outFile)

