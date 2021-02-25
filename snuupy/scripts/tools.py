import pysam
import os
import sh
from io import StringIO
from collections import namedtuple
import pandas as pd
import numpy as np
import glob
from ont_fast5_api.fast5_interface import get_fast5_file
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, zscore
from io import StringIO
from concurrent.futures import ProcessPoolExecutor as Mtp
from loguru import logger
import pickle
import lmdb
import anndata


def getBlock(read, intron):
    """
    @description: get read block
    @param
        read:{pysam.read}
        intron:{pysam.intron}
    @return:
        [(start, end),(start, end)]
    """
    block = []
    preStart, lastEnd = read.reference_start, read.reference_end
    for intronStart, intronEnd in intron:
        block.append((preStart, intronStart))
        preStart = intronEnd
    block.append((preStart, lastEnd))
    return block


def isOne(n, i):
    """
    @description: Binary
    @param
        n:{int}
        i:{int}
    @return:
        bool
    """
    return (n & (1 << i)) != 0


def creatUnexistedDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


class sequence:
    def __init__(self):
        old_chars = "ACGT"
        replace_chars = "TGCA"
        self.tab = str.maketrans(old_chars, replace_chars)

    def original(self, seq):
        return seq

    def complement(self, seq):
        return seq.translate(self.tab)

    def reverse(self, seq):
        return seq[::-1]

    def reverseComplement(self, seq):
        return seq.translate(self.tab)[::-1]


def readFasta(path):
    """
    @description: read fasta
    @param {type} fasta file path
    @return: read generator
    """
    FastaRead = namedtuple("FastaRead", ["name", "seq"])

    def _readFasta(path):
        with open(path, "r") as fh:
            i = 0
            while True:
                lineContent = fh.readline().strip()
                if lineContent == "":
                    break
                if lineContent.startswith(">"):
                    i += 1
                    if i == 1:
                        readName = lineContent[1:].split(" ")[0]
                        readSeq = ""
                    else:
                        read = FastaRead(name=readName, seq=readSeq)
                        yield read
                        readName = lineContent[1:].split(" ")[0]
                        readSeq = ""
                else:
                    readSeq += lineContent
            read = FastaRead(name=readName, seq=readSeq)
            yield read

    return _readFasta(path)


def getAntisense(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    transMap = str.maketrans(old_chars, replace_chars)
    return seq.translate(transMap)[::-1]


class Fasta:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __getitem__(self, key):
        sliceFasta = Fasta(self.name, self.seq[key])
        return sliceFasta

    def __str__(self):
        return f"{self.name}:\n{self.seq}"

    __repr__ = __str__

    def getAnti(self):
        return Fasta(self.name, getAntisense(self.seq))


class FastaContent:
    """
    fasta content
    """

    def __init__(self, path, useIndex=False):
        self.path = path
        self.useIndex = useIndex
        if self.useIndex:
            self.indexPath = f"{self.path}_lmdb/"
            if os.path.exists(self.indexPath):
                pass
            else:
                self.buildIndex()
            self.lmdbEnv = lmdb.Environment(self.indexPath)
            self.lmdbTxn = self.lmdbEnv.begin()
            self.keys = pickle.loads(self.lmdbTxn.get("thisFileIndex".encode()))
            self.__keySt = set(self.keys)

    def __len__(self):
        if not self.useIndex:
            logger.error(f"NO Index mode!")
            0 / 0
        else:
            return len(self.keys)

    def __getitem__(self, readName):
        if not self.useIndex:
            logger.error(f"NO Index mode!")
            0 / 0
        elif readName not in self.__keySt:
            logger.error(f"NOT Found Read name!")
            0 / 0
        else:
            return self.__readFastaReadFromLmdb(readName)

    def buildIndex(self):
        """
        build lmdb database
        """

        def __writeFastaReadToLmdb(lmdbTxn, fastaRead):
            lmdbTxn.put(
                key=f"{fastaRead.name}_name".encode(), value=fastaRead.name.encode()
            )
            lmdbTxn.put(
                key=f"{fastaRead.name}_seq".encode(), value=fastaRead.seq.encode()
            )

        lmdbEnv = lmdb.open(self.indexPath, map_size=1099511627776)
        lmdbTxn = lmdbEnv.begin(write=True)
        readNameLs = []
        for i, fastaRead in enumerate(self.__readFasta()):
            readNameLs.append(fastaRead.name)
            __writeFastaReadToLmdb(lmdbTxn, fastaRead)
            if i % 1e5 == 0:
                logger.info(f"{i} reads processed")
        readNameMergedPk = pickle.dumps(readNameLs)
        lmdbTxn.put(key="thisFileIndex".encode(), value=readNameMergedPk)
        lmdbTxn.commit()
        lmdbEnv.close()

    def __readFasta(self):
        """
        read fasta file
        """
        with open(self.path, "r") as fh:
            i = 0
            while True:
                lineContent = fh.readline().strip()
                if lineContent == "":
                    break
                if lineContent.startswith(">"):
                    i += 1
                    if i == 1:
                        readName = lineContent[1:].split(" ")[0]
                        readSeq = ""
                    else:
                        read = Fasta(name=readName, seq=readSeq)
                        yield read
                        readName = lineContent[1:].split(" ")[0]
                        readSeq = ""
                else:
                    readSeq += lineContent
            read = Fasta(name=readName, seq=readSeq)
            yield read

    def __readFastaReadFromLmdb(self, readName):
        """
        read fasta form lmdb
        """
        read = Fasta(
            name=self.lmdbTxn.get(f"{readName}_name".encode()).decode(),
            seq=self.lmdbTxn.get(f"{readName}_seq".encode()).decode(),
        )
        return read

    def __readLmdb(self):
        """
        read lmdb file
        """
        for readName in self.keys:
            yield self.__readFastaReadFromLmdb(readName)

    def iter(self):
        """
        iter
        """
        if self.useIndex:
            return self.__readLmdb()
        else:
            return self.__readFasta()

    def close(self):
        self.lmdbEnv.close()


class Jinterval:
    def __init__(self, lower, upper, overlapLimit=0.5):
        self.lower, self.upper = lower, upper
        self.interval = [lower, upper]
        self.overlapLimit = overlapLimit

    def __repr__(self):
        return f"Jinterval{self.interval}"

    def __str__(self):
        return f"Jinterval{self.interval}"

    def __and__(self, otherInterval):
        minn = max(self.lower, otherInterval.lower)
        maxn = min(self.upper, otherInterval.upper)
        if (maxn - minn) / (self.upper - self.lower) > self.overlapLimit:
            return [minn, maxn]
        else:
            return False

    def getOverlapRatio(self, otherInterval):
        minn = max(self.lower, otherInterval.lower)
        maxn = min(self.upper, otherInterval.upper)
        return max((maxn - minn) / (self.upper - self.lower), 0)


def bedtoolsGetIntersect(inBam, bedAnno, bedtoolsPath):
    intersectBuff = StringIO()
    sh.Command(bedtoolsPath).intersect(
        "-abam", inBam, "-b", bedAnno, "-wo", "-s", "-split", "-bed", _out=intersectBuff
    )
    intersectBuff.seek(0)
    return intersectBuff


def extract_read_data(fast5_filepath, read_id):
    """
    @description:
        It can handle fast5 basecalled with flip flop model.

    @param
        fast5_filepath.
        read_id.

    @return:
        raw_data
        event_data
        fastq
        start
        stride
        samples_per_nt
    """
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        read = f5.get_read(read_id)
        # compute event length vector
        model_type = read.get_analysis_attributes("Basecall_1D_000")["model_type"]
        if model_type == "flipflop":
            # get the data
            raw_data = read.get_raw_data()
            fastq = read.get_analysis_dataset(
                group_name="Basecall_1D_000/BaseCalled_template", dataset_name="Fastq"
            )
            fastq = fastq.split("\n")[1]
            start = read.get_summary_data("Segmentation_000")["segmentation"][
                "first_sample_template"
            ]
            move = read.get_analysis_dataset(
                group_name="Basecall_1D_000/BaseCalled_template", dataset_name="Move"
            )
            stride = read.get_summary_data("Basecall_1D_000")["basecall_1d_template"][
                "block_stride"
            ]
            start_col = np.arange(start, start + stride * (len(move) - 1) + 1, stride)
            event_data = pd.DataFrame(
                {"move": move, "start": start_col, "move_cumsum": np.cumsum(move)}
            )
            event_data["model_state"] = event_data["move_cumsum"].map(
                lambda x: fastq[x - 1 : x]
            )
            called_events = len(event_data)

            # create event length data for tail normalization
            event_length_vector = np.empty(called_events)
            event_length_vector[:] = np.nan
            count = 0
            for i in range(called_events - 1, -1, -1):
                if event_data["move"][i] == 1:
                    event_length_vector[i] = count + 1
                    count = 0
                else:
                    count += 1
            # multiply moves by length of the event
            event_length_vector = event_length_vector * stride
            event_data["event_length_vector"] = event_length_vector
            # del event_data['move_cumsum']
            # remove NAs
            event_length_vector = event_length_vector[~np.isnan(event_length_vector)]
            # Normalizer for flip-flop based data
            samples_per_nt = np.mean(
                event_length_vector[
                    event_length_vector <= np.quantile(event_length_vector, 0.95)
                ]
            )
        else:
            raise ValueError("model type is not flipflop")

    return raw_data, event_data, fastq, start, stride, samples_per_nt


def _singleApplyFunc(subDtframe, func):
    subResults = subDtframe.apply(func, axis=1)
    return subResults


def multfunc_dtframe(func, data, threading, use_iter=False, use_threading=False, *args):
    if use_threading:
        from multiprocessing.dummy import Pool
    else:
        from multiprocessing import Pool
    result = []
    if not use_iter:
        pool = Pool(threading)
        span = (len(data) // threading) + 1
        for _ in range(threading):
            sub_data = data.iloc[_ * span : (_ + 1) * span]
            result.append(pool.apply_async(func, args=(sub_data, *args)))
        pool.close()
        pool.join()
        result = [x.get() for x in result]
    else:
        forward_chunk_result = []
        latter_chunk_result = []
        chunk_data = next(data)
        while True:

            pool = Pool(threading)
            if chunk_data.empty:
                break
            else:
                span = (len(chunk_data) // (threading - 1)) + 1
                for _ in range(threading - 1):
                    sub_data = chunk_data.iloc[_ * span : (_ + 1) * span]
                    if not sub_data.empty:
                        latter_chunk_result.append(
                            pool.apply_async(func, args=(sub_data, *args))
                        )
                try:
                    chunk_data = next(data)
                    if forward_chunk_result:
                        forward_chunk_result = [x.get() for x in forward_chunk_result]
                        result.extend(forward_chunk_result)
                    pool.close()
                    pool.join()
                    forward_chunk_result = latter_chunk_result
                    latter_chunk_result = []
                except:
                    if forward_chunk_result:
                        forward_chunk_result = [x.get() for x in forward_chunk_result]
                        result.extend(forward_chunk_result)
                    pool.close()
                    pool.join()
                    forward_chunk_result = latter_chunk_result
                    latter_chunk_result = []
                    break
        forward_chunk_result = [x.get() for x in forward_chunk_result]
        result.extend(forward_chunk_result)

    return result


def multiApplyFunc(allDtframe, func, threads):
    allResults = multfunc_dtframe(
        _singleApplyFunc, allDtframe, threads, False, False, func
    )
    allResults = pd.concat(allResults)
    return allResults


def transformExpressionMatrixTo10XMtx(inputPath, outputDir):
    """
    input:
        path or dataframe


    column: gene name
    index: barcode名 ( without -1 )
    """

    try:
        sh.mkdir(outputDir)
    except:
        sh.rm("-rf", outputDir)
        sh.mkdir(outputDir)

    if isinstance(inputPath, str):
        expressionMtx = pd.read_table(
            inputPath,
            index_col=0,
        )
    else:
        expressionMtx = inputPath
        expressionMtx.rename_axis("index", inplace=True)
    expressionMtx = expressionMtx.loc[:, expressionMtx.sum(0) != 0]
    barcodes = pd.Series(expressionMtx.index + "-1")
    barcodes.to_csv(f"{outputDir}barcodes.tsv", header=None, index=None)

    feature = pd.DataFrame(expressionMtx.columns)
    feature[1] = feature.iloc[:, 0]
    feature[2] = "Gene Expression"
    feature.to_csv(f"{outputDir}features.tsv", sep="\t", header=None, index=None)

    indexMap = {
        i: k
        for i, k in zip(expressionMtx.index, range(1, 1 + len(expressionMtx.index)))
    }

    featureMap = {
        i: k
        for i, k in zip(expressionMtx.columns, range(1, 1 + len(expressionMtx.columns)))
    }

    expressionMtx.index = expressionMtx.index.map(indexMap)
    expressionMtx.columns = expressionMtx.columns.map(featureMap)
    expressionMtx = expressionMtx.astype(int)
    expressionMtx.reset_index(inplace=True)
    expressionMtx = expressionMtx.melt(id_vars="index")

    expressionMtx.columns = ["barcode", "feature", "count"]
    expressionMtx = expressionMtx.query("count != 0")
    expressionMtx = expressionMtx.reindex(["feature", "barcode", "count"], axis=1)
    expressionMtx.sort_values(
        ["barcode", "feature"], ascending=[True, False], inplace=True
    )
    featureCounts, barcodeCounts, rowCounts = (
        max(expressionMtx["feature"]),
        max(expressionMtx["barcode"]),
        len(expressionMtx),
    )
    with open(f"{outputDir}matrix.mtx", "w") as fh:
        fh.write(
            f'%%MatrixMarket matrix coordinate integer general\n%metadata_json: {{"format_version": 2, "software_version": "X.X.0"}}\n{featureCounts} {barcodeCounts} {rowCounts}'
        )
        for line in expressionMtx.itertuples():
            fh.write(f"\n{line.feature} {line.barcode} {line.count}")

    sh.gzip(glob.glob(f"{outputDir}*"))


def updateOldMultiAd(adata):
    """
    update MultiAd from old version (all data deposit in X) to the 1.0 version (data deposit in obsm)
    """
    adata = adata.copy()

    def __addMatToObsm(adata, keyword):
        """
        read var name of adata, and add data matched the keyword to uns of adata
        """
        if keyword == "Abundance":
            subIndex = ~adata.var.index.str.contains("APA|Spliced")
        else:
            subIndex = adata.var.index.str.contains(keyword)
        subAd = adata[:, subIndex]
        adata.obsm[keyword] = subAd.X
        adata.uns[f"{keyword}_label"] = subAd.var.index

    __addMatToObsm(adata, "APA")
    __addMatToObsm(adata, "Spliced")
    __addMatToObsm(adata, "Abundance")
    adata = adata[:, ~adata.var.index.str.contains("APA|Spliced")]
    return adata


def getMatFromObsm(
    adata, keyword, useGeneLs=[], normalize=True, logScale=True, ignoreN=False
):
    """
    use MAT deposited in obsm replace the X MAT

    params:
        adata:
            version 1.0 multiAd
        keyword:
            stored in obsm
        useGeneLs:
            if not specified useGeneLs, all features will be output, otherwise only features associated with those gene will be output
        normalize:
            normalize the obtained Mtx or not
        logScale:
            log-transformed or not
        ignoreN:
            ignore ambiguous APA/Splice info
    return:
        anndata
    """
    transformedAd = anndata.AnnData(
        X=adata.obsm[keyword].copy(),
        obs=adata.obs,
        var=adata.uns[f"{keyword}_label"].to_frame().drop(0, axis=1),
        obsp=adata.obsp,
        obsm=adata.obsm,
        uns=adata.uns,
    )
    if normalize:
        sc.pp.normalize_total(transformedAd, target_sum=1e4)
    if logScale:
        sc.pp.log1p(transformedAd)

    useGeneLs = list(useGeneLs)
    if not useGeneLs:
        transformedAd = transformedAd
    else:
        transformedAdFeatureSr = transformedAd.var.index
        transformedAdFeatureFilterBl = (
            transformedAdFeatureSr.str.split("_").str[0].isin(useGeneLs)
        )
        transformedAd = transformedAd[:, transformedAdFeatureFilterBl]

    if ignoreN:
        transformedAdFeatureSr = transformedAd.var.index
        transformedAdFeatureFilterBl = (
            ~transformedAdFeatureSr.str.split("_").str[1].isin(["N", "Ambiguous"])
        )

        transformedAd = transformedAd[:, transformedAdFeatureFilterBl]

    return transformedAd