# This method is a pre-release version
# we have not verified this method
# This method is inspired by scLAPA (https://github.com/BMILAB/scLAPA)
# In short, we first separately calculate correlation mat  for dimension reduced expression mat, APA mat and spliced mat.
# then use SNF method fusion these matrices, the resulting mat can be load to clustering algorithms, such as leiden.

import numpy as np
import pandas as pd
import scanpy as sc
import snf
from .tools import updateOldMultiAd, getMatFromObsm


def main(multiMatPath, useGenePath, outPath):
    """
    In short, we first separately calculate euclidean mat for dimension reduced gene expression mat, APA mat and spliced mat.
    then use SNF method fusion these matrices, the resulting mat can be load to clustering algorithms, such as leiden.
    This method is inspired by scLAPA (https://github.com/BMILAB/scLAPA)

    we have not verified this method and it is a pre-release version.

    multiMatPath:
        multilayer mat path
    useGenePath:
        gene used for calculating correlation matrix. if not provided, all gene will be used. NO HEADER, NO INDEX
        e.g:
            AT1G01010
            AT1G01020
            AT1G01030
    outPath:
        prefix of output file containing fused connectivities matrix and leiden clustering result.
        matrix: fused eulidean mat npy format, could be loaded by numpy.load function
    """

    def __calSimMat(adata):
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
        sc.pp.neighbors(adata, n_pcs=30)

    adata = updateOldMultiAd(sc.read_10x_mtx(multiMatPath))

    if useGenePath:
        useGeneLs = pd.read_table(useGenePath, header=None, names=["gene"])["gene"]
        adata = adata[:, useGeneLs]

    spliceAd = getMatFromObsm(
        adata, "Spliced", adata.var.index, ignoreN=True, clear=True
    )

    apaAd = getMatFromObsm(adata, "APA", adata.var.index, ignoreN=True, clear=True)

    abunAd = getMatFromObsm(
        adata, "Abundance", adata.var.index, ignoreN=True, clear=True
    )

    [__calSimMat(x) for x in [spliceAd, apaAd, abunAd]]

    similarityMat = [x.obsp["connectivities"].A for x in [abunAd, apaAd, spliceAd]]

    fusedMat = snf.snf(similarityMat, K=20, alpha=0.5, t=10)
    np.save(f"{outPath}_fusedMat.npy", fusedMat)

    sc.tl.leiden(adata, adjacency=fusedMat, key_added="leiden_fused")

    clusterDf = adata.obs[["leiden_fused"]]
    clusterDf.to_csv(f"{outPath}_leiden_resolution_1.0.tsv", sep="\t")
