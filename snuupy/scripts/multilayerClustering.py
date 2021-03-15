# This method is a pre-release version
# we have not verified this method
# This method is inspired by scLAPA (https://github.com/BMILAB/scLAPA)
# In short, we first separately calculate correlation mat  for dimension reduced expression mat, APA mat and spliced mat.
# then use SNF method fusion these matrices, the resulting mat can be load to clustering algorithms, such as leiden.

import numpy as np
import pandas as pd
import scanpy as sc
from .tools import updateOldMultiAd, getMatFromObsm, normalizeByScran, transformEntToAd


def useSNF(apaAd, spliceAd, abunAd, adata, outPath):
    import snf
    [sc.pp.scale(x, max_value=10) for x in [apaAd, spliceAd, abunAd]]
    [sc.tl.pca(x, svd_solver="arpack", n_comps=50) for x in [apaAd, spliceAd, abunAd]]
    [sc.pp.neighbors(x, n_pcs=30) for x in [apaAd, spliceAd, abunAd]]
    similarityMat = [x.obsp["connectivities"].A for x in [abunAd, apaAd, spliceAd]]

    fusedMat = snf.snf(similarityMat, K=20, alpha=0.5, t=10)
    np.save(f"{outPath}_fusedMat.npy", fusedMat)
    sc.tl.leiden(adata, adjacency=fusedMat, key_added="leiden_fused")
    adata.write(f"{outPath}_SNF_leiden_resolution_1.h5ad")


def useMOFA(apaAd, spliceAd, abunAd, adata, outPath):
    from mofapy2.run.entry_point import entry_point
    normalizedAd = sc.concat([apaAd, spliceAd, abunAd], axis=1)

    ent = entry_point()
    mofaUseDf = (
        normalizedAd.to_df()
        .reset_index()
        .melt("index")
        .rename(dict(index="sample", variable="feature"), axis=1)
        .assign(
            group="group_0",
            view=lambda df: np.select(
                [df["feature"].str.contains("_")],
                [df["feature"].str.split("_").str[-1]],
                "abundance",
            ),
        )
    ).reindex(["sample", "group", "feature", "view", "value"], axis=1)
    ent.set_data_options(
        scale_groups=False,
        scale_views=True
    )
    ent.set_data_df(mofaUseDf, likelihoods=["gaussian", "gaussian", "gaussian"])
    ent.set_model_options(
        factors=30,
        spikeslab_weights=True,
        ard_factors=True,
        ard_weights=True
    )
    ent.set_train_options(
        iter=1000,
        convergence_mode="fast",
        startELBO=1,
        freqELBO=1,
        dropR2=0.001,
        gpu_mode=True,
        verbose=False,
        seed=1
    )
    ent.build()
    ent.run()
    mofaAd = transformEntToAd(ent)
    mofaAd.write(f"{outPath}_MOFA_result.h5ad")

    abunAd.obsm['X_mofa'] = mofaAd.X.T
    sc.pp.neighbors(abunAd, 28, use_rep="X_mofa", key_added='mofa')
    sc.tl.umap(abunAd, neighbors_key='mofa')
    abunAd.obsm['X_umap_mofa'] = abunAd.obsm['X_umap']
    sc.tl.leiden(abunAd, neighbors_key='mofa', key_added="leiden_mofa")
    adata.write(f"{outPath}_MOFA_leiden_resolution_1.h5ad")


def main(multiMatPath, outPath, method):
    """
    SNF:
        In short, we first separately calculate euclidean mat for dimension reduced gene expression mat, APA mat and spliced mat.
        then use SNF method fusion these matrices, the resulting mat can be load to clustering algorithms, such as leiden.
        This method is inspired by scLAPA (https://github.com/BMILAB/scLAPA)

        we have not verified this method and it is a pre-release version.
    
    MOFA:
        In short, we use MOFA+ (https://biofam.github.io/MOFA2) method to APA, Splice and abundance matrix .

    method:
        snf|mofa
    multiMatPath:
        multilayer mat path
    outPath:
        prefix of output file containing fused connectivities matrix(snf only), mofa results (mofa only) and leiden clustering result.
    """

    adata = updateOldMultiAd(sc.read_10x_mtx(multiMatPath))


    apaAd, spliceAd, abunAd = [
        getMatFromObsm(
            adata, x, strCommand="Nc"
        )
        for x in ["APA", "Spliced", "Abundance"]
    ]


    apaAd, spliceAd, abunAd = [
        normalizeByScran(x) for x in [apaAd, spliceAd, abunAd]
    ]

    [
        sc.pp.highly_variable_genes(x, flavor="cell_ranger", n_top_genes=3000)
        for x in [apaAd, spliceAd, abunAd]
    ]

    apaAd, spliceAd, abunAd = [
        x[:, x.var["highly_variable"]] for x in [apaAd, spliceAd, abunAd]
    ]

    useMethod = {
        "snf": useSNF,
        "mofa": useMOFA
    }[method]

    useMethod(apaAd, spliceAd, abunAd, adata, outPath)
