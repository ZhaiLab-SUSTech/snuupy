
import pickle
import pandas as pd
import numpy as np
from functools import reduce
import click

@click.command()
@click.option('-i', 'INTRON_INFO')
@click.option('-g', 'GENE_NAME_INFO')
@click.option('-o', 'OUT_PATH')
def main(INTRON_INFO, GENE_NAME_INFO, OUT_PATH):
    with open(GENE_NAME_INFO, 'rb') as fh:
        geneNameInfo = pickle.load(fh)
    intronInfo = pd.read_table(INTRON_INFO)
    intronInfo["geneIdNonTrans"] = intronInfo["GeneId"].str.split(".").str[0]
    intronInfo["baselineGeneId"] = intronInfo["Name"].map(
        lambda x: geneNameInfo.get(x, {"gene_id": 0})["gene_id"]
    )
    intronInfo.query("geneIdNonTrans == baselineGeneId", inplace=True)
    intronInfo.drop(["baselineGeneId", "GeneId"], axis=1, inplace=True)
    intronInfo.rename({"geneIdNonTrans": "geneId"}, axis=1, inplace=True)
    intronInfo.to_csv(OUT_PATH, sep="\t", index=False)

main()