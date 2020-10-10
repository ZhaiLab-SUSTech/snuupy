'''
@Date: 2020-07-02 15:47:49
LastEditors: Liuzj
LastEditTime: 2020-09-18 20:58:33
@Description: file content
@Author: liuzj
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/00_pipeline/finalVersion/step12_generateExpressionMatrix.py
'''
import pickle
import click
import pandas as pd
from jpy_tools.ReadProcess import transformExpressionMatrixTo10XMtx

@click.command()
@click.option('-i', 'PARSED_RESULT')
@click.option('-o', 'OUT_PATH')
def main(PARSED_RESULT, OUT_PATH):
    with open(PARSED_RESULT, 'rb') as fh:
        parsedResult = pickle.load(fh)
    parsedResultDt = pd.DataFrame.from_dict(parsedResult, 'index')

    parsedResultDt['barcode'] = parsedResultDt.index.str.split('_').str[0]
    parsedGeneResultDt = parsedResultDt.groupby(['barcode', 'gene_id']).count()
    parsedIsoformResultDt = parsedResultDt.groupby(['barcode', 'isoform_id']).count()

    parsedGeneResultDt.rename_axis(['barcode', 'varId'], inplace=True)
    parsedIsoformResultDt.rename_axis(['barcode', 'varId'], inplace=True)

    parsedGeneResultDt = parsedGeneResultDt['cov']
    parsedIsoformResultDt = parsedIsoformResultDt['cov']

    parsedResult = pd.concat([parsedGeneResultDt, parsedIsoformResultDt])



    writeExpressionInfo = lambda x: (x.pipe(lambda x: x.unstack()).pipe(lambda x: x.rename_axis('')).pipe(
        lambda x: x.rename_axis('',axis=1)).pipe(lambda x :x.fillna(0)))

    transformExpressionMatrixTo10XMtx(writeExpressionInfo(parsedResult), f'{OUT_PATH}allExpressionMtx/')
    transformExpressionMatrixTo10XMtx(writeExpressionInfo(parsedGeneResultDt), f'{OUT_PATH}geneExpressionMtx/')
    transformExpressionMatrixTo10XMtx(writeExpressionInfo(parsedIsoformResultDt), f'{OUT_PATH}isoformExpressionMtx/')
main()