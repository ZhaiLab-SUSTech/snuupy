'''
@Date: 2020-06-10 14:06:50
LastEditors: Liuzj
LastEditTime: 2020-09-11 19:07:41
@Description: 目前使用的策略 总体mismatch不大于6 分别不大于4
@Author: liuzj
@FilePath: /liuzj/projects/split_barcode/01_20200507/01_pipeline/00_pipeline/finalVersion/step07_parseMismatchResult.py
'''
import click
import pandas as pd


@click.command()
@click.option('-i', 'MISMATCH_RESULT', help='step06 output')
@click.option('-o', 'OUTPUT_FEATHER', help='nanopore read id with barcode and umi; feather format')
def main(MISMATCH_RESULT, OUTPUT_FEATHER):
    mismatchResult = pd.read_feather(MISMATCH_RESULT)

    mismatchResult.drop_duplicates(['name','qseqid'], inplace=True)

    mismatchResult['barcodeUmiMismatch'] = mismatchResult['barcodeUmiMismatch'].astype(int)
    mismatchResult['barcodeMismatch'] = mismatchResult['barcodeMismatch'].astype(int)
    mismatchResult['umiMismatch'] = mismatchResult['umiMismatch'].astype(int)

    mismatchResult['hitContent'] = mismatchResult['qseqid'] + ',' + mismatchResult['barcodeUmiMismatch'].astype(str) + ','\
        + mismatchResult['barcodeMismatch'].astype(str) + ',' + mismatchResult['umiMismatch'].astype(str) + ';'
    
    mismatchResult['allHitContents'] = mismatchResult.groupby('name')['hitContent'].transform('sum')
    mismatchResult['hitCounts'] = mismatchResult.groupby('name')['hitContent'].transform('count')

    mismatchResult = mismatchResult.loc[(mismatchResult['barcodeUmiMismatch'] <= 6) & \
        (mismatchResult['barcodeMismatch'] <= 3) & \
            (mismatchResult['umiMismatch'] <= 3)]
    mismatchResult.sort_values(['name','barcodeUmiMismatch','barcodeMismatch','umiMismatch'], inplace=True)
    mismatchResult.drop_duplicates('name', inplace=True)
    
    mismatchResult.reset_index(drop=True, inplace=True)
    mismatchResult.to_feather(OUTPUT_FEATHER)


main()