import pandas as pd



def barcodeAssignment(MISMATCH_RESULT, OUTPUT_FEATHER, MAX_BARCODE_ED, MAX_UMI_ED):
    mismatchResult = pd.read_feather(MISMATCH_RESULT)

    mismatchResult.drop_duplicates(['name','qseqid'], inplace=True)

    mismatchResult['barcodeUmiMismatch'] = mismatchResult['barcodeUmiMismatch'].astype(int)
    mismatchResult['barcodeMismatch'] = mismatchResult['barcodeMismatch'].astype(int)
    mismatchResult['umiMismatch'] = mismatchResult['umiMismatch'].astype(int)

    mismatchResult['hitContent'] = mismatchResult['qseqid'] + ',' + mismatchResult['barcodeUmiMismatch'].astype(str) + ','\
        + mismatchResult['barcodeMismatch'].astype(str) + ',' + mismatchResult['umiMismatch'].astype(str) + ';'
    
    mismatchResult['allHitContents'] = mismatchResult.groupby('name')['hitContent'].transform('sum')
    mismatchResult['hitCounts'] = mismatchResult.groupby('name')['hitContent'].transform('count')

    mismatchResult = mismatchResult.loc[(mismatchResult['barcodeMismatch'] <= MAX_BARCODE_ED) & \
            (mismatchResult['umiMismatch'] <= MAX_UMI_ED)]
    mismatchResult.sort_values(['name','barcodeUmiMismatch','barcodeMismatch','umiMismatch'], inplace=True)
    mismatchResult.drop_duplicates('name', inplace=True)
    
    mismatchResult.reset_index(drop=True, inplace=True)
    mismatchResult.to_feather(OUTPUT_FEATHER)