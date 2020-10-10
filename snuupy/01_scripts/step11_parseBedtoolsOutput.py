
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Description: 添加isoform注释解析 <naive>
Date: 2020-09-18 19:34:16
LastEditTime: 2020-09-25 13:29:43
LastEditors: Liuzj
'''
'''
@Author       : windz
@Date         : 2020-06-04 11:53:34
LastEditTime: 2020-09-18 19:34:23
@Description  : 
'''

from collections import defaultdict
import pickle
import numpy as np
import pandas as pd
import joblib
import click


import logging
logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s: %(message)s',  
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    )


NAMES=[
    'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 
    'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts', 
    'geneChromosome', 'geneStart', 'geneEnd', 'geneName', 'geneScore', 'geneStrand', 
    'geneThickStart', 'geneThickEnd', 'geneItemRGB', 'geneBlockCount', 'geneBlockSizes',   'geneBlockStarts', 'cov'
    ]
USECOLS = [
    'Chromosome', 'Start', 'End', 'Name', 'Strand',
    'geneStart', 'geneEnd', 'geneName', 'geneStrand', 'geneBlockSizes', 'geneBlockStarts', 'cov'
    ]



@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
def main(infile, outfile):
    logging.info('Start read csv')
    df = pd.read_csv(
        infile, 
        sep='\t', 
        names=NAMES,
        usecols=USECOLS,
        header=None
        )
    logging.info('Read csv Done!')

    logging.info('Start find Splice Sites')
    df['geneBlockSizes'] = df['geneBlockSizes'].map(lambda x: np.fromstring(x, sep=','))
    df['geneBlockStarts'] = df['geneBlockStarts'].map(lambda x: np.fromstring(x, sep=','))
    df['five_ss'] = (df['geneStart']+df['geneBlockSizes']+df['geneBlockStarts']).map(lambda x: x[:-1])
    df['three_ss'] = (df['geneStart']+df['geneBlockStarts']).map(lambda x: x[1:])
    logging.info('Find Splice Sites Done!')

    logging.info('Main function')

    df['geneName'] = df['geneName'].str.split('\|\|').str[0]
    specificLine = df['geneName'].str.endswith('_specific')
    specificDf = df.loc[specificLine]
    isoformDf = df.loc[~specificLine]

    results = defaultdict(
    lambda : {
        'cov': 0,
        'is_splicing_intermediates': False,
        'gene_len': 0,
        'gene_id': None,
        'isoform_mapping_ratio':{},
        }
    )
    for item in isoformDf.itertuples():
        # 判断是否为剪切中间体
        if item.Strand == '+':
            is_splicing_intermediates = (abs(item.End-item.five_ss)<=10).any()
        else:
            is_splicing_intermediates = (abs(item.Start-item.three_ss)<=10).any()
        results[item.Name]['is_splicing_intermediates'] = results[item.Name]['is_splicing_intermediates'] or is_splicing_intermediates
        results[item.Name]['isoform_mapping_ratio'][item.geneName] = item.cov / int(item.geneBlockSizes.sum())


        # 取exon覆盖度最大的作为这条read的gene注释
        if results[item.Name]['cov'] < item.cov:
            results[item.Name]['cov'] = item.cov
            results[item.Name]['gene_id'] = item.geneName.split('.')[0]
            results[item.Name]['gene_len'] = item.geneEnd - item.geneStart

    logging.info('Gene Assign Done!')

    specificDf['geneNameTrue'] = specificDf['geneName'].str.split('.').str[0]
    specificDf['geneAssignGene'] = specificDf.pipe(lambda x: x['Name'].map(lambda y:results[y]['gene_id']))
    specificDf = specificDf.query("geneAssignGene == geneNameTrue")
    specificDf['isoformName'] = specificDf['geneName'].str.split('_').str[0]
    specificDf = specificDf.loc[:, ['Name', 'cov', 'isoformName']]
    specificDfGroup = specificDf.groupby('Name').apply(lambda z: {x:y for x,y in zip(z['isoformName'], z['cov'])})
    specificDfGroup = specificDfGroup.to_dict()

    for readName, geneBedInfo in results.items():
        isoformMappingRatio = geneBedInfo['isoform_mapping_ratio']
        if len(isoformMappingRatio) == 1:
            isoformName = list(isoformMappingRatio.keys())[0]
        else:
            readSpecificDfGroup = specificDfGroup.get(readName, {})
            readSpecificDfGroup = {x:readSpecificDfGroup.get(x, 0) for x in isoformMappingRatio.keys()}
            maxIsoformCoverageLength = max(readSpecificDfGroup.values())
            cutLength = maxIsoformCoverageLength - 15
            putativeIsoforms = [x for x,y in readSpecificDfGroup.items() if y > cutLength]
            if len(putativeIsoforms) == 1:
                isoformName = putativeIsoforms[0]
            else:
                putativeIsoforms.sort(key = lambda x: isoformMappingRatio[x])
                putativeIsoformMappingRatio = np.array([isoformMappingRatio[x] for x in putativeIsoforms])
                if putativeIsoformMappingRatio[-1] - putativeIsoformMappingRatio[-2] > 0.1:
                    isoformName = putativeIsoforms[-1]
                else:
                    isoformName = geneBedInfo['gene_id'] + '.N'
        results[readName]['isoform_id'] = isoformName



    logging.info('Main function Done!')
    
    with open(outfile, 'wb') as o:
        pickle.dump(dict(results), o)
    

if __name__ == "__main__":
    main()