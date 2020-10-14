'''
Description: 
Author: Liuzj
Date: 2020-10-14 18:18:06
LastEditTime: 2020-10-14 18:18:17
LastEditors: Liuzj
'''
import subprocess
from io import StringIO
import pickle
import click
import pandas as pd
import os


def read_mm2_output(mm2out):
    '''
    read the tsv file format produced by minimap2
    '''
    MM2_COLUMNS = ['qname', 'qstart', 'qend']
    df = pd.read_csv(mm2out,
                     sep='\t',
                     header=None,
                     names=MM2_COLUMNS,
                     usecols=[0, 2, 3])
    # remove multi-alignment

    df = df.drop_duplicates(subset=['qname'], keep='first')
    print(len(df))
    return df


def build_index(infile, genome, threads, f5dir, f5summary, bed, outfile, minimapPath):
    '''Build index mapping for basecalled reads'''
    # get config
    ref = genome
    mm2thread = str(threads)
    fast5_dir = f5dir
    sequencing_summary = f5summary
    bed = bed
    f5dir = os.popen(f'realpath {f5dir}').read()
    f5dir = f5dir.rstrip() + '/'
    proc = subprocess.run([
        minimapPath, '-x', 'splice', '-t', mm2thread, '--junc-bed', bed, '-uf',
        '-k14', '--secondary=no', '-G', '10000', ref, infile
    ],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    idx_data = read_mm2_output(StringIO(proc.stdout.decode()))
    summary_dt = pd.read_table(f5summary, low_memory=False)
    summary_dt = summary_dt.loc[:, ['filename', 'read_id']]
    summary_dt.set_index('read_id', inplace=True)
    idx_data.set_index('qname', inplace=True)
    idx_data = pd.merge(idx_data,
                        summary_dt,
                        how='left',
                        left_index=True,
                        right_index=True)
    idx_data.filename = f5dir + idx_data.filename
    idx_data.rename({'filename': 'fast5_filepath'}, axis=1, inplace=True)
    idx_data.to_hdf(outfile, key='default')