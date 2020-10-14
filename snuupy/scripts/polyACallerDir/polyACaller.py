from Bio import Align
from Bio.Seq import Seq
import numpy as np
import pickle
import sys
import joblib
import click
import pandas as pd
import sys 
sys.path.append("..") 
from scripts.tools import extract_read_data
from ont_fast5_api.fast5_interface import get_fast5_file
from scripts.tools import multiApplyFunc

def polyACaller(infile, outfile, threads, fp, ep):
    """Estimate poly(A) tail lengths from nanopore cDNA reads"""
    global fp_
    global ep_

    fp_ = fp
    ep_ = ep
    allIndex = pd.read_hdf(infile)
    allResults = multiApplyFunc(allIndex, estimate_tail_length, threads)
    allResults.columns = ['readId','readType','tailStart','tailEnd','tailBaseStart','tailBaseEnd','primerStart','primerEnd','samplePerNt','tailLength','fast5FilePath']
    allResults.to_hdf(outfile,key='defaults')

def find_dna_tailtype(
    fast5_filepath,
    read_id,
    threshold = .6):

    global fp_
    global ep_

    fp = fp_
    ep = ep_

    raw_data, event_data, fastq, start, stride, samples_per_nt = extract_read_data(fast5_filepath, read_id)
    fastq = Seq(fastq)
    # set alignment parameter
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -1.0
    aligner.mode = 'local'
    # get primer
    # fp = Seq(fp)
    ep = ep[-12:] #
    ep = Seq(ep)
    rc_ep = ep.reverse_complement()
    # alignments were performed using Smith–Waterman local alignments
    # as_fp = aligner.align(fp, fastq[:100])
    query_len = 100
    as_ep = aligner.align(ep, fastq[ : query_len])
    as_rc_ep = aligner.align(rc_ep, fastq[-query_len: ])
    # normalized align score
    #nas_fp = as_fp.score/len(fp)
    nas_ep = as_ep.score/len(ep)
    nas_rc_ep = as_rc_ep.score/len(rc_ep)

    # has_precise_boundary = False
    # bases_to_match = 3
    if nas_rc_ep > nas_ep and nas_rc_ep > threshold:
        # check if there is the end primer at the end of the polyA read
        # adjacent to the polyA tail
        # If it is a valid polyA tail, then find the rough starting site of the tail
        read_type = 'polyA'
        tail_is_valid = True

        query_start = as_rc_ep[0].aligned[1][0][0]  # query为read
        target_start = as_rc_ep[0].aligned[0][0][0]  # target为primer

        # 因为引物序列是准确的，所以如果比对结果引物序列有soft-clip，得加上soft-clip序列长度
        polya_end_fastq =  query_start + (len(fastq) - query_len) - target_start  # 引物5‘end 位置
        polyt_start_fastq = None

    elif nas_rc_ep < nas_ep and nas_ep > threshold:
        read_type = 'polyT'
        tail_is_valid = True

        query_end = as_ep[0].aligned[1][-1][-1]  # query为read
        target_end = as_ep[0].aligned[0][-1][-1]  # target为primer

        polyt_start_fastq = query_end + (len(ep) - target_end) + 1  # 引物3’end 位置+1就是polyt起始位置
        polya_end_fastq = None
    else:
        read_type = 'invalid'
        tail_is_valid = False
        polyt_start_fastq = None
        polya_end_fastq = None

    return read_type, tail_is_valid, polya_end_fastq, polyt_start_fastq, event_data, samples_per_nt

def find_tail(tail_seq, tail_length_vector, base_type):
    '''
    Dynamic Programming:
        d[i,j]: 从第i到第j个碱基的tail score, i<=j
        d[i,j]=score, if i==j
        d[i,j]=d[i,j-1]+score
    '''
    MATCH = 1
    MISMATCH = -1.5

    tail_length = len(tail_seq)
    results = np.zeros([tail_length, tail_length])
    score = lambda x : MATCH if tail_seq[x] == base_type else MISMATCH
    maxscore = float('-inf')
    tail_boundary = None
    for i in range(tail_length):
        for j in range(i, tail_length):
            if i == j:
                results[i, j] = score(j) * tail_length_vector[j]
            else:
                results[i, j] = results[i, j-1] + score(j) * tail_length_vector[j]
            if results[i, j] > maxscore:
                maxscore = results[i, j]
                tail_boundary = (i, j)

    return tail_boundary

def NumericOutlier(value):
    '''
    1、计算第一四分位数（Q1）及第三四分位数（Q3）
    2、计算IQR （IQR = Q3 - Q1）
    3、输出正常区间[Q1-1.5IQR，Q3+1.5IQR]
    '''
    iqr = np.quantile(value,0.75) - np.quantile(value,0.25)
    quan_down = np.quantile(value,0.25)-1.5*iqr
    quan_up = np.quantile(value,0.75)+1.5*iqr
    return float(quan_down),float(quan_up)

def estimate_tail_length(idx_data):
    read_id = idx_data.name
    fast5_filepath = idx_data['fast5_filepath']
    read_type, tail_is_valid, polya_end_fastq, polyt_start_fastq, event_data, samples_per_nt = find_dna_tailtype(fast5_filepath, read_id)
    # 由于polyA/T之间还隔着28nt，所以加上18nt。
    shift = 18
    if tail_is_valid:
        if read_type == 'polyA':
            base_type = 'A'
            polya_start_fastq = idx_data['qend']+1
            primer_start = polya_end_fastq
            primer_end = np.nan
            shift_start = polya_start_fastq-shift if polya_start_fastq-shift > 0 else 0
            shift_end = polya_end_fastq-shift if polya_end_fastq+shift < len(event_data) else len(event_data)
        elif read_type == 'polyT':
            base_type = 'T'
            polyt_end_fastq = idx_data['qstart']+1
            primer_start = np.nan
            primer_end = polyt_start_fastq
            shift_start = polyt_start_fastq+shift if polyt_start_fastq-shift > 0 else 0
            shift_end = polyt_end_fastq+shift if polyt_end_fastq+shift < len(event_data) else len(event_data)

        tail_event = event_data[(event_data['move_cumsum']>=shift_start)&(event_data['move_cumsum']<=shift_end)].dropna()

        # 计算event_length_vector正常值范围
        event_length_vector_outlier = NumericOutlier(event_data['event_length_vector'].dropna())
        # 异常值event_length_vector的碱基默认为base_type
        tail_event['bases'] = tail_event.index.map(
            lambda x: base_type
            if tail_event['event_length_vector'].loc[x] < event_length_vector_outlier[0]
            or tail_event['event_length_vector'].loc[x] > event_length_vector_outlier[1]
            else tail_event['model_state'].loc[x])
        tail_seq = ''.join(tail_event['bases'])
        tail_length_vector = tail_event['event_length_vector'].tolist()
        # 计算tail区间
        tail_boundary = find_tail(tail_seq, tail_length_vector, base_type)

        # 有polyA尾巴的才输出
        if tail_boundary:
            tail_event_length = sum(tail_event['event_length_vector'].iloc[tail_boundary[0]:tail_boundary[1]+1])
            tail_length = tail_event_length/samples_per_nt
            tail_start = tail_event['start'].iloc[tail_boundary[0]]
            tail_end = tail_event['start'].iloc[tail_boundary[1]]
            tail_base_start = tail_event['move_cumsum'].iloc[tail_boundary[0]]
            tail_base_end = tail_event['move_cumsum'].iloc[tail_boundary[1]]
        else :
            read_type = 'non-polyA/T'
            tail_start,tail_end,tail_base_start,tail_base_end,tail_length = [0]*5
    else:
        tail_start,tail_end,tail_base_start,tail_base_end,primer_start,primer_end,tail_length = [0] * 7


    return pd.Series([read_id,read_type,tail_start,tail_end,tail_base_start,tail_base_end,primer_start,primer_end,samples_per_nt,tail_length,fast5_filepath])

