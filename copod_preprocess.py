#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 16:02:38 2022

@author: yyj
"""
import sys
import numpy as np
import pysam
import os
import gc
import pandas as pd
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
import datetime
from sklearn.cluster import KMeans

from myModel import COPOD


#read the bam file，return the list of chromosomes
def read_bam_file(filename):
    samfile = pysam.AlignmentFile(filename, "rb")
    chrList = samfile.references
    return chrList

#read the fasta file，return a two-dim list including all information from the bam file
def read_ref_file(filename, chr_num, ref):
    # read reference file
    if os.path.exists(filename):
        print("Read reference file: " + str(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            for line in f:
                linestr = line.strip()
                ref[chr_num] += linestr
    else:
        print("Warning: can not open " + str(filename) + '\n')
    return ref

#
def Binning(ref, binSize, chrLen, filename):
    chrTag = np.full(23, 0)
    chrList = np.arange(23)
    maxNum = int(chrLen.max() / binSize) + 1
    InitRD = np.full((23, maxNum), 0.0)
    # read bam file and get bin rd
    print("Read bam file: " + str(filename))
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        idx = int(line.pos / binSize)
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                InitRD[int(chr)][idx] += 1
                chrTag[int(chr)] = 1

    chrList = chrList[chrTag > 0]
    chrNum = len(chrList)
    RDList = [[] for i in range(chrNum)]
    PosList = [[] for i in range(chrNum)]
    InitGC = np.full((chrNum, maxNum), 0)
    pos = np.full((chrNum, maxNum), 0)

    # initialize bin_data and bin_head
    count = 0
    for i in range(len(chrList)):
        chr = chrList[i]
        binNum = int(chrLen[chr] / binSize) + 1
        for j in range(binNum):
            pos[i][j] = j
            cur_ref = ref[chr][j * binSize:(j + 1) * binSize]
            N_count = cur_ref.count('N') + cur_ref.count('n')
            if N_count == 0:
                gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
            else:
                gc_count = 0
                InitRD[chr][j] = -1000000
                count = count + 1
            InitGC[i][j] = int(round(gc_count / binSize, 3) * 1000)

        # delete
        cur_RD = InitRD[chr][:binNum]
        cur_GC = InitGC[i][:binNum]
        cur_pos = pos[i][:binNum]
        cur_RD = cur_RD / 1000
        index = cur_RD >= 0
        RD = cur_RD[index]
        GC = cur_GC[index]
        cur_pos = cur_pos[index]
        RD[RD == 0] = modeRD(RD)
        RD = gc_correct(RD, GC)
        PosList[i].append(cur_pos)
        RDList[i].append(RD)
    del InitRD, InitGC, pos
    gc.collect()
    return RDList, PosList, chrList


def modeRD(RD):
    newRD = np.full(len(RD), 0)
    for i in range(len(RD)):
        newRD[i] = int(round(RD[i], 3) * 1000)

    count = np.bincount(newRD)
    countList = np.full(len(count) - 49, 0)
    for i in range(len(countList)):
        countList[i] = np.mean(count[i:i + 50])
    modemin = np.argmax(countList)
    modemax = modemin + 50
    mode = (modemax + modemin) / 2
    mode = mode / 1000
    return mode


def gc_correct(RD, GC):
    # correcting gc bias
    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean = np.mean(RD[GC == GC[i]])
        RD[i] = global_rd_ave * RD[i] / mean
    return RD


def scaling_RD(RD, mode):
    posiRD = RD[RD > mode]
    negeRD = RD[RD < mode]
    if len(posiRD) < 50:
        mean_max_RD = np.mean(posiRD)
    else:
        sort = np.argsort(posiRD)
        maxRD = posiRD[sort[-50:]]
        mean_max_RD = np.mean(maxRD)

    if len(negeRD) < 50:
        mean_min_RD = np.mean(negeRD)
    else:
        sort = np.argsort(negeRD)
        minRD = negeRD[sort[:50]]
        mean_min_RD = np.mean(minRD)
    scaling = mean_max_RD / (mode + mode - mean_min_RD)
    for i in range(len(RD)):
        if RD[i] < mode:
            RD[i] /= scaling
    return RD


def plot(pos, data):
    pos1 = np.arange(1, len(data)+1)
    plt.scatter(pos1, data, s=5, c="black")
    plt.xlabel("pos")
    plt.ylabel("scalRD")
    plt.show()


def seg_RD(RD, binHead, seg_start, seg_end, seg_count, binSize):
    seg_RD = np.full(len(seg_count), 0.0)
    for i in range(len(seg_RD)):
        seg_RD[i] = np.mean(RD[seg_start[i]:seg_end[i]])
        seg_start[i] = binHead[seg_start[i]] * binSize + 1
        if seg_end[i] == len(binHead):
            seg_end[i] = len(binHead) - 1
        seg_end[i] = binHead[seg_end[i]] * binSize + binSize

    return seg_RD, seg_start, seg_end


def Write_data_file(chr, seg_start, seg_end, seg_count, scores, p_value_file):
    """
    write data file
    pos_start, pos_end, lof_score, p_value
    """
    output = open(p_value_file, "w")
    output.write("start" + '\t' + "end" + '\t' + "read depth" + '\t' + "lof score" + '\t' + "p value" + '\n')
    for i in range(len(scores)):
        output.write(
            str(chr[i]) + '\t' + str(seg_start[i]) + '\t' + str(seg_end[i]) +
            '\t' + str(seg_count[i]) + '\t' + str(scores[i]) + '\n')


def Write_CNV_File(chr, CNVstart, CNVend, CNVtype, CN, filename):
    """
    write cnv result file
    pos start, pos end, type, copy number
    """
    output = open(filename, "w")
    for i in range(len(CNVtype)):
        if CNVtype[i] == 2:
            output.write("chr" + str(chr[i]) + '\t' + str(CNVstart[i]) + '\t' + str(
                CNVend[i]) + '\t' + str("gain") + '\t' + str(CN[i]) + '\n')
        else:
            output.write("chr" + str(chr[i]) + '\t' + str(CNVstart[i]) + '\t' + str(
                CNVend[i]) + '\t' + str("loss") + '\t' + str(CN[i]) + '\n')


def Read_seg_file(num_col, num_bin, seg_file):
    """
    read segment file (Generated by DNAcopy.segment)
    seg file: col, chr, start, end, num_mark, seg_mean
    """
    seg_start = []
    seg_end = []
    seg_count = []
    seg_len = []
    with open(seg_file, 'r') as f:
        for line in f:
            linestrlist = line.strip().split('\t')
            start = (int(linestrlist[0]) - 1) * num_col + int(linestrlist[2]) - 1
            end = (int(linestrlist[0]) - 1) * num_col + int(linestrlist[3]) - 1
            if start < num_bin:
                if end > num_bin:
                    end = num_bin - 1
                seg_start.append(start)
                seg_end.append(end)
                seg_count.append(float(linestrlist[5]))
                seg_len.append(int(linestrlist[4]))
    seg_start = np.array(seg_start)
    seg_end = np.array(seg_end)

    return seg_start, seg_end, seg_count, seg_len


def calculating_CN(mode, CNVRD, CNVtype):

    CN = np.full(len(CNVtype), 0)
    index = CNVtype == 1
    lossRD = CNVRD[index]
    if len(lossRD) > 2:
        data = np.c_[lossRD, lossRD]
        del_type = KMeans(n_clusters=2, random_state=9).fit_predict(data)
        CNVtype[index] = del_type
        if np.mean(lossRD[del_type == 0]) < np.mean(lossRD[del_type == 1]):
            homoRD = np.mean(lossRD[del_type == 0])
            hemiRD = np.mean(lossRD[del_type == 1])
            for i in range(len(CN)):
                if CNVtype[i] == 0:
                    CN[i] = 0
                elif CNVtype[i] == 1:
                    CN[i] = 1
        else:
            hemiRD = np.mean(lossRD[del_type == 0])
            homoRD = np.mean(lossRD[del_type == 1])
            for i in range(len(CN)):
                if CNVtype[i] == 1:
                    CN[i] = 0
                elif CNVtype[i] == 0:
                    CN[i] = 1
        purity = 2 * (homoRD - hemiRD) / (homoRD - 2 * hemiRD)

        for i in range(len(CNVtype)):
            if CNVtype[i] == 2:
                CN[i] = int(2*CNVRD[i] / (mode * purity) - 2 * (1-purity) / purity)
    return CN


def boxplot(scores):
    four = pd.Series(scores).describe()
    Q1 = four['25%']
    Q3 = four['75%']
    IQR = Q3 - Q1
    upper = Q3 + 0.5 * IQR
    lower = Q1 - 0.5 * IQR
    return upper


def combiningCNV(seg_chr, seg_start, seg_end, seg_count, scores, upper, mode):

    index = scores > upper
    CNV_chr = seg_chr[index]
    CNVstart = seg_start[index]
    CNVend = seg_end[index]
    CNVRD = seg_count[index]

    type = np.full(len(CNVRD), 1)
    for i in range(len(CNVRD)):
        if CNVRD[i] > mode:
            type[i] = 2

    for i in range(len(CNVRD) - 1):
        if CNVend[i] + 1 == CNVstart[i + 1] and type[i] == type[i + 1]:
            CNVstart[i + 1] = CNVstart[i]
            type[i] = 0

    index = type != 0
    CNVRD = CNVRD[index]
    CNV_chr = CNV_chr[index]
    CNVstart = CNVstart[index]
    CNVend = CNVend[index]
    CNVtype = type[index]

    return CNV_chr, CNVstart, CNVend, CNVRD, CNVtype

starttime = datetime.datetime.now()
# get params
bam_path = sys.argv[1]
bam = sys.argv[2]
refpath = sys.argv[3]
outpath = sys.argv[4]
segpath = sys.argv[5]
binSize = int(sys.argv[6])
col = int(sys.argv[7])
    
p_value_file = outpath + '/' + bam + ".score.txt"
outfile = outpath + '/' + bam + ".result.txt"

    
ref = [[] for i in range(23)]
refList = read_bam_file(bam_path + bam)
for i in range(len(refList)):
    chr = refList[i]
    chr_num = chr.strip('chr')
    if chr_num.isdigit():
        chr_num = int(chr_num)
        reference = refpath + '/chr' + str(chr_num) + '.fa'
        ref = read_ref_file(reference, chr_num, ref)
    
chrLen = np.full(23, 0)
for i in range(1, 23):
    chrLen[i] = len(ref[i])
RDList, PosList, chrList = Binning(ref, binSize, chrLen, bam_path + bam)
all_chr = []
all_RD = []
all_start = []
all_end = []
modeList = np.full(len(chrList), 0.0)
for i in range(len(chrList)):
    print("analyse " + str(chrList[i]))
    RD = np.array(RDList[i][0])
    pos = np.array(PosList[i][0])
    numbin = len(RD)
    modeList[i] = modeRD(RD)
    scalRD = scaling_RD(RD, modeList[i])
    #plot(pos, RD)
    print("segment count...")
    v = robjects.FloatVector(scalRD)
    m = robjects.r['matrix'](v, ncol=col)
    robjects.r.source("CBS_data.R")
    robjects.r.CBS_data(m, segpath)
    
    num_col = int(numbin / col) + 1
    seg_start, seg_end, seg_count, seg_len = Read_seg_file(num_col, numbin, segpath)
    seg_count = np.array(seg_count)
    
    seg_count = seg_count[:-1]
    seg_start = seg_start[:-1]
    seg_end = seg_end[:-1]
    
    seg_count, seg_start, seg_end = seg_RD(RD, pos, seg_start, seg_end, seg_count, binSize)
    all_RD.extend(seg_count)
    all_start.extend(seg_start)
    all_end.extend(seg_end)
    all_chr.extend(chrList[i] for j in range(len(seg_count)))
    
all_chr = np.array(all_chr)
all_start = np.array(all_start)
all_end = np.array(all_end)
all_RD = np.array(all_RD)
for i in range(len(all_RD)):
    if np.isnan(all_RD[i]).any():
        all_RD[i] = (all_RD[i - 1] + all_RD[i + 1]) / 2
    
# COPOD_CNV
print("calculating scores...")
# train COPOD detector
clf_name = 'COPOD'
clf = COPOD()
    
# you could try parallel version as well.
all_RD = all_RD.reshape(-1,1)
clf.fit(all_RD)
#scores = clf.decision_function(all_RD)
pred = clf.labels_  # binary labels (0: inliers, 1: outliers)
scores = clf.decision_scores_  # raw outlier scores
print("all_RD")
print(all_RD)
print("scores:")
print(scores)
mode = np.mean(modeList)
print("mode:")
print(mode)
#print("pred:")
print(pred)
    

Write_data_file(all_chr, all_start, all_end, all_RD, scores, p_value_file)
upper = boxplot(scores)
CNV_chr, CNVstart, CNVend, CNVRD, CNVtype = combiningCNV(all_chr, all_start, all_end, all_RD, scores, upper, mode)
CN = calculating_CN(mode, CNVRD, CNVtype)
Write_CNV_File(CNV_chr, CNVstart, CNVend, CNVtype, CN, outfile)
    
    
endtime = datetime.datetime.now()
print("running time: " + str((endtime - starttime).seconds) + " seconds")
