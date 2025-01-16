#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 17:09:24 2022

@author: yyj
"""
import copod_preprocess

if __name__ == "__main__":
    
    binSize = 1000
    col = 50
    path = 'simu_chr21_0.2_4x'
    bam_path = '../CNV_data/simu_chr21_0.2_4x/'
    bam = 'sim1_4_4100_read.sort.bam'
    refpath ='../CNV_data'   
    outpath = '../CNV_data/simu_chr21_0.2_4x/result/COPOD_step1'
    segpath = outpath +str("/seg_sim1")
    params = (bam_path, bam, refpath, outpath, segpath, binSize, col)
    copod_preprocess.main(params)


