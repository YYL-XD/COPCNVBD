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
    k = 10
    
    for i in range(1,2):
        num = str(i)
        path = 'simu_chr21_0.2_4x'
        bam_name = '_4_4100_'
        bam_path = '../CNV_data/'+path + '/'
        bam = 'sim' + num + bam_name+'read.sort.bam'
        refpath ='../CNV_data'
        
        outpath = '../CNV_data/'+path+'/result/COPOD_step1'
        segpath = outpath +str("/seg_")+'sim' + num
        p_value_file = outpath + '/' + bam + ".score.txt"
        outfile = outpath + '/' + bam + ".result.txt"
        
        params = (bam_path, bam, refpath, outpath, segpath, binSize, col, k)
        copod_preprocess.main(params)


