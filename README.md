# COPCNVBD
An Integrated Approach for Somatic Copy Number Variation and Breakpoint Detection using Next-generation Sequencing Data
v1.0

Yaoyao Li
Hangzhou Institute of Technology, Xidian University, 
Hangzhou, China

liyaoyao@xidian.edu.cn

###############################################################################



******************************************************************************
INTRODUCTION
******************************************************************************


Description
===========

COPCNVBD is a python application for somatic copy number variation and breakpoint 
detection from whole genome sequencing data.

Requirements
============

Python 3.8, or newer, is recommended. 
pysam 0.22.1
pandas 2.2.3
numpy 2.0.2
scikit-learn 1.5.2
rpy2 3.5.8
cvxpy
R 4.2.1 or newer

Program descriptions:
====================
before run this algorithm, you should prepare the file of referance (*.fa), the bam file of tested sample(*.bam)

# step1 use the improved COPOD algorithm to acquire the approximate CNV regions
"""
set the input parameters:  
	refpath: the path of reference data
 	outpath: the path of outfile by this algorithm
 	segpath: the path of segfile by segmentation algorithm 
 	binsize: the size of non-overlapping window in steps of preprocessing to divide the genome into segments.
 	col: the size of window when using total variation methods.

eg: 	path = 'simu_chr21_0.2_4x'
        bam_name = '_4_4100_'
        bam_path = '../CNV_data/'+path + '/'
        bam = 'sim' + num + bam_name+'read.sort.bam'
        refpath ='../CNV_data'
        
        outpath = '../CNV_data/'+path+'/result/COPOD_step1'
        segpath = outpath +str("/seg_")+'sim' + num
        p_value_file = outpath + '/' + bam + ".score.txt"
        outfile = outpath + '/' + bam + ".result.txt"
"""

"""
the columns of result file of approximate CNV region represent: 
	chromosome name, 
	start position of the CNV region, 
	end position of the CNV region, 
	the type of CNV region, 
	the copy number
"""
	
run python copod_main.py



# step2: use the PEM information and image edge detection algorithm to identify CNV breakpoint
"""
input parameters:
	refpath: the path of reference data
 	outpath: the path of outfile by this algorithm
	step1_path: the path of approximate CNV region profile

"""
eg: 
	path = 'simu_chr21_0.2_4x'
        bam_name = '_4_4100_'
        bam_path = '../CNV_data/' + path + '/'
        bam = 'sim' + num + bam_name + 'read.sort.bam'
        refpath = '../CNV_data/'
        outpath = '../CNV_data/' + path + '/result/COPOD_step2/'
        step1_path = '../CNV_data/' + path + '/result/COPOD_step1/'

"""
"""
the columns of result file of approximate CNV region represent: 
	chromosome name, 
	start position of the CNV region, 
	end position of the CNV region, 
	the type of CNV region
"""

run python copcnvbd_main.py




