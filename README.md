# COPCNVBD
An Integrated Approach for Somatic Copy Number Variation and Breakpoint Detection using Next-generation Sequencing Data
v1.0

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

Test data
============
users can download the test data sim1_4_4100_read.sort.bam by the master branch.

Program descriptions:
====================
before run this algorithm, you should prepare the file of referance (*.fa), the bam file of tested sample(*.bam)

# step1 use the improved COPOD algorithm to acquire the approximate CNV regions

## set the input parameters: 
	input_bam_path: the path of input bam file
 	bam_name: the name of bam
	refpath: the path of reference data
 	outpath: the path of outfile by this algorithm
 	segpath: the path of segfile by segmentation algorithm 
 	binsize: the size of non-overlapping window in steps of preprocessing to divide the genome into segments.
 	col: the size of window when using total variation methods.
##

"""
the columns of result file of approximate CNV region represent: 
	chromosome name, 
	start position of the CNV region, 
	end position of the CNV region, 
	the type of CNV region, 
	the copy number
"""
	
run python copod_preprocess.py /../CNV_data/simu_chr21_0.2_4x/ /sim1_4_4100_read.sort.bam /../CNV_data /../CNV_data/simu_chr21_0.2_4x/result/COPOD_step1 /../CNV_data/simu_chr21_0.2_4x/result/COPOD_step1/seg_sim1 1000 50



# step2: use the PEM information and image edge detection algorithm to identify CNV breakpoint
"""
input parameters:
	input_bam_path: the path of input bam file
 	bam_name: the name of bam
	refpath: the path of reference data
 	outpath: the path of outfile by this algorithm
	step1_path: the path of approximate CNV region profile
 	alpha: the penalty parameter in total variation method

"""
"""
the columns of result file of approximate CNV region represent: 
	chromosome name, 
	start position of the CNV region, 
	end position of the CNV region, 
	the type of CNV region
"""

run python copcnvbd_main.py /../CNV_data/simu_chr21_0.2_4x/ /sim1_4_4100_read.sort.bam /../CNV_data/ /../CNV_data/simu_chr21_0.2_4x/result/COPOD_step2 /../CNV_data/simu_chr21_0.2_4x/result/COPOD_step1 0.3

# Note: use a parameter sweep to choose the penalty paramter alpha
run perl copcnvbd_parameter_sweep.pl


