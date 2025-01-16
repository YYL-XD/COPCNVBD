# -*- coding: utf-8 -*-
"""

@author: 24240
"""
import datetime
import sys

import pysam
import pandas as pd
import numpy as np
import cvxpy as cp
from scipy.stats import alpha


# 输出每个read的ID 和 开始和终止位置
def func1(bam_path, bam, outpath):
    bam_file = pysam.AlignmentFile(bam_path + bam, "rb")
    query_positions = {}
    # for read in bam_file.fetch():
    for read in bam_file:
        query_positions[read.query_name] = (read.reference_start, read.reference_end)

    with open(outpath + bam + "_query_positions.txt", "w") as f:
        for query_name in sorted(query_positions.keys()):
            start, end = query_positions[query_name]
            f.write(f"{query_name}\t{start}\t{end}\n")

'''
# 读取数据
with open("query_positions.txt", "r") as f:
    lines = f.readlines()
# 将每行数据拆分为列表
data = [line.strip().split('\t') for line in lines]
'''
def func2(outpath,bam):
    data = pd.read_csv(outpath+bam+"_query_positions.txt", sep='\t',header = None,names=['readID','start', 'end'])
    data['end'] = data['end'].astype('Int64')

    # 将相邻的单数行和偶数行两行合并成为一行，并生成DataFrame
    df = pd.DataFrame([[data['readID'][i], data['readID'][i+1], data['start'][i],data['start'][i+1], data['end'][i], data['end'][i+1]] for i in range(0, len(data), 2)],
                    columns=['readID_odd', 'readID_even', 'start_odd', 'end_odd', 'start_even', 'end_even'])

    #df.fillna(111111111111111111111111111111111111111111111111111111111111111,inplace=True)
    # 生成两列，分别为后四列的最小值和最大值
    df['min_pos'] = df[['start_odd', 'end_odd', 'start_even', 'end_even']].min(axis=1)
    df['max_pos'] = df[['start_odd', 'end_odd', 'start_even', 'end_even']].max(axis=1)

    df.to_csv(outpath+bam+"_merged_data.txt", sep="\t", index=False,header=False)


def func3(step1_path, bam, outpath):
    # 读取第一个txt文件
    with open(step1_path + bam + ".result.txt") as f1:
        lines1 = f1.readlines()

    # 读取第二个txt文件
    with open(outpath + bam + "_merged_data.txt") as f2:
        lines2 = f2.readlines()

    # 新的列表
    new_lines = []

    # 遍历第一个txt文件每一行
    for line1 in lines1:
        # 提取第一个txt文件的起始位置和终止位置
        chr1, start1, end1, type1, value1 = line1.strip().split("\t")
        start1 = int(start1)
        end1 = int(end1)
        min_pos = start1
        max_pos = end1

        # 遍历第二个txt文件每一行
        for line2 in lines2:
            # 提取第二个txt文件的起始位置和终止位置
            _, _, _, _, _, _, start2, end2 = line2.strip().split("\t")
            start2 = int(start2)
            end2 = int(end2)

            # 判断第二个txt文件的起始位置和终止位置是否与第一个txt文件有重叠
            if start2 <= end1 and end2 >= start1:
                # 找到最小位置和最大位置
                min_pos = min(min_pos, start2)
                max_pos = max(max_pos, end2)

        # 将新的两列添加到line1中
        new_line = line1.strip() + "\t" + str(min_pos) + "\t" + str(max_pos) + "\n"
        # 将新的行添加到新的列表中
        new_lines.append(new_line)

    # 将新的列表写入到新的txt文件中
    with open(outpath + bam + "_hx_result.txt", "w") as f:
        f.writelines(new_lines)


def readFasta(filename):
    seq = ''
    fread = open(filename)
    # delete row one
    fread.readline()

    line = fread.readline().strip()
    while line:
        seq += line
        line = fread.readline().strip()
    return seq


# function to read readCount file , generate readCount array
def readRd(filename, seqlen):
    #print(seqlen)
    readCount = np.full(seqlen, 0.0)
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                posList = line.positions
                readCount[posList] += 1

    return readCount


def read_step1_result(filename):
    region_array = pd.read_csv(filename,
                               names=['chr_name', 'step1_start', 'step1_end', 'type', 'type_num', 'hx_start', 'hx_end'],
                               header=None, sep='\t')
    region_array = region_array[['chr_name', 'step1_start', 'step1_end', 'type', 'hx_start', 'hx_end']]
    region_array_len = region_array.shape[0]
    jz_region = region_array[['chr_name', 'step1_start', 'step1_end']].values.tolist()
    hx_region = region_array[['chr_name', 'hx_start', 'hx_end']].values.tolist()
    return region_array, region_array_len, jz_region, hx_region


def tv_smooth(y, lam):
    n = len(y)
    x = cp.Variable(n)
    obj = cp.Minimize(0.5 * cp.sum_squares(x - y) + lam * cp.norm(cp.diff(x), 1))
    prob = cp.Problem(obj)
    prob.solve(solver='OSQP')
    return x.value


def find_hx_pos(arr, threshold):
    grad = np.diff(arr)
    indices = np.where(np.abs(grad) > threshold)[0] + 1
    return indices


def find_nearest_boundary(candidates, left_boundary, right_boundary):
    # 如果候选边界为空，则直接使用基准边界
    if len(candidates) == 0:
        return left_boundary, right_boundary

    # 查找最近的左边界
    if left_boundary in candidates:
        cnv_left_boundary = left_boundary
    else:
        left_pointers = [left_boundary - i for i in range(1, len(candidates) + 1)]
        right_pointers = [left_boundary + i for i in range(1, len(candidates) + 1)]
        for l, r in zip(left_pointers, right_pointers):
            if l in candidates:
                cnv_left_boundary = l
                break
            elif r in candidates:
                cnv_left_boundary = r
                break
        else:
            cnv_left_boundary = left_boundary

    # 查找最近的右边界
    if right_boundary in candidates:
        cnv_right_boundary = right_boundary
    else:
        left_pointers = [right_boundary - i for i in range(1, len(candidates) + 1)]
        right_pointers = [right_boundary + i for i in range(1, len(candidates) + 1)]
        for l, r in zip(left_pointers, right_pointers):
            if r in candidates:
                cnv_right_boundary = r
                break
            elif l in candidates:
                cnv_right_boundary = l
                break
        else:
            cnv_right_boundary = right_boundary

    return cnv_left_boundary, cnv_right_boundary


def main(bam, outpath, refpath, bam_path, alpha):
    hx_result = outpath + bam + "_hx_result.txt"
    region_array, region_array_len, jz_region, hx_region = read_step1_result(hx_result)

    binLen = 10
    ref = refpath + 'chr21.fa'
    chrName = ref.split('.')[0]
    rdFile = bam_path + bam
    outputFile = outpath + bam + '_step2_result_' + str(alpha) +'alpha.txt'

    treeNum = 256
    treeSampleNum = 256
    ##alpha = 0.25

    errorNum = 0.005
    seq = readFasta(ref)
    # The length of seq
    seqlen = len(seq)
    print("seqlen:" + str(seqlen))
    rd = readRd(rdFile, seqlen)
    print("rd already finished")
    # fillin n and N position
    rd_mean = np.mean(rd)
    for i in range(seqlen):
        if seq[i] not in ['a', 'A', 't', 'T', 'c', 'C', 'g', 'G']:
            rd[i] = rd_mean
    rd = rd.astype(np.float64)
    print("preprocessing has finished!")

    # hx_rd
    result = [[elem for elem in rd[hx_region[i][1]:hx_region[i][2]]] for i in range(len(hx_region))]
    print("rd_result is saved!")

    # hx_pos
    result_tv = []
    hx_pos = []
    for i in range(len(result)):
        y = np.array(result[i])
        x_smooth = tv_smooth(y, alpha)
        start = int(jz_region[i][1])
        end = int(jz_region[i][2])
        start_left = start - (end - start + 1)
        end_left = start - 1
        rd_left = np.mean(rd[start_left:end_left + 1])
        rd_jz = np.mean(rd[start:end])
        threshold = np.abs(np.diff(np.array([rd_left, rd_jz])))[0]
        #print(threshold)
        hx_pos_i = find_hx_pos(x_smooth, threshold)
        if len(hx_pos_i) == 0:
            print("array is empty")
        else:
            hx_pos_i = hx_pos_i + int(hx_region[i][1]) + 1
        result_tv.append(x_smooth)
        hx_pos.append(hx_pos_i)
    print("hx_pos is following....")
    print("\n")
    #print(hx_pos)

    # CNVbd
    step2_result = region_array[['chr_name', 'type']]
    step2_result['step2_start'] = 0
    step2_result['step2_end'] = 0

    for j in range(len(hx_pos)):
        jz_left = int(jz_region[j][1])
        jz_right = int(jz_region[j][2])
        step2_leftbd, step2_rightbd = find_nearest_boundary(hx_pos[j], jz_left, jz_right)
        # step2_result['step2_start'][j] = step2_leftbd
        step2_result.loc[j, 'step2_start'] = step2_leftbd
        # step2_result['step2_end'][j] = step2_rightbd
        step2_result.loc[j, 'step2_end'] = step2_rightbd

    step2_result = step2_result[['chr_name', 'step2_start', 'step2_end', 'type']]
    step2_result.to_csv(outputFile, sep='\t', index=False)


#get param and run
start_time = datetime.datetime.now()
bam_path = sys.argv[1] #'/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/'
bam = sys.argv[2]  #'sim1_4_4100_read.sort.bam'
refpath = sys.argv[3]#'/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/'
outpath = sys.argv[4]#'/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step2/'
step1_path = sys.argv[5]#'/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step1/'
alpha = float(sys.argv[6])#0.25

func1(bam_path, bam, outpath)
func2(outpath, bam)
func3(step1_path, bam, outpath)
main(bam, outpath, refpath, bam_path, alpha)
end_time = datetime.datetime.now()
print("running time: " + str((end_time - start_time).seconds) + " seconds")

## python copcnvbd_main.py /home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/ /sim1_4_4100_read.sort.bam /home/lyy/PycharmProjects/pythonProjecttest/CNV_data/ /home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step2/ /home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step1/ 0.3