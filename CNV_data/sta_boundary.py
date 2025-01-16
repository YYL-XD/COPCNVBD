import sys
import numpy as np

outfile = sys.argv[1]   #outfile = "sim1_4_4100_read.sort.bam_step2_result_0.35alpha.txt"
ground_file = sys.argv[2] #"groundtruth"

result_start = []
result_end = []
result_type = []
#result_chrom = []
with open(outfile, 'r') as f:
    next(f);
    for line in f:
        linestr = line.strip()
        linestrlist = linestr.split('\t')
        #result_chrom.append(linestrlist[0])
        result_start.append(int(linestrlist[1]))
        result_end.append(int(linestrlist[2]))
        result_type.append(linestrlist[3])


flag = np.full(len(result_type), 1)
for i in range(len(result_type) - 1):
    if result_end[i] == result_start[i+1] and result_type[i] == result_type[i+1]:
        result_start[i+1] = result_start[i]
        flag[i] = 0
result_start = np.array(result_start)
result_end = np.array(result_end)
result_type = np.array(result_type)
index = flag > 0
result_start = result_start[index]
result_end = result_end[index]
result_type = result_type[index]


truth_start = []
truth_end = []
truth_type = []
#truth_chrom = []
with open(ground_file, 'r') as f:
    for line in f:
        linestr = line.strip('\n')
        linestrlist = linestr.split('\t')
        #truth_chrom.append(linestrlist[0])
        truth_start.append(int(linestrlist[0]))
        truth_end.append(int(linestrlist[1]))
        truth_type.append(linestrlist[3])


count = 0
length = 0
flag = [0 for i in range(len(truth_type))]
for i in range(len(result_type)):
    for j in range(len(truth_type)):
        if (truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]
                and flag[j] == 0):
            length += (abs(result_start[i] - truth_start[j]) + abs(result_end[i] - truth_end[j]))/(
                    truth_end[j] - truth_start[j])
            count += 1
            type = 1
            flag[j] = 1
            break
        elif (truth_start[j] <= result_end[i] <= truth_end[j] and truth_type[j] == result_type[i]
              and flag[j] == 0):
            count += 1
            length += (abs(result_start[i] - truth_start[j]) + abs(result_end[i] - truth_end[j])) / (
                        truth_end[j] - truth_start[j])
            type = 2
            flag[j] = 1
            break
        elif (result_start[i] <= truth_start[j] <= result_end[i] and result_start[i] <= truth_end[j] <= result_end[i]
              and flag[j] == 0):
            count += 1
            length += (abs(result_start[i] - truth_start[j]) + abs(result_end[i] - truth_end[j])) / (
                        truth_end[j] - truth_start[j])
            type = 3
            flag[j] = 1
            break

print("the result file is " +  outfile  + ", the proportion of CNV boundary deviation is：" + str(length/count))
