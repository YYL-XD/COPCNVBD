#!/usr/bin/perl

use strict;
use warnings;
use IPC::System::Simple qw(system);

my $python_file1 = "copod_preprocess.py";
my $python_file2 = "copcnvbd_main.py";
my $python_file3 = "sta_boundary.py";

# input param of copodcnv
my $param1_file1 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/";
my $param2_file1 = "/sim1_4_4100_read.sort.bam";
my $param3_file1 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data";
my $param4_file1 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step1";
my $param5_file1 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step1/seg_sim1";
my $param6_file1 = 1000;
my $param7_file1 = 50;

# input param of copodcnvbd
my $param1_file2 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/";
my $param2_file2 = "/sim1_4_4100_read.sort.bam";
my $param3_file2 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/";
my $param4_file2 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step2/";
my $param5_file2 = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step1/";
my @param6_values_file2 = (0.2, 0.25, 0.3, 0.35);  
#my $param6_file2 = 0.25;  


# input param of statistic result
#my $param1_file3  = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step2/sim1_4_4100_read.sort.bam_step2_result_". $param6_file2 ."alpha.txt";
my $param2_file3 = "groundtruth";


print "Running $python_file1 with parameters: $param1_file1 $param2_file1 $param3_file1 $param4_file1 $param5_file1 $param6_file1 $param7_file1\n";
system("python3 $python_file1 $param1_file1 $param2_file1 $param3_file1 $param4_file1 $param5_file1 $param6_file1 $param7_file1") == 0 or die "Failed to run $python_file1: $!";

# run parameter sweep
foreach my $param6 (@param6_values_file2) {
    print "Running $python_file2 with parameters: $param1_file2 $param2_file2 $param3_file2 $param4_file2 $param5_file2 $param6 \n";
    system("python3 $python_file2 $param1_file2 $param2_file2 $param3_file2 $param4_file2 $param5_file2 $param6") == 0 or die "Failed #to run $python_file2 with param6=$param6: $!";
    my $param1_file3  = "/home/lyy/PycharmProjects/pythonProjecttest/CNV_data/simu_chr21_0.2_4x/result/COPOD_step2/sim1_4_4100_read.sort.bam_step2_result_". $param6 ."alpha.txt";
    system("python3 $python_file3 $param1_file3 $param2_file3") == 0 or die "Failed #to run $python_file3: $!";
}


