#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Check outputs of various pipelines.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: check HISAT2, STAR, GATK and other outputs for completion

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants
module load python/3.5.0 py_packages/3.5
python

"""

import glob
import os
import re
import shutil

from var_call_class import bam_gatk

"""Global variables."""
home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
os.chdir(home_dir)
with open(home_dir + 'metadata/fq_subdir_prefix_list.txt', 'r') as in_f:
    all_f = [i.strip() for i in in_f]

"""Hisat2 output check."""
# hs_met_iter = glob.iglob(home_dir + 'FASTQ/*hisat2_metrics.txt')
del_count, uniq_del_ct, tot_count = 0, 0, 0
for f_i in all_f:
    hs_met = f_i + '_hisat2_metrics.txt'
    if not os.path.exists(hs_met):
        print('no metrics for', f_i)
        continue
    if os.stat(hs_met).st_size == 0:
        uniq_del_ct += 1
        hs_f_iter = glob.iglob(re.sub('_metrics.txt', '*', hs_met))
        for hs_f in hs_f_iter:
            print('deleting unfinishing alignment for ' + hs_f)
            # os.remove(hs_f)
            del_count += 1
    tot_count += 1

# number of files deleted:
print(tot_count, uniq_del_ct, del_count)

"""STAR output check."""
true_del = False  # set to True if actually deleting, keep as False if checking
del_count, uniq_del_ct, tot_count, done_ct = 0, 0, 0, 0
for f_i in all_f:
    star_log = f_i + '_starLog.progress.out'
    star_f_final = re.sub('progress.out', 'final.out', star_log)
    if not os.path.exists(star_f_final):
        star_f_iter = glob.iglob(re.sub('Log.progress.out', '*', star_log))
        for star_f in star_f_iter:
            if os.path.isdir(star_f):
                print('deleting directory and contents of ' + star_f)
                if true_del:
                    shutil.rmtree(star_f)
            else:
                print('deleting unfinishing alignment for ' + star_f)
                if true_del:
                    os.remove(star_f)
            del_count += 1
        uniq_del_ct += 1
    else:
        done_ct += 1
    tot_count += 1

# number of files deleted:
print(del_count, uniq_del_ct, done_ct, tot_count)

"""GATK output check.

Overview:
if next file does not exist, delete current file
OR
if filterved VCF does not exist, delete all GATK intermediate files

if stderr has '##### ERROR A USER ERROR has occurred'
grep -c '##### ERROR A USER ERROR has occurred' *stderr
then notify user and examine (and probably delete) the HISAT2 input
Run extra HISAT2 alignments on interactive6
"""

# with open(home_dir + 'metadata/fq_subdir_prefix_list.txt', 'r') as in_f:
#     all_f = [i.strip() for i in in_f]


aligner = 'star'  # 'star' hisat2
true_del = False  # set to True if actually deleting, keep as False if checking
# bam_name_i = all_f[0]
todo_ct, del_ct, all_ct = 0, 0, 0
for bam_name_i in all_f:
    bam_i = bam_gatk(bam_name_i, home_dir, aligner=aligner)
    # Only delete intermediate files if VCF was not completed
    print(bam_i.vcf)
    break
    if not os.path.exists(bam_i.vcf):
        print('todo:', bam_i.clean_sam)
        f_list = [bam_i.clean_sam, bam_i.sorted_bam, bam_i.dedup_bam,
                  bam_i.split_trim_bam, bam_i.bqsr_bam, bam_i.vcf_nofilter]
        # print(f_list)
        for f in f_list:
            # print(f)
            if os.path.exists(f):
                print('deleting', f)
                # the 0-sized files might reflect HISAT2 problematic inputs
                print(os.stat(f).st_size)
                if true_del:
                    os.remove(f)
                del_ct += 1
            bai = f[:-1] + 'i'
            if os.path.exists(bai):
                print('deleting', bai)
                # the 0-sized files might reflect HISAT2 problematic inputs
                print(os.stat(bai).st_size)
                if true_del:
                    os.remove(bai)
                del_ct += 1
        todo_ct += 1
    # else:
    #     print('GATK fully completed for', bam_name_i)
    #     # get the file line of the already completed jobs
    #     print(all_ct)
    all_ct += 1


print(all_ct, todo_ct, del_ct)

"""Update FQ prefix list files.

"""

in_f_loc = home_dir + 'metadata/fq_prefix_list.txt'
out_f_loc = home_dir + 'metadata/fq_subdir_prefix_list.txt'
with open(in_f_loc, 'r') as in_f, open(out_f_loc, 'w') as out_f:
    for line in in_f:
        sub_line = line.strip()
        sub_line += '/' + sub_line.split('/')[-1] + '\n'
        print(sub_line)
        _ = out_f.write(sub_line)


#
