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

"""Hisat2 output check."""
hs_met_iter = glob.iglob(home_dir + 'FASTQ/*hisat2_metrics.txt')
del_count, uniq_del_ct, tot_count = 0, 0, 0
for hs_met in hs_met_iter:
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
star_log_iter = glob.iglob(home_dir + 'FASTQ/*starLog.progress.out')
del_count, uniq_del_ct, tot_count = 0, 0, 0
for star_log in star_log_iter:
    star_f_final = re.sub('progress.out', 'final.out', star_log)
    if not os.path.exists(star_f_final):
        star_f_iter = glob.iglob(re.sub('Log.progress.out', '*', star_log))
        for star_f in star_f_iter:
            if os.path.isdir(star_f):
                print('deleting directory and contents of ' + star_f)
                shutil.rmtree(star_f)
            else:
                print('deleting unfinishing alignment for ' + star_f)
                os.remove(star_f)
            del_count += 1
        uniq_del_ct += 1
    # else:
    #     print(star_f_final)
    tot_count += 1

# number of files deleted:
print(del_count, uniq_del_ct, tot_count)

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

with open(home_dir + 'metadata/fq_prefix_list.txt', 'r') as in_f:
    all_f = [i.strip() for i in in_f]

aligner = 'hisat2'  # 'star'
# bam_name_i = all_f[0]
todo_ct, del_ct, all_ct = 0, 0, 0
for bam_name_i in all_f:
    bam_i = bam_gatk(bam_name_i, home_dir, aligner=aligner)
    if not os.path.exists(bam_i.vcf):
        f_list = [bam_i.clean_sam, bam_i.sorted_bam, bam_i.dedup_bam,
                  bam_i.split_trim_bam, bam_i.bqsr_bam, bam_i.vcf_nofilter]
        # print(f_list)
        for f in f_list:
            # print(f)
            if os.path.exists(f):
                print('deleting', f)
                # the 0-sized files might reflect HISAT2 problematic inputs
                print(os.stat(f).st_size)
                # os.remove(f)
                del_ct += 1
        # break
        todo_ct += 1
    else:
        print('GATK fully completed for', bam_name_i)
        # get the file line of the already completed jobs so they aren't rerun
        print(all_ct)
    all_ct += 1


print(todo_ct, del_ct)
"""
Failed files
Exception in thread "main" htsjdk.samtools.SAMFormatException:
Error parsing text SAM file. Not enough fields

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/
tail 90099_C4_THS_010_BxE10_3_1_17_S2_L001_hisat2.sam
Line 96058617
SND00311:234:CC56TANXX:1:2305:13307:65206

rm 90099_C4_THS_010_BxE10_3_1_17_S2_L001_hisat2_rg_sorted.bam
# STAR number of input reads 45255155

"""

#
