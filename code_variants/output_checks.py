#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Check outputs of various pipelines.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: check HISAT2, STAR, GATK and other outputs for completion

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/
module load python/3.5.0 py_packages/3.5
python

"""

import glob
import os
import re
import shutil

"""Global variables."""
home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'

"""Hisat2 output check."""
hs_met_iter = glob.iglob(home_dir + 'FASTQ/*hisat2_metrics.txt')
count, uniq_ct = 0, 0
for hs_met in hs_met_iter:
    if os.stat(hs_met).st_size == 0:
        uniq_ct += 1
        hs_f_iter = glob.iglob(re.sub('_metrics.txt', '*', hs_met))
        for hs_f in hs_f_iter:
            print('deleting unfinishing alignment for ' + hs_f)
            # os.remove(hs_f)
            count += 1

# number of files deleted:
print(uniq_ct, count)

"""STAR output check."""
star_log_iter = glob.iglob(home_dir + 'FASTQ/*starLog.progress.out')
count, uniq_ct, tot_count = 0, 0, 0
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
            count += 1
        uniq_ct += 1
    else:
        print(star_f_final)
    tot_count += 1

# number of files deleted:
print(count, uniq_ct, tot_count)

"""GATK output check."""

#
