#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run data that failed BSUB.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: Run HISAT2, STAR, GATK and other pipelines that may have
    slipped through the cracks on BSUB


## for HISAT2
module purge
module load fastqc/0.11.7
module load hisat2/2.0.5
module load python/3.5.0

## for STAR
module purge
module load python/3.5.0
module load star/2.6.1d

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/
python

"""

import glob
import re
import subprocess
import os

"""Global variables."""
home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'

"""Files with completed runs."""
hs_met_iter = glob.iglob(home_dir + 'FASTQ/*hisat2_metrics.txt')
hs_f = [re.sub('_hisat2.*', '', i) for i in hs_met_iter]
star_log_iter = glob.iglob(home_dir + 'FASTQ/*starLog.final.out')
star_f = [re.sub('_star.*', '', i) for i in star_log_iter]

"""All files."""
with open(home_dir + 'metadata/fq_prefix_list.txt', 'r') as in_f:
    all_f = [i.strip() for i in in_f]

"""Get to-do file prefixes."""
hs_todo = [i for i in all_f if i not in hs_f]
star_todo = [i for i in all_f if i not in star_f]
len(hs_todo)
len(star_todo)

"""Loop over to-do files."""
os.chdir(home_dir + 'code_variants')
hs_cmd = 'python align_qc.py --hisat2 --fq {} --n 10'
for hs_todo_i in hs_todo:
    subprocess.call(hs_cmd.format(hs_todo_i), shell=True)

star_cmd = 'python align_qc.py --star --fq {} --n 10'
for star_todo_i in star_todo:
    subprocess.call(star_cmd.format(star_todo_i), shell=True)

#
