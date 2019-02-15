#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Reorganize files by moving to sample-specific subdirectories.

Felix Richter
felix.richter@icahn.mssm.edu
2/9/2019

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants
module load python/3.5.0 py_packages/3.5
python

"""


import os
import glob
import re
import subprocess

from var_call_class import bam_gatk

home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
os.chdir(home_dir)

with open(home_dir + 'metadata/fq_subdir_prefix_list.txt', 'r') as in_f:
    all_f = [i.strip() for i in in_f]

"""Move intermediate files to sub-directory
"""
with open(home_dir + 'metadata/fq_prefix_list.txt', 'r') as in_f:
    all_f = [i.strip() for i in in_f]

both_done_ct = 0
for bam_name_i in all_f:
    bam_i_star = bam_gatk(bam_name_i, home_dir, aligner='star')
    bam_i_hisat2 = bam_gatk(bam_name_i, home_dir, aligner='hisat2')
    star_done = os.path.exists(bam_i_star.vcf)
    hisat2_done = os.path.exists(bam_i_hisat2.vcf)
    # if not star_done:
    #     continue
    # if not hisat2_done:
    #     continue
    new_dir = re.sub('_star', '', bam_i_star.prefix)
    id = re.sub('.*ASTQ/|.*trim_q20/', '', new_dir)
    f_glob = glob.iglob(new_dir + '*')
    if not os.path.exists(new_dir):
        print('making', new_dir)
        # os.mkdir(new_dir)
    for f_i in f_glob:
        if f_i == new_dir:
            print(f_i)
            print('skipping directory')
            continue
        if f_i == bam_i_hisat2.vcf:
            # print(f_i)
            continue
        if f_i == bam_i_star.vcf:
            # print(f_i)
            continue
        new_f_i = re.sub(id, id + '/' + id, f_i)
        mv_cmd = 'mv {} {}'.format(f_i, new_f_i)
        print(mv_cmd)
        # subprocess.call(mv_cmd, shell=True)
    both_done_ct += 1

"""Fixing a mistake: moving the unfinished data back into top-level."""
for bam_name_i in all_f:
    bam_i_star = bam_gatk(bam_name_i, home_dir, aligner='star')
    bam_i_hisat2 = bam_gatk(bam_name_i, home_dir, aligner='hisat2')
    new_dir = re.sub('_star', '', bam_i_star.prefix)
    id = re.sub('.*ASTQ/|.*trim_q20/', '', new_dir)
    vcf_sub_f_hisat2 = re.sub(id, id + '/' + id, bam_i_hisat2.vcf)
    vcf_sub_f_star = re.sub(id, id + '/' + id, bam_i_star.vcf)
    if os.path.exists(vcf_sub_f_hisat2):
        mv_cmd = 'mv {} {}'.format(vcf_sub_f_hisat2, bam_i_hisat2.vcf)
        print(mv_cmd)
        subprocess.call(mv_cmd, shell=True)
    if os.path.exists(vcf_sub_f_star):
        mv_cmd = 'mv {} {}'.format(vcf_sub_f_star, bam_i_star.vcf)
        print(mv_cmd)
        subprocess.call(mv_cmd, shell=True)

#
