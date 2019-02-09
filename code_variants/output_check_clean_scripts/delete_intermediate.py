#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Delete intermediate files.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: check HISAT2, STAR, GATK and other outputs for completion

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants
module load python/3.5.0 py_packages/3.5
python

"""

import os

from var_call_class import bam_gatk

home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
os.chdir(home_dir)

with open(home_dir + 'metadata/fq_subdir_prefix_list.txt', 'r') as in_f:
    all_f = [i.strip() for i in in_f]

aligner = 'hisat2'  # 'star' hisat2

# f_list = [bam_i.clean_sam, bam_i.sorted_bam, bam_i.dedup_bam,
#           bam_i.split_trim_bam, bam_i.bqsr_bam, bam_i.vcf_nofilter]
del_true = False
done_id_ct, del_ct, all_ct = 0, 0, 0
for bam_name_i in all_f:
    bam_i = bam_gatk(bam_name_i, home_dir, aligner=aligner)
    if os.path.exists(bam_i.vcf):
        # keep the compressed BAM for HISAT2 instead of SAM
        if aligner == 'hisat2':
            bam_list = [bam_i.in_sam]
        else:
            bam_list = [bam_i.sorted_bam]
        bam_list.extend([bam_i.dedup_bam, bam_i.split_trim_bam])
        # delete the bam and bai
        for bam in bam_list:
            if os.path.exists(bam):
                print('deleting', bam)
                if del_true:
                    os.remove(bam)
                del_ct += 1
            bai = bam[:-1] + 'i'
            if os.path.exists(bai):
                print('deleting', bai)
                if del_true:
                    os.remove(bai)
                del_ct += 1
            # else:
            #     print('Does not exist, check why', bai)
            #     raise ValueError
        done_id_ct += 1
    all_ct += 1

print(del_ct, done_id_ct, all_ct)


"""
"""
