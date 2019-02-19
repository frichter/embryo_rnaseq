#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Identify regions with enough info for variant calling.

Felix Richter
felix.richter@icahn.mssm.edu
2/9/2019
Description: Identify regions with enough info for variant calling

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants
module load gatk/3.6-0
module load picard/2.7.1
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
python

EXAMPLE RUNS:
ID="/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/75888_C4_THS_014_BxE8_2_28_17_S18_L004/75888_C4_THS_014_BxE8_2_28_17_S18_L004"
time python callable_regions.py --aligner "star" --bam $ID
time python callable_regions.py --aligner "hisat2" --bam $ID

"""


import os
import argparse

from var_call_class import bam_gatk


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--aligner', choices=["star", "hisat2"], default=None)
    parser.add_argument('--bam', type=str, help='BAM prefix')
    args = parser.parse_args()
    return args


def prepare_known_regions_files(home_dir):
    """Prepare known regions files."""
    known_dict = {'clinvar_P': 'clinvar_hg38_P_sorted.bed',
                  'clinvar_P_LP': 'clinvar_hg38_P_LP_sorted.bed',
                  'exons': 'exon_locs_hg38_sorted.bed',
                  'genome': 'grch38_ens94_sorted.bed'}
    for known_folder, known_f in known_dict.items():
        known_dict[known_folder] = home_dir + 'known_region/' + known_f


def main():
    """Run GATK."""
    args = _get_args()
    home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
    os.chdir(home_dir)
    bam_i = bam_gatk(args.bam, home_dir, aligner=args.aligner)
    if not os.path.exists(bam_i.bqsr_bam):
        return 'Final GATK BAM not ready ' + bam_i.bqsr_bam
    bam_i.run_callable_loci_gatk()
    return 'Callable loci done for ' + args.bam


if __name__ == "__main__":
    done_msg = main()
    print(done_msg)

"""Resources: GATK callableloci or mosdepth

https://github.com/brentp/mosdepth
https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_coverage_CallableLoci.php

Function: filter for callable regions, pipe to new directory tree

import os

from var_call_class import bam_gatk
from callable_class import call_loci

home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
os.chdir(home_dir)
bam_f = ('/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/' +
         '75888_C4_THS_014_BxE8_2_28_17_S18_L004/75888_C4_THS_014_BxE8_2_28_17_S18_L004')
bam_i = bam_gatk(bam_f, home_dir, aligner='hisat2')
if not os.path.exists(bam_i.bqsr_bam):
    print('Final GATK BAM not ready ' + bam_i.bqsr_bam)


bam_i.run_callable_loci_gatk()
bam_i.id

# create intersections and unions of callable regions
call_i = call_loci(bam_i.id, home_dir)
call_i.subset_callable_loop()
call_i.intersect_callable()
call_i.union_callable()
# overlap with known regions
for known_folder, known_f in known_dict.items():
    print(known_f)
    call_i.intersect_w_known_loci(known_f, known_folder)
# calculate lengths
len_dict_i = call_i.calc_all_bed_lengths()
# write length dictionary to summary file in top-level directory
# (as final output of this entire pipeline)
len_loc = call_i.subdir[:-1] + '_lengths.txt'
with open(len_loc, 'w') as out_f:
    for k, v in len_dict_i.items():
        _ = out_f.write('\t'.join([k, v]) + '\n')


"""
#
