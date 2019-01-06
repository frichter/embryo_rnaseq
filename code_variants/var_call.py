#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Embryo data variant calling with GATK and samtools mpileup.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: variant calling algorithms!

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants
module load gatk/3.6-0
module load picard/2.7.1
module load python/3.5.0 py_packages/3.5
python

"""

import os

from var_call_class import bam_gatk

home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
file_prefix = home_dir + 'FASTQ/95724_C1_THS_024_BxE2_3_13_17_S20_L004'
os.chdir(home_dir)
bam_i = bam_gatk(file_prefix, home_dir, aligner='hisat2')
bam_i.run_picard_rg()  # 12m
bam_i.run_picard_md()
bam_i.run_gatk_split_trim()
# Indel Realignment (optional): not doing indels currently
bam_i.run_bqsr()
bam_i.run_gatk_hc()
bam_i.run_gatk_var_filter()

"""Future GATK directions (move to gatk_call_notes.py when implementing):

- Consider parallelism or other methods of speeding up
https://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
- Consider adding in the -Xmx and -Xms options
- Consider other HaplotypeCaller arguments (e.g., -maxAltAlleles)
https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

"""

#

"""Running Samtools

# notes from HISAT2 docs
samtools view -bS eg2.sam > eg2.bam
samtools sort eg2.bam -o eg2.sorted.bam
samtools mpileup -uf $HISAT2_HOME/example/reference/22_20-21M.fa
eg2.sorted.bam | bcftools view -bvcg - > eg2.raw.bcf

# URLS
Run samtools mpileup, http://www.htslib.org/doc/samtools.html
Then apparently run bcftools mpileup or call,
http://www.htslib.org/doc/bcftools.html#mpileup
Example workflow:
http://www.htslib.org/workflow/#mapping_to_variant

"""
