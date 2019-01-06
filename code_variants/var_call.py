#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Embryo data variant calling with GATK and samtools mpileup.

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
bam_i.run_picard_rg()
bam_i.run_picard_md()
# Indel Realignment (optional): not doing indels currently
bam_i.run_bqsr()
bam_i.run_gatk_hc()
bam_i.run_gatk_var_filter()

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
