#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Embryo data variant calling with GATK and samtools mpileup.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: variant calling algorithms!

module load gatk/3.6-0
module load picard/2.7.1
module load python/3.5.0 py_packages/3.5
python

"""

import re
import subprocess
import os

"""Running GATK3 RNAseq variant calling
https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

Running GATK4
https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels

"""

os.chdir('/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/')
prefix = 'FASTQ/trim_q20/75888_C4_THS_014_BxE8_2_28_17_S18_L004_hisat2'
id = re.sub('.*ASTQ/|.*trim_q20/', '', prefix)

# Add read groups, sort, mark duplicates, and create index
rg_cmd = ('time java -jar $PICARD AddOrReplaceReadGroups I={} ' +
          'O={} SO=coordinate RGID=id RGLB=nextera ' +
          'RGPL=ilmn RGPU=machine RGSM={}').format(
    prefix + '.sam', prefix + '_rg_sorted.bam', id)
print(rg_cmd)
subprocess.call(rg_cmd, shell=True)
# time 60 mins
# mark MarkDuplicates docs:
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php
dup_cmd = ('time java -jar $PICARD MarkDuplicates I={} ' +
           'O={}_dedupped.bam CREATE_INDEX=true ' +
           'VALIDATION_STRINGENCY=SILENT ' +
           'M={}_gatk3_md_metrics.txt').format(
    prefix + '_rg_sorted.bam', prefix, prefix)
print(dup_cmd)
subprocess.call(dup_cmd, shell=True)

# Split'N'Trim and reassign mapping qualities
ref_fa = ('/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/' +
          'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
split_cmd = ('java -jar $GATK_JAR -T SplitNCigarReads -R {} ' +
             '-I {}_dedupped.bam -o {}_split.bam -rf ' +
             'ReassignOneMappingQuality -RMQF 255 -RMQT 60 ' +
             '-U ALLOW_N_CIGAR_READS').format(
    ref_fa, prefix, prefix)
print(split_cmd)
subprocess.call(split_cmd, shell=True)

# Indel Realignment (optional)

# 5. Base Recalibration

# 6. Variant calling
hc_cmd = ('java -jar $GATK_JAR -T HaplotypeCaller -R {}' +
          '-I {}_split.bam -dontUseSoftClippedBases -stand_call_conf 20.0' +
          '-o {}_gatk3.vcf').format(
    ref_fa, prefix, prefix)
print(hc_cmd)
# subprocess.call(hc_cmd, shell=True)

# 7. Variant filtering
filter_cmd = ('java -jar $GATK_JAR -T VariantFiltration ' +
              '-R {} -V {}_gatk3.vcf -window 35 -cluster 3 ' +
              '-filterName FS -filter "FS > 30.0" -filterName QD ' +
              '-filter "QD < 2.0" -o {}_gatk3_filtered.vcf ').format(
    ref_fa, prefix, prefix)
print(filter_cmd)
# subprocess.call(filter_cmd, shell=True)

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
