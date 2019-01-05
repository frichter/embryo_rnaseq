#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Embryo data variant calling with GATK and samtools mpileup.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: variant calling algorithms!

env | grep -i 'gatk\|picard'

module load gatk/3.6-0
module load picard/2.7.1

echo $GATK_JAR
echo $PICARD

"""

import re

"""GATK quick start guide
https://software.broadinstitute.org/gatk/documentation/quickstart?v=3
java -version
java -jar $GATK_JAR -h
java -jar $PICARD -h
"""


"""Running GATK3 RNAseq variant calling
https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

Running GATK4
https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels

"""
prefix = 'FASTQ/75888_C4_THS_014_BxE8_2_28_17_S18_L004_hisat2'
id = re.sub('.*ASTQ/', '', prefix)

# Add read groups, sort, mark duplicates, and create index
rg_cmd = ('java -jar $PICARD AddOrReplaceReadGroups I={} ' +
          'O={} SO=coordinate RGID=id RGLB=nextera ' +
          'RGPL=ilmn RGPU=machine RGSM={}').format(
    prefix + '.sam', prefix + '_rg_sorted.bam', id)
print(rg_cmd)
# subprocess.call(rg_cmd, shell=True)
# mark MarkDuplicates docs:
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php
dup_cmd = ('java -jar $PICARD MarkDuplicates I={} ' +
           'O={}_dedupped.bam CREATE_INDEX=true ' +
           'VALIDATION_STRINGENCY=SILENT ' +
           'M={}.metrics').format(
    prefix + '_rg_sorted.bam', prefix, prefix)
print(dup_cmd)
# subprocess.call(dup_cmd, shell=True)

# Split'N'Trim and reassign mapping qualities
ref_fa = 'ref.fasta'
split_cmd = ('java -jar $GATK_JAR -T SplitNCigarReads -R {} ' +
             '-I {}_dedupped.bam -o {}_split.bam -rf ' +
             'ReassignOneMappingQuality -RMQF 255 -RMQT 60 ' +
             '-U ALLOW_N_CIGAR_READS').format(
    ref_fa, prefix, prefix)
print(split_cmd)
# subprocess.call(split_cmd, shell=True)

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

"""RNA cocktail options used (all GATK options have been x-checked)

HaplotypeCaller: -stand call conf 20.0
-stand emit conf 20.0 -A StrandBiasBySample -A StrandAlleleCountsBySampleGATK

VariantFiltration: -window 35 -cluster 3 -filterName FS -filter
"FS >30.0" -filterName QD -filter "QD <2.0"


v3.5-0-g36282e4
(picard 1.129)

samtools mpileup -C50 -d 100000
bcftools filter -s LowQual -e ‘%QUAL<20 —— DP>10000’


"""

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
