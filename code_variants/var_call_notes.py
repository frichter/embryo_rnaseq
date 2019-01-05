#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Notes for variant calling from RNAseq in embryo data."""

import re
import subprocess
import os

"""GATK quick start guide
https://software.broadinstitute.org/gatk/documentation/quickstart?v=3
java -version
java -jar $GATK_JAR -h
java -jar $PICARD -h
"""

"""PREPARING THE FASTA INDEX:
module load gatk/3.6-0
module load picard/2.7.1
module load samtools/1.9

ENSEMBL_DIR="/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/"
ENSEMBL_FA=$ENSEMBL_DIR"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GATK_DICT=$ENSEMBL_DIR"Homo_sapiens.GRCh38.dna.primary_assembly.dict"

cd $ENSEMBL_DIR
time samtools faidx $ENSEMBL_FA
time java -jar $PICARD CreateSequenceDictionary R=$ENSEMBL_FA O=$GATK_DICT
# Tutorial:
https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
# Reference:
https://software.broadinstitute.org/gatk/documentation/article.php?id=1601

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
subprocess.call(dup_cmd, shell=True)  # 64 minutes

# Split'N'Trim and reassign mapping qualities
ref_fa = ('/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/' +
          'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
split_cmd = ('time java -jar $GATK_JAR -T SplitNCigarReads -R {} ' +
             '-I {}_dedupped.bam -o {}_split.bam -rf ' +
             'ReassignOneMappingQuality -RMQF 255 -RMQT 60 ' +
             '-U ALLOW_N_CIGAR_READS').format(
    ref_fa, prefix, prefix)
print(split_cmd)
subprocess.call(split_cmd, shell=True)

# Indel Realignment (optional)
# 5. Base Recalibration
bqsr_cmd = ('time java -jar $GATK_JAR -T PrintReads ' +
            '-R {} -I {} -BQSR {}_bqsr.grp -o {}').format(
    ref_fa, prefix + '_split.bam', prefix, prefix + '_bqsr.bam')
print(bqsr_cmd)
subprocess.call(bqsr_cmd, shell=True)
# 6. Variant calling
hc_cmd = ('time java -jar $GATK_JAR -T HaplotypeCaller -R {} ' +
          '-I {}_bqsr.bam -dontUseSoftClippedBases -stand_call_conf 20.0 ' +
          '-o {}_gatk3.vcf').format(
    ref_fa, prefix, prefix)
print(hc_cmd)
# subprocess.call(hc_cmd, shell=True)
# 7. Variant filtering
filter_cmd = ('time java -jar $GATK_JAR -T VariantFiltration ' +
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
