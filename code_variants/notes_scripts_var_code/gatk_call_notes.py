#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Notes for variant calling from RNAseq in embryo data."""

import re
import subprocess
import os

"""Other variant calling pipelines:
rnacocktail
https://github.com/bioinform/rnacocktail/blob/master/src/run_variant.py

review all options here:
https://bioinform.github.io/rnacocktail/

Running GATK3 RNAseq variant calling
https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

Running GATK4
https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels

"""

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


"""
home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
file_prefix = home_dir + 'FASTQ/95724_C1_THS_024_BxE2_3_13_17_S20_L004'
os.chdir(home_dir)
bam_i = bam_gatk(file_prefix, home_dir, aligner='hisat2')
# bam_i.run_picard_cs()  # used in rnacocktail but purpose unclear so skipping
bam_i.run_picard_rg()  # 12m
bam_i.run_picard_md()  # also 12m
bam_i.run_gatk_split_trim()  # 32mins
# Indel Realignment (optional): not doing indels currently
bam_i.run_bqsr()
bam_i.run_gatk_hc()
bam_i.run_gatk_var_filter()

# trying with STAR
home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
file_prefix = home_dir + 'FASTQ/75888_C4_THS_014_BxE8_2_28_17_S18_L004'
os.chdir(home_dir)
bam_i = bam_gatk(file_prefix, home_dir, aligner='star')
bam_i.prefix
# bam_i.run_picard_rg()
star_bam = file_prefix + '_starAligned.sortedByCoord.out.bam'
rg_cmd = ('time java -Djava.io.tmpdir={} ' +
          '-jar $PICARD AddOrReplaceReadGroups I={} ' +
          'O={} SO=coordinate RGID=FelixRichter RGLB=Nextera ' +
          'RGPL=ILLUMINA RGPU=machine RGSM={}').format(
    bam_i.tmp_dir,
    # use clean_sam if using run_picard_cs, otherwise in_sam:
    star_bam,
    bam_i.sorted_bam, bam_i.id)
print(rg_cmd)
subprocess.call(rg_cmd, shell=True)

bam_i.run_picard_md()  # also 12m
bam_i.run_gatk_split_trim()  # 32mins
# Indel Realignment (optional): not doing indels currently
bam_i.run_bqsr()
bam_i.run_gatk_hc()
bam_i.run_gatk_var_filter()

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
"""
# 5. Base Recalibration

"""

"""1. Analyze patterns of covariation in the sequence dataset."""
ref_gatk_dir = '/sc/orga/projects/chdiTrios/Felix/dbs/gatk_resources/'
bqsr_mk_tbl_cmd = (
    'time java -jar $GATK_JAR -T BaseRecalibrator ' +
    '-R {} -I {} -knownSites {} -knownSites {} -o {}').format(
    ref_fa,
    prefix + '_split.bam',
    ref_gatk_dir + 'dbsnp_146.grch38.sorted.vcf',
    ref_gatk_dir + 'Mills_and_1000G_gold_standard.indels.grch38.sorted.vcf',
    prefix + '.recal_data.table')
print(bqsr_mk_tbl_cmd)
subprocess.call(bqsr_mk_tbl_cmd, shell=True)
"""2. Do a second pass to analyze covariation remaining after recalibration
"""
bqsr_mk_tbl_p2_cmd = (
    'time java -jar $GATK_JAR -T BaseRecalibrator ' +
    '-R {} -I {} -knownSites {} -knownSites {} ' +
    '-BQSR {} -o {}').format(
    ref_fa,
    prefix + '_split.bam',
    ref_gatk_dir + 'dbsnp_146.grch38.sorted.vcf',
    ref_gatk_dir + 'Mills_and_1000G_gold_standard.indels.grch38.sorted.vcf',
    prefix + '.recal_data.table',
    prefix + '.recal_data_pass2.table')
print(bqsr_mk_tbl_p2_cmd)
subprocess.call(bqsr_mk_tbl_p2_cmd, shell=True)
"""3. Generate before/after plots
java -jar GenomeAnalysisTK.jar \
    -T AnalyzeCovariates \
    -R reference.fa \
    -L 20 \
    -before recal_data.table \
    -after post_recal_data.table \
    -plots recalibration_plots.pdf
"""
# bqsr_cmd = ('time java -jar $GATK_JAR -T PrintReads ' +
#             '-R {} -I {} -BQSR {}_bqsr.grp -o {}').format(
#     ref_fa, prefix + '_split.bam', prefix, prefix + '_bqsr.bam')
# print(bqsr_cmd)
# subprocess.call(bqsr_cmd, shell=True)

"""4. Apply the recalibration to your sequence data
java -jar GenomeAnalysisTK.jar \
    -T PrintReads \
    -R reference.fa \
    -I input_reads.bam \
    -L 20 \
    -BQSR recal_data.table \
    -o recal_reads.bam
"""
bqsr_cmd = ('time java -jar $GATK_JAR -T PrintReads ' +
            '-R {} -I {} -BQSR {} -o {}').format(
    ref_fa, prefix + '_split.bam',
    prefix + '.recal_data_pass2.table',
    prefix + '_bqsr.bam')
print(bqsr_cmd)
subprocess.call(bqsr_cmd, shell=True)
# 6. Variant calling
tmp_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/tmp_dir/'
hc_cmd = ('time java -Djava.io.tmpdir={} -jar $GATK_JAR ' +
          '-T HaplotypeCaller -R {} ' +
          '-I {}_split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 ' +
          '-dbsnp {} ' +
          # dbsnp uses known sites for variant annotation
          '-o {}_gatk3.vcf').format(
    tmp_dir, ref_fa, prefix,
    ref_gatk_dir + 'dbsnp_146.grch38.sorted.vcf',
    prefix)
print(hc_cmd)
subprocess.call(hc_cmd, shell=True)
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
