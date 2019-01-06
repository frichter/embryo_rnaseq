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

import re
import subprocess
import os

"""Running GATK3 RNAseq variant calling
https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

Running GATK4
https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels

"""


class bam_gatk(object):
    """Create an object and functions for GATK variant calling from a BAM."""

    def __init__(self, file_prefix, home_dir, aligner='hisat2'):
        """Create the BAM object."""
        self.home_dir = home_dir
        self.ref_fa = (
            '/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/' +
            'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
        aligner_types = ['hisat2', 'star']
        if aligner not in aligner_types:
            raise ValueError("Invalid aligner. Expected one: {}".format(
                aligner_types))
        self.prefix = file_prefix + '_' + aligner
        self.id = re.sub('.*ASTQ/|.*trim_q20/', '', self.prefix)
        # intermediate bams from the GATK steps
        self.init_file_names()

    def init_file_names(self):
        """Initialize file names for intermediate and final GATK steps."""
        self.in_sam = self.prefix + '.sam'
        self.sorted_bam = self.prefix + '_rg_sorted.bam'
        self.dedup_bam = self.prefix + '_dedupped.bam'
        self.split_trim_bam = self.prefix + '_split.bam'
        self.bqsr_bam = self.prefix + '_bqsr.bam'
        self.vcf_nofilter = self.prefix + '_gatk3_nofilter.vcf'
        self.vcf = self.prefix + '_gatk3.vcf'

    def run_picard_cmds(self):
        """Run PICARD commands."""
        # Add read groups, sort, mark duplicates, and create index
        if not os.path.exists(self.sorted_bam):
            rg_cmd = ('time java -jar $PICARD AddOrReplaceReadGroups I={} ' +
                      'O={} SO=coordinate RGID=FelixRichter RGLB=Nextera ' +
                      'RGPL=Illumina RGPU=machine RGSM={}').format(
                self.in_sam, self.sorted_bam, self.id)
            print(rg_cmd)
            subprocess.call(rg_cmd, shell=True)
        # time 60 mins
        # mark MarkDuplicates docs:
        # https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php
        if not os.path.exists(self.dedup_bam):
            dup_cmd = ('time java -jar $PICARD MarkDuplicates I={} ' +
                       'O={} CREATE_INDEX=true ' +
                       'VALIDATION_STRINGENCY=SILENT ' +
                       'M={}_gatk3_md_metrics.txt').format(
                self.sorted_bam, self.dedup_bam, self.prefix)
            print(dup_cmd)
            subprocess.call(dup_cmd, shell=True)

    def run_gatk_cmds(self):
        """Run GATK3 commands for variant calling."""
        # Split'N'Trim and reassign mapping qualities
        if not os.path.exists(self.split_trim_bam):
            split_cmd = ('time java -jar $GATK_JAR -T SplitNCigarReads ' +
                         '-R {} -I {} -o {} -rf ' +
                         'ReassignOneMappingQuality -RMQF 255 -RMQT 60 ' +
                         '-U ALLOW_N_CIGAR_READS').format(
                self.ref_fa, self.dedup_bam, self.split_trim_bam)
            print(split_cmd)
            subprocess.call(split_cmd, shell=True)
        # Indel Realignment (optional): not doing indels anyway
        # 5. Base Recalibration. Docs:
        # https://software.broadinstitute.org/gatk/documentation/article?id=44
        if not os.path.exists(self.bqsr_bam):
            bqsr_cmd = ('time java -jar $GATK_JAR -T PrintReads ' +
                        '-R {} -I {} -BQSR {}_bqsr.grp -o {}').format(
                self.ref_fa, self.split_trim_bam, self.prefix, self.bqsr_bam)
            print(bqsr_cmd)
            subprocess.call(bqsr_cmd, shell=True)
        # 6. Variant calling
        if not os.path.exists(self.vcf_nofilter):
            hc_cmd = ('java -jar $GATK_JAR -T HaplotypeCaller -R {} ' +
                      '-I {} -dontUseSoftClippedBases ' +
                      '-stand_call_conf 20.0 ' +
                      '-o {}').format(
                self.ref_fa, self.bqsr_bam, self.vcf_nofilter)
            print(hc_cmd)
            subprocess.call(hc_cmd, shell=True)
        # 7. Variant filtering
        if not os.path.exists(self.vcf):
            filter_cmd = ('java -jar $GATK_JAR -T VariantFiltration ' +
                          '-R {} -V {} -window 35 -cluster 3 ' +
                          '-filterName FS -filter "FS > 30.0" ' +
                          '-filterName QD ' +
                          '-filter "QD < 2.0" -o {}').format(
                self.ref_fa, self.vcf_nofilter, self.vcf)
            print(filter_cmd)
            subprocess.call(filter_cmd, shell=True)

    # def run_bqsr(self):
    #     """Run all the BQSR commands."""
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
