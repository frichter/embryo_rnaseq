#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""CLASSES for Embryo data variant calling with GATK and samtools mpileup.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: variant calling algorithms!

"""

import re
import subprocess
import os


class bam_gatk(object):
    """Create an object and functions for GATK variant calling from a BAM."""

    def __init__(self, file_prefix, home_dir, aligner='hisat2'):
        """Create the BAM object."""
        # attach various reference file and directory locations to the obj
        self.home_dir = home_dir
        self.ref_fa = (
            '/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/' +
            'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
        ref_gatk_dir = '/sc/orga/projects/chdiTrios/Felix/dbs/gatk_resources/'
        self.ks_dbsnp = ref_gatk_dir + 'dbsnp_146.grch38.sorted.vcf'
        self.ks_indels = (
            ref_gatk_dir +
            'Mills_and_1000G_gold_standard.indels.grch38.sorted.vcf')
        self.tmp_dir = ('/sc/orga/projects/chdiTrios/Felix/' +
                        'embryo_rnaseq/tmp_dir/')
        # include aligner information in file name
        aligner_types = ['hisat2', 'star']
        self.aligner = aligner
        if self.aligner not in aligner_types:
            raise ValueError("Invalid aligner. Expected {} but got {}".format(
                aligner_types, self.aligner))
        self.prefix = file_prefix + '_' + aligner
        rm_txt = '.*ASTQ/|.*trim_q20/|_' + self.aligner
        # splitting in case file is in subdirectory w its own name
        self.id = re.sub(rm_txt, '', self.prefix).split('/')[0]
        # if file is in sub-directory then get final ID including aligner name
        # intermediate bams from the GATK steps
        self.init_file_names()

    def init_file_names(self):
        """Initialize file names for intermediate and final GATK steps."""
        self.in_sam = self.prefix + '.sam'
        if self.aligner == 'star':
            self.in_sam = self.prefix + 'Aligned.sortedByCoord.out.bam'
        elif self.aligner == 'hisat2':
            self.in_sam = self.prefix + '.sam'
        self.clean_sam = self.prefix + '_clean.bam'
        self.sorted_bam = self.prefix + '_rg_sorted.bam'
        self.dedup_bam = self.prefix + '_dedupped.bam'
        self.split_trim_bam = self.prefix + '_split.bam'
        self.bqsr_bam = self.prefix + '_bqsr.bam'
        self.vcf_nofilter = self.prefix + '_gatk3_nofilter.vcf'
        self.vcf = self.prefix + '_gatk3.vcf'
        # if file is in sub-directory but final VCF is in
        # top level directory, then update final VCF to be top-level
        if self.vcf.count(self.id) > 1:
            self.vcf = re.sub('FASTQ/' + self.id, 'FASTQ', self.vcf)

    def clean_known_sites(self):
        """Delete the automatically made index for knownsites."""
        if os.path.exists(self.ks_dbsnp + '.idx'):
            os.remove(self.ks_dbsnp + '.idx')
        if os.path.exists(self.ks_indels + '.idx'):
            os.remove(self.ks_indels + '.idx')

    def run_picard_cs(self):
        """Run Picard clean sam command (used in RNAcocktail)."""
        if not os.path.exists(self.clean_sam):
            rg_cmd = ('time java -Djava.io.tmpdir={} ' +
                      '-jar $PICARD CleanSam I={} O={}').format(
                self.tmp_dir, self.in_sam, self.clean_sam)
            print(rg_cmd)
            subprocess.call(rg_cmd, shell=True)
        else:
            print('already ran PICARD CleanSam for', self.id)

    def run_picard_rg(self):
        """Run PICARD read group commands."""
        # Add read groups, sort, mark duplicates, and create index
        if not os.path.exists(self.sorted_bam):
            rg_cmd = ('time java -Djava.io.tmpdir={} ' +
                      '-jar $PICARD AddOrReplaceReadGroups I={} ' +
                      'O={} SO=coordinate RGID=FelixRichter RGLB=Nextera ' +
                      'RGPL=ILLUMINA RGPU=machine RGSM={}').format(
                self.tmp_dir,
                # use clean_sam if using run_picard_cs, otherwise in_sam:
                self.in_sam,
                self.sorted_bam, self.id)
            print(rg_cmd)
            subprocess.call(rg_cmd, shell=True)
        else:
            print('already ran PICARD AddOrReplaceReadGroups for', self.id)

    def run_picard_md(self):
        """Run PICARD mark duplicates commands."""
        # mark MarkDuplicates docs:
        # https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php
        if not os.path.exists(self.dedup_bam):
            dup_cmd = ('time java -jar -Djava.io.tmpdir={} ' +
                       '$PICARD MarkDuplicates I={} ' +
                       'O={} CREATE_INDEX=true ' +
                       'VALIDATION_STRINGENCY=SILENT ' +
                       'M={}_gatk3_md_metrics.txt').format(
                self.tmp_dir, self.sorted_bam, self.dedup_bam, self.prefix)
            print(dup_cmd)
            subprocess.call(dup_cmd, shell=True)
        else:
            print('already ran PICARD MarkDuplicates for', self.id)

    def run_gatk_split_trim(self):
        """Run GATK3 Split'N'Trim and reassign mapping qualities."""
        if not os.path.exists(self.split_trim_bam):
            split_cmd = ('time java -Djava.io.tmpdir={} ' +
                         '-jar $GATK_JAR -T SplitNCigarReads ' +
                         '-R {} -I {} -o {} -rf ' +
                         'ReassignOneMappingQuality -RMQF 255 -RMQT 60 ' +
                         '-U ALLOW_N_CIGAR_READS').format(
                self.tmp_dir, self.ref_fa, self.dedup_bam, self.split_trim_bam)
            print(split_cmd)
            subprocess.call(split_cmd, shell=True)
        else:
            print('already ran GATK SplitNCigarReads for', self.id)

    def run_bqsr(self):
        """Run all the BQSR commands."""
        """1. Analyze patterns of covariation in the sequence dataset."""
        if not os.path.exists(self.bqsr_bam):
            # self.clean_known_sites()
            bqsr_mk_tbl_cmd = (
                'time java  -Djava.io.tmpdir={} ' +
                '-jar $GATK_JAR -T BaseRecalibrator ' +
                '-R {} -I {} -knownSites {} -knownSites {} -o {}').format(
                self.tmp_dir, self.ref_fa,
                self.split_trim_bam,
                self.ks_dbsnp, self.ks_indels,
                self.prefix + '.recal_data.table')
            print(bqsr_mk_tbl_cmd)
            subprocess.call(bqsr_mk_tbl_cmd, shell=True)
            """2. Do a second pass to analyze covariation
            remaining after recalibration
            """
            self.clean_known_sites()
            bqsr_mk_tbl_p2_cmd = (
                'time java  -Djava.io.tmpdir={} ' +
                '-jar $GATK_JAR -T BaseRecalibrator ' +
                '-R {} -I {} -knownSites {} -knownSites {} ' +
                '-BQSR {} -o {}').format(
                self.tmp_dir, self.ref_fa,
                self.split_trim_bam,
                self.ks_dbsnp, self.ks_indels,
                self.prefix + '.recal_data.table',
                self.prefix + '.recal_data_pass2.table')
            print(bqsr_mk_tbl_p2_cmd)
            subprocess.call(bqsr_mk_tbl_p2_cmd, shell=True)
            """3. Generate before/after plots (skipping).
            java -jar GenomeAnalysisTK.jar \
                -T AnalyzeCovariates \
                -R reference.fa \
                -L 20 \
                -before recal_data.table \
                -after post_recal_data.table \
                -plots recalibration_plots.pdf
            """
            """4. Apply the recalibration to your sequence data."""
            bqsr_cmd = ('time java -Djava.io.tmpdir={} ' +
                        '-jar $GATK_JAR -T PrintReads ' +
                        '-R {} -I {} -BQSR {} -o {}').format(
                self.tmp_dir, self.ref_fa, self.split_trim_bam,
                self.prefix + '.recal_data_pass2.table', self.bqsr_bam)
            print(bqsr_cmd)
            subprocess.call(bqsr_cmd, shell=True)
        else:
            print('already ran GATK BQSR for', self.id)

    def run_gatk_hc(self, gvcf=False):
        """Run GATK3 commands for variant calling."""
        if gvcf:
            gvcf_opt = '-ERC GVCF '
            self.vcf_nofilter = re.sub('.vcf$', '.g.vcf', self.vcf_nofilter)
        else:
            gvcf_opt = ''
        if not os.path.exists(self.vcf_nofilter):
            # self.clean_known_sites()
            hc_cmd = ('time java -Djava.io.tmpdir={} ' +
                      '-jar $GATK_JAR -T HaplotypeCaller -R {} ' +
                      '-I {} -dontUseSoftClippedBases ' +
                      # dbsnp uses known sites for variant annotation:
                      '--dbsnp {} ' +
                      '-stand_call_conf 20.0 ' +
                      '{}' +  # gVCF option
                      '-o {}').format(
                self.tmp_dir,
                self.ref_fa,
                self.bqsr_bam,  # switch to self.bqsr_bam
                self.ks_dbsnp,
                gvcf_opt,
                self.vcf_nofilter)
            print(hc_cmd)
            subprocess.call(hc_cmd, shell=True)
        else:
            print('already ran GATK HaplotypeCaller for', self.id)

    def run_gatk_var_filter(self):
        """Run GATK3 commands for variant filtering."""
        if not os.path.exists(self.vcf):
            filter_cmd = ('time java -Djava.io.tmpdir={} ' +
                          '-jar $GATK_JAR -T VariantFiltration ' +
                          '-R {} -V {} -window 35 -cluster 3 ' +
                          '-filterName FS -filter "FS > 30.0" ' +
                          '-filterName QD ' +
                          '-filter "QD < 2.0" -o {}').format(
                self.tmp_dir, self.ref_fa, self.vcf_nofilter, self.vcf)
            print(filter_cmd)
            subprocess.call(filter_cmd, shell=True)
        else:
            print('already ran GATK VariantFiltration for', self.id)

    def run_callable_loci_gatk(self):
        """Run GATK3's callabaleloci command."""
        self.callable_bed = self.prefix + '_callable.bed'
        self.callable_sum = self.prefix + '_callable_tbl.txt'
        if not os.path.exists(self.callable_bed):
            cl_cmd = ('time java -Djava.io.tmpdir={} ' +
                      '-jar $GATK_JAR -T CallableLoci -R {} ' +
                      '-I {} -summary {} ' +
                      '-o {}').format(
                self.tmp_dir,
                self.ref_fa,
                self.bqsr_bam,
                self.callable_sum,
                self.callable_bed)
            print(cl_cmd)
            subprocess.call(cl_cmd, shell=True)

#
