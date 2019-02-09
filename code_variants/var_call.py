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

# HISAT2 not ready for this ID:
ID="/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/95614_C1_THS_009_BxE9_2_20_17_S8_L002"

# already run (1593m22.789s + 158m41.209s = ~29.2hrs)
ID="/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/96134_C7_THS_030_BxE1_3_20_17_S31_L007"

time python var_call.py --aligner "hisat2" --bam $ID

# running on interactive2
ID="/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/96134_C7_THS_030_BxE1_3_20_17_S31_L007"
time python var_call.py --aligner "star" --bam $ID

# next:
ID="/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/92182_C5_THS_011_E13_2_23_17_S31_L007"


"""

import os
import argparse

from var_call_class import bam_gatk


def main():
    """Run GATK."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--aligner', choices=["star", "hisat2"], default=None)
    parser.add_argument('--bam', type=str, help='BAM prefix')
    args = parser.parse_args()
    home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
    os.chdir(home_dir)
    bam_i = bam_gatk(args.bam, home_dir, aligner=args.aligner)
    if not os.path.exists(bam_i.in_sam):
        return 'BAM/SAM not available ' + bam_i.in_sam
    bam_i.run_picard_rg()  # 12m
    bam_i.run_picard_md()  # also 12m
    bam_i.run_gatk_split_trim()
    # Indel Realignment (optional): not doing indels currently
    bam_i.run_bqsr()
    bam_i.run_gatk_hc()
    bam_i.run_gatk_var_filter()
    return 'GATK done for ' + args.bam


if __name__ == "__main__":
    done_msg = main()
    print(done_msg)


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
