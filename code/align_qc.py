#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Align and QC check for Embryo data.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: run alignment commands, FastQC scripts and anything else fastq
             related (e.g., trimming)
"""

import os
# import glob
# import re
import argparse

from align_qc_class import fq_pair_qc


def main():
    """Run main function to of analysis."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--trimonly', default=False, action='store_true')
    parser.add_argument('--fq', type=str, help='FQ prefix')
    args = parser.parse_args()
    db_loc = '/sc/orga/projects/chdiTrios/Felix/dbs/'
    fq_i = fq_pair_qc(
        pair_file_loc=args.fq,
        home_dir='/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/',
        hisat2_idx=db_loc + 'grch38_snp_tran/genome_snp_tran')
    os.chdir(fq_i.home_dir)
    fq_i.TrimAdapters()  # takes about 200 minutes
    if args.trimonly:
        return args.fq + ' done (trimming only!)'
    print(fq_i.r1, fq_i.r2)
    fq_i.RunHISAT2()
    # fq_i.RunSTAR()
    fq_i.FastQC()
    return args.fq + ' done'


if __name__ == "__main__":
    done_msg = main()
    print(done_msg)
