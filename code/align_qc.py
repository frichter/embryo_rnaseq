#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Align and QC check for Embryo data.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: run alignment commands, FastQC scripts and anything else fastq
             related (e.g., trimming)

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code

# trim_galore uses python/2.7.14 and py_packages/2.7 (loaded automatically)
module purge
module load trim_galore/0.4.5 fastqc/0.11.7 trimmomatic/0.36
module load hisat2/2.0.5 star/2.6.1d

python
"""

import os
import glob
import re

from align_qc_class import fq_pair_qc

# fq_file_loc = 'FASTQ/75888_C4_THS_014_BxE8_2_28_17_S18_L004'
fq_file_loc = 'FASTQ/76448_C7_THS_025_BxE3_3_13_17_S30_L006'
db_loc = '/sc/orga/projects/chdiTrios/Felix/dbs/'

fq_i = fq_pair_qc(
    pair_file_loc=fq_file_loc,
    home_dir='/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/',
    hisat2_idx=db_loc + 'grch38_snp_tran/genome_snp_tran')
os.chdir(fq_i.home_dir)
fq_i.TrimAdapters()
print(fq_i.r1, fq_i.r2)
fq_i.FastQC()
fq_i.RunHISAT2()

star_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/grch38_star/'
star_cmd = ('time STAR --runThreadN 24 --genomeDir {} ' +
            '--readFilesCommand zcat --outFileNamePrefix {} ' +
            '--outSAMtype BAM SortedByCoordinate ' +
            '--twopassMode Basic ' +
            '--readFilesIn {} {}').format(
    star_dir, fq_i.prefix + '_star', fq_i.r1, fq_i.r2)
print(star_cmd)
# fq_i.RunSTAR()

""" Looping over multiple FastQ files (use this or BSUB)."""

home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
os.chdir(home_dir)
fq_file_iter = glob.iglob(home_dir + '*.fastq.gz')
fq_file_list = [re.sub('_R[12].*fastq.gz', '', i) for i in fq_file_iter]
# keep uniques
print(len(fq_file_list))
fq_file_list = list(set(fq_file_list))
print(len(fq_file_list))
for fq_i_file_loc in fq_file_list:
    fq_i = fq_pair_qc(
        pair_file_loc=fq_i_file_loc,
        home_dir=home_dir,
        hisat2_idx='grch38_snp_tran/genome_snp_tran')
    fq_i.TrimAdapters()
    fq_i.FastQC()
    fq_i.RunHISAT2()
    fq_i.RunSTAR()


"""Notes for preparing

# after creating github locally/online, syncing with exising
# directory on minerva:
https://stackoverflow.com/a/29078055/10688049

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ

# testing trim_galore
# FastQC and trim_galore
time trim_galore --gzip --paired \
75888_C4_THS_014_BxE8_2_28_17_S18_L004_R1_001.fastq.gz \
75888_C4_THS_014_BxE8_2_28_17_S18_L004_R2_001.fastq.gz

# download HISAT2 index
cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/

cd /sc/orga/projects/chdiTrios/Felix/dbs/
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz
time tar -xvzf grch38_snp_tran.tar.gz

module purge
module load python/3.5.0 py_packages/3.5

# generating your own STAR index


cd /sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/\
Homo_sapiens.GRCh38.94.gtf.gz

time gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
time gunzip Homo_sapiens.GRCh38.94.gtf.gz


module purge
module load star/2.6.1d

ENSEMBL_DIR="/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/"
ENSEMBL_FA=$ENSEMBL_DIR"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ENSEMBL_GTF=$ENSEMBL_DIR"Homo_sapiens.GRCh38.94.gtf"

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/

# for sjdbOverhang ReadLength-1: got ReadLength=126 from fastq files outputs

# testing STAR
time STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir grch38_star/ \
--genomeFastaFiles $ENSEMBL_FA \
--sjdbGTFfile $ENSEMBL_GTF \
--sjdbOverhang 125

"""
