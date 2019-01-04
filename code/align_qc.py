"""Align and QC check for Embryo data.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: run alignment commands, FastQC scripts and anything else fastq
             related (e.g., maybe trimming in the future)

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code

 # trim_galore uses python/2.7.14 and py_packages/2.7 (loaded automatically)
 module purge
 module load trim_galore/0.4.5 fastqc/0.11.7 trimmomatic/0.36

 module purge
 module load python/3.5.0 py_packages/3.5
 module load hisat2/2.0.5 star/2.6.1d

 source ~/venv_ore/bin/activate

 # decompose file name into:
 96295_C1_THS_029_E9_3_15_17_S38_L008_R2_001

python
"""

import os
# import glob

from align_qc_class import fq_pair_qc

# fq_file_loc = 'FASTQ/75888_C4_THS_014_BxE8_2_28_17_S18_L004'
fq_file_loc = 'FASTQ/76448_C7_THS_025_BxE3_3_13_17_S30_L006'

home_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/'
os.chdir(home_dir)

fq_i = fq_pair_qc(
    pair_file_loc=fq_file_loc,
    home_dir=home_dir,
    hisat2_idx='grch38_snp_tran/genome_snp_tran')
fq_i.TrimAdapters()
fq_i.FastQC()
fq_i.RunHISAT2()
fq_i.RunSTAR()

""" Looping over multiple FastQ files (use this or BSUB)
fq_file_iter = glob.iglob(home_dir + '*.fastq.gz')
fq_file_list = [re.sub('_R[12].*fastq.gz', '', i) for i in fq_file_iter]
keep uniques
len(fq_file_list)
fq_file_list = list(set(fq_file_list))
len(fq_file_list)
for fq_i_file_loc in fq_file_list:
    fq_i = fq_pair_qc(
        pair_file_loc=fq_i_file_loc,
        home_dir=home_dir,
        hisat2_idx='grch38_snp_tran/genome_snp_tran')
    fq_i.TrimAdapters()
    fq_i.FastQC()
    fq_i.RunHISAT2()
    fq_i.RunSTAR()

"""

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
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz
"""
