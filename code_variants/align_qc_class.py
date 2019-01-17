#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Classes to align and QC check for Embryo data.

Felix Richter
felix.richter@icahn.mssm.edu
12/26/2018
Description: run alignment commands, FastQC scripts and anything else fastq
             related (e.g., maybe trimming in the future)

"""

import re
import subprocess
import os


class fq_pair(object):
    """Create a FastQ pair object."""

    def __init__(self, pair_file_loc, home_dir, hisat2_idx, star_idx, n):
        """Create an object for the FASTQ file."""
        # assign the filenames to the instance
        self.pair_file_loc = pair_file_loc
        self.r1 = pair_file_loc + '_R1_001.fastq.gz'
        self.r2 = pair_file_loc + '_R2_001.fastq.gz'
        # obtain the Lane from the filename
        self.lane = re.search('L[0-9]{3}', pair_file_loc)
        # obtain blinded ID
        self.id = re.search('^[0-9]{6}', pair_file_loc)
        # create the output sam file
        self.prefix = pair_file_loc
        # go to home directory
        self.home_dir = home_dir
        # hisat2 index
        self.hisat2_idx = hisat2_idx
        # STAR aligner index
        self.star_idx = star_idx
        self.n_cores = n

    def RunHISAT2(self):
        """Run HISAT2.

        Manual: https://ccb.jhu.edu/software/hisat2/manual.shtml

        Notes on HISAT2:
            -p is number of cores to use
        """
        # confirm sam file isn't already made, then run hisat2
        if not os.path.exists(self.prefix + '_hisat2.sam'):
            hisat2_cmd = ('time hisat2 --time -x {} -1 {} -2 {} -S {}.sam ' +
                          '--met-file {}_metrics.txt ' +
                          '--dta --un-conc {}_noPEalign -p {}').format(
                self.hisat2_idx, self.r1, self.r2, self.prefix + '_hisat2',
                self.prefix + '_hisat2', self.prefix + '_hisat2', self.n_cores)
            print(hisat2_cmd)
            subprocess.call(hisat2_cmd, shell=True)
        else:
            print(self.prefix + '.sam already made')
        """
        other HISAT2 options to consider implementing
        trim 5' or 3' ends (--trim3 <int> and --trim5 <int>)
        --phred33 or --phred64
        --known-splicesite-infile, --novel-splicesite-outfile,
             --novel-splicesite-infile
        reads that fail to align: --un <path> --un-conc <path>
        alignment metrics: --met-file <path> --met-stderr
        readgroup IDs with --rg-id <text>

        note that hisat2-build can be used for custom transcriptome
        indices (if you want to use those instead of a genome index)
        good examples:
        https://bioinformatics.sciberg.com/2018/11/01/mapping-with-hisat2/
        """

    def RunSTAR(self):
        """Run STAR aligner.

        Manual:
        https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
        """
        out_bam = self.prefix + '_starAligned.sortedByCoord.out.bam'
        if not os.path.exists(out_bam):
            star_cmd = ('time STAR --runThreadN {} --genomeDir {} ' +
                        '--readFilesCommand zcat --outFileNamePrefix {} ' +
                        '--outSAMtype BAM SortedByCoordinate ' +
                        '--twopassMode Basic ' +
                        '--readFilesIn {} {}').format(
                self.n_cores, self.star_idx, self.prefix + '_star',
                self.r1, self.r2)
            print(star_cmd)
            subprocess.call(star_cmd, shell=True)
        else:
            print(self.prefix + ' STAR alignment already done')


class fq_pair_qc(fq_pair):
    """Create a FastQ object with QC commands."""

    def TrimAdapters(self):
        """Trime adapters.

        Manual/methods:
        https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

        Use the val_1 and val_2 files
        (Source: https://www.biostars.org/p/256388/#256612)

        Checking ASCII scores:
        https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
        """
        trimmed_r1 = re.sub('.fastq.gz', '_val_1.fq.gz', self.r1)
        trimmed_r2 = re.sub('.fastq.gz', '_val_2.fq.gz', self.r2)
        # remove files if not completed previously
        # if not self.check_trim_complete(self.r1):
        #     os.remove(trimmed_r1)
        # if not self.check_trim_complete(self.r2):
        #     os.remove(trimmed_r2)
        if not os.path.exists(trimmed_r1):
            trim_cmd = ('time trim_galore -o {} --gzip ' +
                        '--quality 0 --paired {} {}').format(
                self.home_dir + 'FASTQ/', self.r1, self.r2)
            print(trim_cmd)
            subprocess.call(trim_cmd, shell=True)
            # any way to check if successful? If report is >15 lines
        else:
            print(self.prefix + ' already trimmed')
        self.r1 = trimmed_r1
        self.r2 = trimmed_r2

    def FastQC(self):
        """Run FastQC command.

        Manual:
        https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
        """
        out_f = re.sub('.fq.gz|.fastq.gz', '_fastqc.html', self.r1)
        if not os.path.exists(out_f):
            fastqc_cmd = 'time fastqc %s'
            print('running fastqc for %s' % self.r1)
            subprocess.call(fastqc_cmd % self.r1, shell=True)
            print('running fastqc for %s' % self.r2)
            subprocess.call(fastqc_cmd % self.r2, shell=True)
        else:
            print('fastqc plots already made for', self.r1)


"""

# convert SAM to BAM
time samtools view -u D1_D2_odd_repeat/D1_CTRL1.sam |
samtools sort -T D1_D2_odd_repeat/D1_CTRL1.temp -o \
D1_D2_odd_repeat/D1_CTRL1.sorted.bam -@ 12


"""

#
