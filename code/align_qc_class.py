"""Classes to align and QC check for Embryo data.

Felix Richter
felix.richter@icahn.mssm.edu
12/26/2018
Description: run alignment commands, FastQC scripts and anything else fastq
             related (e.g., maybe trimming in the future)


cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ

# trim_galore uses python/2.7.14 and py_packages/2.7 (loaded automatically)
module purge
module load trim_galore/0.4.5 fastqc/0.11.7 trimmomatic/0.36

module purge
module load python/3.5.0 py_packages/3.5
module load fastqc/0.11.7
module load trimmomatic/0.36
module load hisat2/2.0.5 star/2.6.1d

source ~/venv_ore/bin/activate

# decompose file name into:
96295_C1_THS_029_E9_3_15_17_S38_L008_R2_001

# consider multiQC
# Were adatpers trimmed? No!
python
"""

import re
import subprocess
import os


class fq_pair(object):
    """Create a FastQ pair object."""

    def __init__(self, pair_file_loc, home_dir, hisat2_idx):
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

    def RunHISAT2(self):
        """Run HISAT2.

        Manual: https://ccb.jhu.edu/software/hisat2/manual.shtml

        Notes on HISAT2:
            -p is number of cores to use
        """
        # confirm sam file isn't already made, then run hisat2
        if not os.path.exists(self.prefix + '.sam'):
            hisat2_cmd = ('time hisat2 --time -x %s -1 %s -2 %s -S %s.sam ' +
                          '--un-conc %s_noPEalign -p 24') % \
                (self.hisat2_idx, self.r1, self.r2, self.prefix, self.prefix)
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
        """

    def RunSTAR(self):
        """Run STAR aligner.

        Manual:
        https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
        """
        pass


class fq_pair_qc(fq_pair):
    """Create a FastQ object with QC commands."""

    def TrimAdapters(self):
        """Trime adapters.

        Manual/methods:
        https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
        """
        trimmed_r1 = re.sub('.fastq.gz', '_trimmed.fq.gz', self.r1)
        trimmed_r2 = re.sub('.fastq.gz', '_trimmed.fq.gz', self.r2)
        if not os.path.exists(trimmed_r1):
            trim_cmd = 'time trim_galore --gzip --paired {} {}'.format(
                self.r1, self.r2)
            # print(trim_cmd)
            subprocess.call(trim_cmd, shell=True)
            # any way to check if successful?
        self.r1 = trimmed_r1
        self.r2 = trimmed_r2

    def FastQC(self):
        """Run FastQC command.

        Manual:
        https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
        """
        fastqc_cmd = 'time fastqc %s'
        print('running fastqc for %s' % self.r1)
        subprocess.call(fastqc_cmd % self.r1, shell=True)
        print('running fastqc for %s' % self.r2)
        subprocess.call(fastqc_cmd % self.r2, shell=True)


"""

## what is the --dta tag? did I use this? Need this for stringtie

STAR info:
## not sure which of these:
module load star/2.5.3a
module load rna-star/2.4.0d

# convert SAM to BAM
time samtools view -u D1_D2_odd_repeat/D1_CTRL1.sam |
samtools sort -T D1_D2_odd_repeat/D1_CTRL1.temp -o \
D1_D2_odd_repeat/D1_CTRL1.sorted.bam -@ 12


"""

#
