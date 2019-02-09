#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Notes for prepping known sites for BQSR."""


"""Docs for known sites:
https://software.broadinstitute.org/gatk/documentation/article.php?id=1247

The most recent dbSNP release (build ID > 132)
Mills_and_1000G_gold_standard.indels.b37.vcf
cd /sc/orga/projects/chdiTrios/Felix/dbs/gatk_resources

# two sources (FTP and google cloud)
ftp://ftp.broadinstitute.org/bundle/hg38/
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0

# indels
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# dbsnp
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi

# need to use GATK bundle reference genome I guess
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
time wget \
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz

# according to this, the google cloude files are the most uptodate:
https://gatkforums.broadinstitute.org/gatk/discussion/23170/gatk-resource-bundle

time gunzip Homo_sapiens_assembly38.fasta.gz

# Docs for running BQSR
https://software.broadinstitute.org/gatk/documentation/article?id=2801

convert from hg38 to GRCh38

# solution to this problem:
https://gatkforums.broadinstitute.org/gatk/discussion/8895/picard-sort-vcf-error

cd /sc/orga/projects/chdiTrios/Felix/dbs/gatk_resources
module load gatk/3.6-0
module load picard/2.7.1
module load samtools/1.9
module load python/3.5.0 py_packages/3.5
python

import gzip
import re

ref_gatk_dir = '/sc/orga/projects/chdiTrios/Felix/dbs/gatk_resources/'
ref_dict = ('/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/' +
            'Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai')

ref_list = []
with open(ref_dict, 'r') as f:
    for line in f:
        ref_list.append(line.split('\t')[0])


f = ref_gatk_dir + 'dbsnp_146.hg38.vcf.gz'
# f = ref_gatk_dir + 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
f_out_loc = re.sub('hg38|assembly38', 'grch38', f)
f_out_loc = re.sub('.gz', '', f_out_loc)
print(f_out_loc)
count = 0
contig_header_dict = {}
with gzip.open(f, 'r') as f_in, open(f_out_loc, 'w') as f_out:
    for line in f_in:
        line = str(line, 'utf-8')
        # print(line)
        # convert from hg38 to GRCh38
        line = re.sub('contig=<ID=chrM', 'contig=<ID=MT', line)
        line = re.sub('contig=<ID=chr', 'contig=<ID=', line)
        # exclude line if it has assembly information
        if 'contig=<ID=' in line:
            contig = line.split(',')[0].split('=')[-1]
            # print(contig)
            contig_header_dict[contig] = line
            continue
        line = re.sub('^chrM', 'MT', line)
        line = re.sub('^chr', '', line)
        line = re.sub('assembly=20', 'assembly=null', line)
        # only keep line if it has a contig in the reference dictionary
        if not line.startswith('#'):
            contig = line.split('\t')[0]
            if contig in ref_list:
                _ = f_out.write(line)
        else:
            _ = f_out.write(line)
        count += 1
        if count % 10000 == 0:
            # print(count/1275224) # for Mills_and_1000G_gold_standard
            # 149125263 for dbsnp 146
            print(count)

print(count)


## Trying just converting the knownSites files to GRCH38 and running the
SortVcf command

ENSEMBL_DIR="/sc/orga/projects/chdiTrios/Felix/dbs/grch38_ens94/"
ENSEMBL_FA=$ENSEMBL_DIR"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GATK_DICT=$ENSEMBL_DIR"Homo_sapiens.GRCh38.dna.primary_assembly.dict"

cd /sc/orga/projects/chdiTrios/Felix/dbs/gatk_resources

time java -jar $PICARD SortVcf \
    I=Mills_and_1000G_gold_standard.indels.grch38.vcf \
    O=Mills_and_1000G_gold_standard.indels.grch38.sorted.vcf \
    SEQUENCE_DICTIONARY=$GATK_DICT

time java -jar $PICARD SortVcf \
    I=dbsnp_146.grch38.vcf \
    O=dbsnp_146.grch38.sorted.vcf \
    SEQUENCE_DICTIONARY=$GATK_DICT

real    87m30.823s
user    212m42.771s

BQSR gives this error:
##### ERROR MESSAGE: Input files knownSites and reference
have incompatible contigs.

# actually was a problem with the idx file so I deleted that

# removing the contig=<ID= lines
f = 'Mills_and_1000G_gold_standard.indels.grch38.sorted.vcf'
f_out_loc = re.sub('sorted', 'sorted_noassemb', f)
with open(f, 'r') as f_in, open(f_out_loc, 'w') as f_out:
    for line in f_in:
        line = str(line, 'utf-8')
        # exclude line if it has assembly information
        if 'contig=<ID=' in line:
            continue
        _ = f_out.write(line)

Soooo many bugs in GATK. That I think could be fixed very easily
and make the software soooo much more robust

1. no matching contigs. Even though GATK misleads by saying hg38/GRCh38
are the same, they are not (chr vs not, MT vs M, random bits not named
the same). Fine I'll focus on the main chr by renaming those and
throw out all my other data

2. Even though I ran GATK PICARD sort as the first step of
my analysis, somehow my BAMs are not sorted appropriately! Luckily
the error page is helpful
https://gatkforums.broadinstitute.org/gatk/discussion/1328/errors-about-contigs-in-bam-or-vcf-files-not-being-properly-ordered-or-sorted

3. Even though I renamed and re-ordered, I still got an error
bc an alternate contig was not an exact match, even though
the reorder page specifically says
"Be aware that this tool will drop reads that don't have
equivalent contigs in the new reference"
So I scroll down to the comments which say I should use
ALLOW_INCOMPLETE_DICT_CONCORDANCE=TRUE

This worries me because I don't know what other side effects this
setting might have... but nevertheless, onward

4. New error
Exception in thread "main" htsjdk.samtools.SAMFormatException: SAM
validation error: ERROR: Read name SND00311:233:CC5P3ANXX:4:1301:12998:5699,
No real operator (M|I|D|N) in CIGAR
wtf? I created these BAMs/SAMs with PICARD effectively from scratch
and I still get this error? Looks like I'm not the first one...
https://gatkforums.broadinstitute.org/gatk/discussion/6571/picard-reordersam-error-read-cigar-m-operator-maps-off-end-of-reference
"""
