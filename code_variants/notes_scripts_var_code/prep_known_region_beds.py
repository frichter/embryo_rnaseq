#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Prepare known regions for benchmarking with callable loci.

Felix Richter
felix.richter@icahn.mssm.edu
2/18/2019
Description: Prepare known regions


cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/known_region
module purge
module load bedtools/2.27.0
module load python/3.5.0
module load py_packages/3.5
module load ucsc-utils/2015-04-07
python

"""


import re
import subprocess

from pybedtools import BedTool


"""Download from local on data4
ssh data4
# preparing CCDS exonic regions
time wget -O clinvar_locs_hg38.txt \
https://www.dropbox.com/s/ocfvu82s4ondlcv/clinvar_locs_hg38.txt?dl=1
time wget -O exon_locs_hg38.txt \
https://www.dropbox.com/s/65fczfzhqgl7alr/exon_locs_hg38.txt?dl=1
"""


"""Exonic CCDS regions from UCSC"""

in_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/known_region/'
in_loc = in_dir + 'exon_locs_hg38.txt'
out_loc = in_dir + 'exon_locs_hg38.bed'
with open(in_loc, 'r') as in_f, open(out_loc, 'w') as out_f:
    header = next(in_f)
    for line in in_f:
        line_list = line.strip().split('\t')
        # print(line_list)
        # need to remove chr and convert chrM to MT
        chrom = re.sub('chrM', 'MT', line_list[0])
        chrom = re.sub('chr', '', chrom)
        start_list = line_list[1].split(',')[:-1]
        end_list = line_list[2].split(',')[:-1]
        for start_i, end_i in zip(start_list, end_list):
            out_line = '\t'.join([chrom, start_i, end_i]) + '\n'
            _ = out_f.write(out_line)
        # break


"""preparing Clinvar pathogenic regions/variants."""

in_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/known_region/'
in_loc = in_dir + 'clinvar_locs_hg38.txt'
out_loc = in_dir + 'clinvar_hg38_P_LP.bed'
with open(in_loc, 'r') as in_f, open(out_loc, 'w') as out_f:
    header = next(in_f)
    for line in in_f:
        line_list = line.strip().split('\t')
        # print(line_list)
        # need to remove chr and convert chrM to MT
        chrom = re.sub('chrM', 'MT', line_list[0])
        chrom = re.sub('chr', '', chrom)
        # if line_list[3] == 'Pathogenic':
        if 'athogenic' in line_list[3]:
            out_line = '\t'.join([chrom, line_list[1], line_list[2]]) + '\n'
            # print(out_line)
            _ = out_f.write(out_line)


"""Sort and merge output bed files."""

out_loc_list = ['exon_locs_hg38.bed', 'clinvar_hg38_P_LP.bed',
                'clinvar_hg38_P.bed']
in_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/known_region/'
for bed_loc in out_loc_list:
    print(bed_loc)
    bed_loc = in_dir + bed_loc
    bed_sorted_loc = re.sub('.bed', '_sorted.bed', bed_loc)
    sort_cmd = 'sort -V -k1,1 -k2,2n {} > {}'.format(bed_loc, bed_sorted_loc)
    subprocess.call(sort_cmd, shell=True)
    bed_a = BedTool(bed_sorted_loc)
    bed_merge = bed_a.merge()
    bed_merge.saveas(bed_sorted_loc)


"""Downloading known coverage tracks. THESE ARE NOT ACTUALLY AVAILABLE...

UCSC > hg19 (not available for hg38) > Variation > Genome Coverage %
or > Exome Coverage %
Only using 30x. Note that Genome Coverage % is a dead link
# possibly submit bug here: https://genome.ucsc.edu/contacts.html
# See:
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=gnomadGenomesAvgCoverage&hgta_table=gnomadGenomesMedianCoverage&hgta_doSchema=describe+table+schema

# 1) liftover from hg19 to hg38
in_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/known_region/'
os.chdir(in_dir)
in_f = ['exac_30x.bed', 'gnomad_30x.bed'][0]
out_f = re.sub('.bed', '_hg38.bed', in_f)
out_f_unmapped = re.sub('.bed', '_unmapped.txt', out_f)
chain_f='/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/hg19ToHg38.over.chain'
lift_cmd = 'liftOver -multiple -minMatch=0.1 {} {} {} {}'.format(
    in_f, chain_f, out_f, out_f_unmapped)
subprocess.call(lift_cmd, shell=True)
2) remove chr prefix, convert chrM to MT

"""

""" """
