#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""CLASSES for Embryo data callable loci.

Felix Richter
felix.richter@icahn.mssm.edu
2/17/2019
Description: Class/object with functions and attributes relevant to
    determining which regions are callable

"""


import subprocess
import os


import pybedtools
from pybedtools import BedTool


class call_loci(object):
    """Create an object and functions for callable loci."""

    def __init__(self, id, home_dir):
        """Create the callable locus object."""
        self.aligner_ls = ['star', 'hisat2']
        self.id = id
        self.subdir = '{}callable_comparison/{}/'.format(home_dir, id)
        if not os.path.exists(self.subdir):
            os.mkdir(self.subdir)
        # original input file (output from GATK callableloci)
        self.call_loci_ls = [
            '{}FASTQ/{}/{}_{}_callable.bed'.format(home_dir, id, id, i)
            for i in self.aligner_ls]
        # intermediate files:
        self.call_only_ls = [
            '{}callable_{}.bed'.format(self.subdir, i)
            for i in self.aligner_ls]
        self.call_inter = '{}callable_{}.bed'.format(self.subdir, 'inter')
        self.call_union = '{}callable_{}.bed'.format(self.subdir, 'union')
        self.callable_ls = self.call_only_ls + [
            self.call_inter, self.call_union]
        callable_fs = ['star', 'hisat2', 'intersect', 'union']
        self.callable_dict = dict(zip(callable_fs, self.callable_ls))
        self.len_dict = dict(zip(self.callable_ls, [0]*len(self.callable_ls)))
        pybedtools.set_tempdir(home_dir + '/tmp_dir/')

    def subset_callable_per_f(self, call_loci, call_only):
        """Create the callable region subset bed file for each aligner."""
        if os.path.exists(call_only):
            return 'already created ' + call_only
        with open(call_loci, 'r') as in_f, open(call_only, 'w') as out_f:
            for line in in_f:
                if 'CALLABLE' in line:
                    out_f.write(line)
        return call_only + ' done'

    def subset_callable_loop(self):
        """Loop over the callable loci and subset."""
        for i, j in zip(self.call_loci_ls, self.call_only_ls):
            # print(i, j)
            status = self.subset_callable_per_f(i, j)
            print(status)

    def intersect_callable(self):
        """Calculate intersection of callable loci.

        Source:
        https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

        """
        if os.path.exists(self.call_inter):
            return 'intersection already done for ' + self.call_inter
        bed_a = BedTool(self.call_only_ls[0])
        # turn into for loop if using more aligners
        bed_b = BedTool(self.call_only_ls[1])
        bed_inter = bed_a.intersect(bed_b)
        # see if it needs to be merged. Don't think so but double check
        # bed_inter_merged = bed_inter.merge()
        bed_inter.saveas(self.call_inter)

    def union_callable(self):
        """Get the union of callable loci."""
        if os.path.exists(self.call_union):
            return 'Union already done for ' + self.call_union
        cat_sort_cmd = 'cat {} | sort -V -k1,1 -k2,2n > {}'.format(
            ' '.join([i for i in self.call_only_ls]), self.call_union)
        print(cat_sort_cmd)
        subprocess.call(cat_sort_cmd, shell=True)
        bed = BedTool(self.call_union)
        bed_merged = bed.merge()
        bed_merged.saveas(self.call_union)

    def intersect_w_known_loci(self, known_f, known_folder):
        """Intersect with callable regions."""
        # print(known_f, known_folder)
        if not os.path.exists(self.subdir + known_folder):
            os.mkdir(self.subdir + known_folder)
        bed_known = BedTool(known_f)
        for k, v in self.callable_dict.items():
            out_f = '{}{}/{}.bed'.format(self.subdir, known_folder, k)
            self.len_dict[out_f] = 0
            if os.path.exists(out_f):
                continue
            bed_i = BedTool(v)
            bed_inter = bed_known.intersect(bed_i)
            bed_inter.saveas(out_f)


""" """
#
