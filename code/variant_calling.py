#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Embryo data variant calling with GATK and samtools mpileup.

Felix Richter
felix.richter@icahn.mssm.edu
1/4/2019
Description: variant calling algorithms!

"""


"""RNA cocktail options used

samtools mpileup -C50 -d 100000bcftools filter
-s LowQual -e ‘%QUAL<20 —— DP>10000’


v3.5-0-g36282e4
(picard 1.129)

Picard AddOrReplaceReadGroups: SO=coordinate

Picard MarkDuplicates: CREATE INDEX=true VALIDATION STRINGENCY=SILENTGATK

SplitNCigarReads: -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
-U ALLOW N CIGAR READSGATK

HaplotypeCaller: -stand call conf 20.0
-stand emit conf 20.0 -A StrandBiasBySample -A StrandAlleleCountsBySampleGATK

VariantFiltration: -window 35 -cluster 3 -filterName FS -filter
"FS >30.0" -filterName QD -filter "QD <2.0"

"""
