#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Notes for variant calling from RNAseq in embryo data."""

"""GATK quick start guide
https://software.broadinstitute.org/gatk/documentation/quickstart?v=3
java -version
java -jar $GATK_JAR -h
java -jar $PICARD -h
"""


"""RNA cocktail options used (all GATK options have been x-checked)

HaplotypeCaller: -stand call conf 20.0
-stand emit conf 20.0 -A StrandBiasBySample -A StrandAlleleCountsBySampleGATK

VariantFiltration: -window 35 -cluster 3 -filterName FS -filter
"FS >30.0" -filterName QD -filter "QD <2.0"


v3.5-0-g36282e4
(picard 1.129)

samtools mpileup -C50 -d 100000
bcftools filter -s LowQual -e ‘%QUAL<20 —— DP>10000’


"""
