"""Notes."""


import glob
import os
import re


def delete_sam_f(in_sam, hs_iter):
    """Delete the pertinent SAM files."""
    in_sam_f_iter = glob.iglob(re.sub('.sam', '*', in_sam))
    for in_sam_f in in_sam_f_iter:
        just_rerun = any([i in in_sam_f for i in hs_iter])
        if not just_rerun:
            if os.path.exists(in_sam_f):
                print('deleting', in_sam_f)
                # os.remove(in_sam_f)
        else:
            print('KEEPING', in_sam_f)
            # raise ValueError()


# hs_iter = [
#     '90124_C1_THS_015_BxE1_3_1_17_S12_L003',
#     '90124_C1_THS_015_BxE5_3_1_17_S8_L002',
#     '90099_C4_THS_010_BxE9_3_1_17_S36_L008',
#     '90099_C4_THS_010_BxE3_3_1_17_S28_L006',
#     '91122_C1_THS_027_BxE2_3_17_17_S20_L004']
# if '92182_C5_THS_011_E13_2_23_17_S31_L007_hisat2' in bam_i.vcf:
#     print('SKIPPING 92182_C5_THS_011_E13_2_23_17_S31_L007_hisat2')
#     continue
# delete the relevant input SAM files as well!
# delete_sam_f(bam_i.in_sam, hs_iter)

"""
Failed files
Exception in thread "main" htsjdk.samtools.SAMFormatException:
Error parsing text SAM file. Not enough fields

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/
tail 90099_C4_THS_010_BxE10_3_1_17_S2_L001_hisat2.sam
Line 96058617
SND00311:234:CC56TANXX:1:2305:13307:65206

rm 90099_C4_THS_010_BxE10_3_1_17_S2_L001_hisat2_rg_sorted.bam
# STAR number of input reads 45255155

# for HISAT2, check if the GATK issues are resolved after re-running
# HISAT2. e.g:
cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants/logs/gatk_old
grep '90124_C1_THS_015_BxE5_3_1_17_S8_L002_hisat2' *std*
ls *_13.std*

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants
module purge
module load gatk/3.6-0
module load picard/2.7.1
module load python/3.5.0 py_packages/3.5

ID="/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/90124_C1_THS_015_BxE5_3_1_17_S8_L002"
time python var_call.py --aligner "hisat2" --bam $ID
"""
