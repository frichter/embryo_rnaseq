#BSUB -W 144:00
#BSUB -q premium
#BUSB -n 30
#BSUB -R "rusage[mem=34000]"
#BSUB -P acc_schade01a
#BSUB -J "embryo_var_star_gatk[70]"
#BSUB -m manda
#BSUB -o logs/%J_%I.stdout
#BSUB -e logs/%J_%I.stderr


cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants

# submit with 
# bsub < var_call_qc_bsub.sh
# max indices [1-81]

let i=$LSB_JOBINDEX
ID=$(tail -n+$i ../metadata/fq_subdir_prefix_list.txt | head -n1)
echo $ID

#################################### Variant calling
## Settings: -W 144:00, -R mem=34000
## -J embryo_var_hs_gatk or embryo_var_star_gatk
module purge
module load gatk/3.6-0
module load picard/2.7.1
module load python/3.5.0 py_packages/3.5
# time python var_call.py --aligner "hisat2" --bam $ID
time python var_call.py --aligner "star" --bam $ID
# Max time because I don't want to deal with figuring out what
# happens if it crashes

#################################### Callable regions
## probably just the same GATK settings -W 10:00
## -J embryo_callable_star embryo_callable_hs
# module load gatk/3.6-0
# module load picard/2.7.1
# module load python/3.5.0 py_packages/3.5
# 
# time python callable_regions.py --aligner "star" --bam $ID
# time python callable_regions.py --aligner "hisat2" --bam $ID


############ Minerva resources
# decide on partition w https://hpc.mssm.edu/resources/hardware
# 256 Gb/node on Manda (mothra not listed)