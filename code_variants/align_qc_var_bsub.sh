#BSUB -W 12:00
#BSUB -q premium
#BUSB -n 30
#BSUB -R "rusage[mem=45000]"
#BSUB -P acc_schade01a
#BSUB -J "embryo_align_star[2-81]"
#BSUB -m manda
#BSUB -o logs/%J_%I.stdout
#BSUB -e logs/%J_%I.stderr


cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants

# submit with 
# bsub < align_qc_var_bsub.sh
# max indices [1-81]

let i=$LSB_JOBINDEX
ID=$(tail -n+$i ../metadata/fq_prefix_list.txt | head -n1)
echo $ID

#################################### Trimming
## Settings: -W 4:00, -R mem=1000, -J embryo_trim)
# module purge
# module load trim_galore/0.4.5
# python align_qc.py --trimonly --fq $ID


#################################### HISAT2 alignment
## Settings: (-W 144:00, -R mem=8500, -n 30 -J embryo_align_hisat2
# module purge
# module load fastqc/0.11.7
# module load hisat2/2.0.5
# module load python/3.5.0 py_packages/3.5
# python align_qc.py --hisat2 --fq $ID
## HISAT2 is advertised as requiring only 3-4 Gb

#################################### STAR alignment
## Settings: -W 4:00, -R mem=38000, -n 30 -J embryo_align_star
module purge
module load python/3.5.0
module load star/2.6.1d
python align_qc.py --star --fq $ID

#################################### Variant calling
## Settings: -W 6:00, -R mem=35000, -J embryo_var_gatk
# module purge
# module load gatk/3.6-0
# module load picard/2.7.1
# module load python/3.5.0 py_packages/3.5
# time python var_call.py --aligner "hisat2" --bam $ID
# python var_call.py --aligner "star" --bam $ID

############ Minerva resources
# decide on partition w https://hpc.mssm.edu/resources/hardware