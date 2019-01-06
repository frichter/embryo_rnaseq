#BSUB -W 12:00
#BSUB -q premium
#BUSB -n 20
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_schade01a
#BSUB -J "embryo_align_hisat2[11-81]"
#BSUB -m mothra
#BSUB -o logs/out_%J_%I.stdout
#BSUB -e logs/err_%J_%I.stderr


cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code_variants

# submit with 
# bsub < align_qc_var_bsub.sh
# max indices [1-81]

let i=$LSB_JOBINDEX
FQ_ID=$(tail -n+$i ../metadata/fq_prefix_list.txt | head -n1)
echo $FQ_ID

#################################### Trimming
## Settings: -W 4:00, -R mem=1000, -J embryo_trim)
# module purge
# module load trim_galore/0.4.5
# python align_qc.py --trimonly --fq $FQ_ID


#################################### HISAT2 alignment
## Settings: (-W 6:00, -R mem=20000, -J embryo_align_hisat2
module purge
module load fastqc/0.11.7
module load hisat2/2.0.5
module load python/3.5.0 py_packages/3.5
python align_qc.py --hisat2 --fq $FQ_ID
## HISAT2 is advertised as requiring only 3-4 Gb

#################################### STAR alignment
## Settings: -W 1:00, -R mem=40000, -J embryo_align_star
# module purge
# module load python/3.5.0
# module load star/2.6.1d
# python align_qc.py --star --fq $FQ_ID

#################################### Variant calling
## Settings: -W 6:00, -R mem=35000, -J embryo_var_gatk
# python var_call.py --fq $FQ_ID
