#BSUB -W 6:00
#BSUB -q premium
#BUSB -n 3
#BSUB -R "rusage[mem=20000]"
#BSUB -P acc_schade01a
#BSUB -J "embryo_align_hisat2[4-5]"
#BSUB -m mothra
#BSUB -o logs/out_%J_%I.stdout
#BSUB -e logs/err_%J_%I.stderr


cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code

# submit with 
# bsub < align_qc_var_bsub.sh
# max indices [1-81]

let i=$LSB_JOBINDEX
FQ_ID=$(tail -n+$i ../metadata/fq_prefix_list.txt | head -n1)

## For trimming (-W 4:00, -R mem=1000, -J embryo_trim)
# module purge
# module load trim_galore/0.4.5
# python align_qc.py --trimonly --fq $FQ_ID

## For once trimming is done
module purge
module load fastqc/0.11.7
module load hisat2/2.0.5 star/2.6.1d
module load python/3.5.0 py_packages/3.5

## (-W 6:00, -R mem=20000, -J embryo_align_hisat2)
python align_qc.py --hisat2 --fq $FQ_ID
## for HISAT2 use 20Gb as it's advertised as requiring only 3-4 Gb

## (-W 6:00, -R mem=35000, -J embryo_align_hisat2)
# python align_qc.py --star --fq $FQ_ID

## for variant calling
# python var_call.py --fq $FQ_ID
