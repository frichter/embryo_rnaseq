#BSUB -W 4:00
#BSUB -q premium
#BUSB -n 3
#BSUB -R "rusage[mem=5000]"
#BSUB -P acc_schade01a
#BSUB -J "embryo_trim[10-12]"
#BSUB -m mothra
#BSUB -o logs/out_%J_%I.stdout
#BSUB -e logs/err_%J_%I.stderr


# submit with bsub < align_qc_var_bsub.sh
# max indices either [1-81] or [1-82]

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/code

let i=$LSB_JOBINDEX
let i=10
FQ_ID=$(tail -n+$i ../metadata/fq_prefix_list.txt | head -n1)

## For trimming (-W 4:00, -R mem=5000, -J embryo_trim_)
module purge
module load trim_galore/0.4.5
python align_qc.py --trimonly --fq $FQ_ID

## For once trimming is done
## (-W :00, -R mem=50000, -J embryo_align_)
# module purge
# module load fastqc/0.11.7
# module load hisat2/2.0.5 star/2.6.1d
# module load python/3.5.0 py_packages/3.5
# python align_qc.py --fq $FQ_ID

## for variant calling
# python var_call.py --fq $FQ_ID
