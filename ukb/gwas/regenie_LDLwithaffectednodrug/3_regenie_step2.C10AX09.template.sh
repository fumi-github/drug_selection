#!/bin/bash
#SBATCH -p long
#SBATCH -t 30-00:00:00
#SBATCH -c 8
#SBATCH --mem 7G

. /home/ftakeuchi/mambaforge/etc/profile.d/conda.sh
. /home/ftakeuchi/mambaforge/etc/profile.d/mamba.sh
mamba activate regenie_env

c=xxx
regenie \
  --step 2 \
  --pgen /home/ftakeuchi/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr${c}_v3_qc \
  --phenoFile ukb_LDL_QT.txt \
  --covarFile ukb_LDLwithaffectednodrug_covariates.txt \
  --catCovarList C10AA,C10AB,C10AC,C10AD,C10AX06,C10AX09,affectednodrug \
  --qt \
  --pred ukb_LDLwithaffectednodrug_step1_QT_pred.list \
  --bsize 400 \
  --threads 8 \
  --interaction C10AX09[0] \
  --out ukb_LDLwithaffectednodrug.C10AX09_step2_QT_chr${c}
