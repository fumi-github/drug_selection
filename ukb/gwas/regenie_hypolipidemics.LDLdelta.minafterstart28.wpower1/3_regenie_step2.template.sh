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
  --phenoFile ukb_hypolipidemics.LDLdelta_QT.txt \
  --phenoColList=ezetimibe,fibrate,statin \
  --covarFile ukb_hypolipidemics.LDL_covariates.txt \
  --qt \
  --pred ukb_hypolipidemics.LDLdelta_step1_QT_pred.list \
  --bsize 400 \
  --threads 8 \
  --out ukb_hypolipidemics.LDLdelta_step2_QT_chr${c}
