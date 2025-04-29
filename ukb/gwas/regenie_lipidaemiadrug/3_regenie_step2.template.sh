#!/bin/bash
#SBATCH -p long
#SBATCH -t 30-00:00:00
#SBATCH -c 8
#SBATCH --mem 6G

. /home/ftakeuchi/mambaforge/etc/profile.d/conda.sh
. /home/ftakeuchi/mambaforge/etc/profile.d/mamba.sh
mamba activate regenie_env

c=xxx
regenie \
  --step 2 \
  --pgen /home/ftakeuchi/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr${c}_v3_qc \
  --phenoFile ukb_lipidaemiadrug_BT.txt \
  --covarFile ukb_lipidaemiadrug_covariates.txt \
  --bt \
  --firth --approx --pThresh 0.01 \
  --pred ukb_lipidaemiadrug_step1_BT_pred.list \
  --bsize 400 \
  --threads 8 \
  --out ukb_lipidaemiadrug_step2_BT_chr${c}
