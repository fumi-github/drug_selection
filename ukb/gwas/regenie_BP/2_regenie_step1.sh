#!/bin/bash
#SBATCH -p gpu
#SBATCH -t 30-00:00:00
#SBATCH -c 8
#SBATCH --mem 12G

. /home/ftakeuchi/mambaforge/etc/profile.d/conda.sh
. /home/ftakeuchi/mambaforge/etc/profile.d/mamba.sh
mamba activate regenie_env

regenie \
  --step 1 \
  --bed /home/ftakeuchi/ukb55469/Genetic_data/Genotype/Plink/ukb_cal_allChrs \
  --extract qc_pass.snplist \
  --keep qc_pass.id \
  --phenoFile ukb_BP_QT.txt \
  --covarFile ukb_BP_covariates.txt \
  --catCovarList C02,C03,C07,C08,C09 \
  --qt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /labs/sysgen/workspace/users/ftakeuchi/regenie_tmp_preds \
  --threads 8 \
  --out ukb_BP_step1_QT
