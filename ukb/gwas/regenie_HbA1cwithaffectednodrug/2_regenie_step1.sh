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
  --extract ../qc_pass.snplist \
  --keep ../qc_pass.id \
  --phenoFile ukb_HbA1c_QT.txt \
  --covarFile ukb_HbA1cwithaffectednodrug_covariates.txt \
  --catCovarList A10A,A10BA,A10BB,A10BD,A10BF,A10BG,A10BX,affectednodrug \
  --qt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /labs/sysgen/workspace/users/ftakeuchi/regenie_tmp_preds \
  --threads 8 \
  --out ukb_HbA1cwithaffectednodrug_step1_QT
