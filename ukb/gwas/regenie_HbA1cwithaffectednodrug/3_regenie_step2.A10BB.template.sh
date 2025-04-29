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
  --phenoFile ukb_HbA1c_QT.txt \
  --covarFile ukb_HbA1cwithaffectednodrug_covariates.txt \
  --catCovarList A10A,A10BA,A10BB,A10BD,A10BF,A10BG,A10BX,affectednodrug \
  --qt \
  --pred ukb_HbA1cwithaffectednodrug_step1_QT_pred.list \
  --bsize 400 \
  --threads 8 \
  --interaction A10BB[0] \
  --out ukb_HbA1cwithaffectednodrug.A10BB_step2_QT_chr${c}
