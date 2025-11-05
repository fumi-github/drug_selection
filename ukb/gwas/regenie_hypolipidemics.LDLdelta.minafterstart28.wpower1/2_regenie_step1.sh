#!/bin/bash
#SBATCH -p long
#SBATCH -w bhri-hpcn-08
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
  --phenoFile ukb_hypolipidemics.LDLdelta_QT.txt \
  --phenoColList=ezetimibe,fibrate,statin \
  --covarFile ukb_hypolipidemics.LDL_covariates.txt \
  --qt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /labs/workspace/users/ftakeuchi/regenie_tmp3_preds \
  --threads 8 \
  --out ukb_hypolipidemics.LDLdelta_step1_QT
