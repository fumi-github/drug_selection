#!/bin/bash
#SBATCH -p long
#SBATCH -t 30-00:00:00
#SBATCH -c 8
#SBATCH --mem 32G

module load plink2

p=xxx # C02
c=yyy

plink2 --logistic \
       --out ${p}.chr${c} \
       --pfile /home/ftakeuchi/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr${c}_v3_qc \
       --pheno ../regenie_lipidaemiadrug_onlyaffectedexclsupplement/ukb_lipidaemiadrug_onlyaffectedexclsupplement_BT.txt \
       --pheno-name ${p} \
       --1 \
       --covar ../regenie_lipidaemiadrug_onlyaffectedexclsupplement/ukb_lipidaemiadrug_onlyaffectedexclsupplement_covariates.txt \
       --covar-name sexM age bmi agec2 PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20 \
       --covar-variance-standardize \
       --maf 0.01 \
       --threads 8 \
       --memory 32000
