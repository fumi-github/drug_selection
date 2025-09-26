#!/bin/bash
#SBATCH -p standard
#SBATCH -c 4
#SBATCH --mem 32G

module load plink2

c=yyy

plink2 --export A \
       --out ukb_imp_chr${c}_v3_qc.withHF.snps \
       --pfile /home/ftakeuchi/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr${c}_v3_qc \
       --extract snps_rs.chr${c}.withHF.txt \
       --export-allele snps_rs_allele.chr${c}.withHF.txt \
       --threads 4 \
       --memory 32000
