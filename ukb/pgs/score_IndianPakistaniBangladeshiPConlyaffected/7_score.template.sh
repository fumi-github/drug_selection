#!/bin/bash
#SBATCH -p standard
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem 32G

module load plink2

for c in xxx
do
	 
for p in ../coeff/coeff.C10??.txt ../coeff/coeff.C10????.txt
do
pradix=${p/.txt}
pradix=${pradix/..\/coeff\/}
plink2 --score $p list-variants cols=scoresums \
       --keep ukb_lipidaemiadrug_onlyaffectedexclsupplement_BT.IndianPakistaniBangladeshiPC.txt \
       --out ${pradix}.chr$c \
       --pfile /home/ftakeuchi/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr${c}_v3_qc \
       --threads 8 \
       --memory 32000
done

for p in ../coeff/coeff.C0?.txt
do
pradix=${p/.txt}
pradix=${pradix/..\/coeff\/}
plink2 --score $p list-variants cols=scoresums \
       --keep ukb_hypertensiondrug_onlyaffectedDBPge100SBPge160_BT.IndianPakistaniBangladeshiPC.txt \
       --out ${pradix}.chr$c \
       --pfile /home/ftakeuchi/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr${c}_v3_qc \
       --threads 8 \
       --memory 32000
done

done
