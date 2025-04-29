#!/bin/bash
#Partition: 'sysgen' (<24h) or 'sysgen_long' (>24h)
#SBATCH -p long
#Give the job a name:
#SBATCH --job-name="QC-Gen"
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# The amount of memory in megabytes per process in the job (try to put the smaller file first):
#SBATCH --mem=20048
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2-0:0:00
# The job command(s):

/software/plink2/v2.00a2.3LM/plink2 \
  --bfile /home/ftakeuchi/ukb55469/Genetic_data/Genotype/Plink/ukb_cal_allChrs \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --write-snplist --write-samples --no-id-header \
  --out qc_pass
