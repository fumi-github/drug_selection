library(dplyr)
library(data.table)

# Load genetic data for chromosome 1.
genotype_data = read.table("ukb_imp_chr1_v3_qc.snps.raw", header=TRUE, sep="\t")
# genotype_data = read.table("ukb_imp_chr1_v3_qc.withHF.snps.raw", header=TRUE, sep="\t")
genotype_data = genotype_data[, -c(3:6)] # Remove PAT, MAT, SEX, PHENOTYPE columns.

# Loop to load and merge genetic data for chromosomes 2 through 22.
for (i in 2:22) {
  print(i)
  chr_data = read.table(paste0("ukb_imp_chr", i, "_v3_qc.snps.raw"), header=TRUE, sep="\t")
  # chr_data = read.table(paste0("ukb_imp_chr", i, "_v3_qc.withHF.snps.raw"), header=TRUE, sep="\t")
  
  # Verify that sample IDs (FID) are identical before merging.
  print(identical(genotype_data$FID, chr_data$FID))
  
  # Column-bind the new chromosome data, excluding its sample ID columns.
  genotype_data = cbind(genotype_data, chr_data[, -c(1:6)])
}
rsallele = colnames(genotype_data)[-c(1:2)]

# a = readRDS("aucR2ROC_lipidaemiadrug_exclsupplement/auc.C10AA.PC.wowinsorize.rds") %>% filter(glmnets==0.006)
# a = readRDS("aucR2ROC_lipidaemiadrug_exclsupplement/auc.C10AB.PC.wowinsorize.rds") %>% filter(glmnets==0.08)
# a = readRDS("aucR2ROC_lipidaemiadrug_exclsupplement/auc.C10AX09.PC.wowinsorize.rds") %>% filter(glmnets==0.03)
a = readRDS("aucR2ROC_lipidaemiadrug_onlyaffectedexclsupplement/auc.C10AA.PC.wowinsorize.rds") %>% filter(glmnets==0.006)
a = readRDS("aucR2ROC_lipidaemiadrug_onlyaffectedexclsupplement/auc.C10AB.PC.wowinsorize.rds") %>% filter(glmnets==0.09)
a = readRDS("aucR2ROC_lipidaemiadrug_onlyaffectedexclsupplement/auc.C10AX09.PC.wowinsorize.rds") %>% filter(glmnets==0.05) 

# a = readRDS("aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.C03.PC.wowinsorize.rds") %>% filter(glmnets==0.01)
# a = readRDS("aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.C07.PC.wowinsorize.rds") %>% filter(glmnets==0.005)
# a = readRDS("aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.C08.PC.wowinsorize.rds") %>% filter(glmnets==0.004)
# a = readRDS("aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.C09.PC.wowinsorize.rds") %>% filter(glmnets==0.002)
a = readRDS("aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.C03.PC.wowinsorize.rds") %>% filter(glmnets==0.01)
a = readRDS("aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.C07.PC.wowinsorize.rds") %>% filter(glmnets==0.002)
a = readRDS("aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.C08.PC.wowinsorize.rds") %>% filter(glmnets==0.004)
a = readRDS("aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.C09.PC.wowinsorize.rds") %>% filter(glmnets==0.003)

i="C09"

x = a %>% select(stest, coef2, pgscoeff)
x = x %>% group_by(stest) %>%
  summarize(coef2avg = list(colMeans(do.call(rbind, coef2))),
            pgscoeffavg = mean(pgscoeff, na.rm=TRUE)) #2025.09.29; sum(pgscoeff, na.rm=TRUE) / n())
x = do.call(rbind, x$coef2avg) * x$pgscoeffavg
x = colMeans(x)
res = as.data.frame(x)
res$rs = sub("_.*", "", rsallele)
res$allele = sub(".*_", "", rsallele)
write.table(res[, c(2:3, 1)], file=paste0("coeff.", i, ".txt"),
            row.names=F, col.names=F, sep=" ", quote=FALSE)

x = a %>% select(coef1)
x = colMeans(do.call(rbind, x[[1]]))
res = as.data.frame(x)
write.table(res, file=paste0("coeff_basal.", i, ".txt"),
            row.names=T, col.names=F, sep=" ", quote=FALSE)

### aou

snpsukbaou = read.table("snps_ukb_aou.txt", sep="\t", header=TRUE)
# snpsukbaou = read.table("snps_ukb_aou.withHF.txt", sep="\t", header=TRUE)

res = read.table(paste0("coeff.", i, ".txt"), sep=" ")

identical(sort(snpsukbaou$ID), sort(res$V1))
res$IDaou = snpsukbaou$IDaou[match(res$V1, snpsukbaou$ID)]

# Invert rs7948851 G/C !!!
x = (res$V1 %in% c("rs7948851"))
res$V3[x] = - res$V3[x] 

write.table(res[, c(4, 2:3)], file=paste0("coeff.", i, ".hg38.txt"),
            row.names=F, col.names=F, quote=F)
