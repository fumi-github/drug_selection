library(tidyverse)
library(qqman)

data = data.table::fread(
  # "C10AA.delta_1.tbl"
  # "C10AA.logratio.100x_1.tbl"
  # "C10AA.myboxcox01.100x_1.tbl"
  # "C10AA.delta.SAMPLESIZE_1.tbl" # DON'T USE SAMPLESIZE; ignore heterogeneity
  # "C10AA.logratio.SAMPLESIZE_1.tbl"
  # "C10AA.myboxcox01.SAMPLESIZE_1.tbl"
  "geno01maf001mac0hwedev02/C10AA.delta_1.tbl"
  # "geno01maf001mac0hwedev02/C10AA.logratio.100x_1.tbl"
  # "geno01maf001mac0hwedev02/C10AA.myboxcoxv2.100x_1.tbl"
  # "geno01maf001mac0hwedev02/C03.myboxcoxv2.100x_1.tbl"
  # "geno01maf001mac0hwedev02/C07.myboxcoxv2.100x_1.tbl"
  # "geno01maf001mac0hwedev02/C08.myboxcoxv2.100x_1.tbl"
  # "geno01maf001mac0hwedev02/C09.myboxcoxv2.100x_1.tbl"
)
data <- data %>%
  rename(rsAminAmax = MarkerName,
         P = `P-value`) %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2),
         Z = Effect / StdErr,
         ALLELEMINOR = ifelse(Freq1 < 0.5, Allele1, Allele2),
         ZMINOR = ifelse(ALLELEMINOR == Allele1, Z, -Z),
         significant = P < 5 * 10^-8) %>%
  filter(HetDf == 1)

ukb = data.table::fread(
"../../ukb/gwas/regenie_hypolipidemics.LDLdelta.minafterstart28.wpower1/ukb_hypolipidemics.LDLdelta_step2_QT_chrall_statin.regenie.rs"
# "../../ukb/gwas/regenie_hypolipidemics.LDLlogratio.minafterstart28.wpower1/ukb_hypolipidemics.LDLlogratio_step2_QT_chrall_statin.regenie.rs"
# "../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox01.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox01_step2_QT_chrall_statin.regenie.rs"
)
ukb <- ukb %>%
  mutate(P = 10^(-LOG10P),
         Z = BETA / SE,
         ALLELEMINOR = ifelse(A1FREQ < 0.5, ALLELE1, ALLELE0),
         ZMINOR = ifelse(ALLELEMINOR == ALLELE1, Z, -Z),
         significant = P < 5 * 10^-8)

glgc = data.table::fread("jointGwasMc_LDL.txt.gz")
glgc <- glgc %>%
  mutate(A1 = toupper(A1),
         A2 = toupper(A2),
         `P-value` = as.numeric(`P-value`),
         Z = beta / se,
         rsAminAmax = paste0(rsid, ":", pmin(A1, A2), ":", pmax(A1, A2)),
         chr = as.numeric(sub("chr", "", sub(":.*", "", SNP_hg19))),
         pos = as.numeric(sub(".*:", "", SNP_hg19))) %>%
  arrange(chr, pos)

glgc <- glgc %>%
  filter(`P-value` < 5e-8) %>%
  filter(rsAminAmax %in% data$rsAminAmax)

### from locus.R
distance = 500 * 1000
snpstowindow = function(chr, pos, pow) {
  if (length(chr)==0) {
    return(data.frame(row.names=c("chr","start","end","peak","pow")))
  }
  result = c(chr[1], pos[1]);
  cprev = chr[1];
  pprev = pos[1];
  peak = pos[1];
  powmax  = pow[1];
  if (length(chr) > 1) {
    for (i in 2:length(chr)) {
      if (chr[i]==cprev & (pos[i]-pprev <= distance)) {
        pprev = pos[i];
        if (powmax < pow[i]) { peak = pos[i] }
        powmax = max(powmax, pow[i])
      } else {
        result = c(result, pprev, peak, powmax, chr[i], pos[i]);
        cprev = chr[i];
        pprev = pos[i];
        peak = pos[i];
        powmax = pow[i];
      }
    }
  }
  result = c(result, pprev, peak, powmax);
  result = data.frame(matrix(result, ncol=5, byrow=TRUE));
  names(result) = c("chr", "start", "end", "peak", "pow");
  result
}

x = snpstowindow(glgc$chr, glgc$pos, abs(glgc$Z))
glgc = glgc[
  match(paste0(x$chr, " ", x$peak), paste0(glgc$chr, " ", glgc$pos)), ]

x = intersect(data$rsAminAmax, glgc$rsAminAmax)
data2 = data[match(x, data$rsAminAmax), ]
glgc2 = glgc[match(x, glgc$rsAminAmax), ]
glgc2$ZMINOR = ifelse(data2$ALLELEMINOR == glgc2$A1, glgc2$Z, - glgc2$Z)
ggplot() +
  geom_point(
    aes(x = glgc2$ZMINOR,
        y = data2$ZMINOR),
    size = 0.75) +
  scale_color_manual(values=c("black", "red")) +
  xlab("Effect on LDL (Z-score)") +
  ylab("Effect on LDL change by statin\n(Z-score)")
  # ylab("Effect on log(LDL) change by statin\n(Z-score)")
data2[data2$significant, ]
glgc2[data2$significant, ]

### Manhattan
output = 
  # "C10AA.delta.tiff"
  "C10AA.logratio.100x.tiff"
  # "C10AA.myboxcoxv2.100x.tiff"
  # "C03.myboxcoxv2.100x.tiff"
  # "C07.myboxcoxv2.100x.tiff"
  # "C08.myboxcoxv2.100x.tiff"
  # "C09.myboxcoxv2.100x.tiff"
  
x = match(data$rsAminAmax, ukb$rsAminAmax)
data$CHROM = ukb$CHROM[x]
data$GENPOS = ukb$GENPOS[x]

tiff(filename=output, width=960, height=240)
f = manhattan(data, snp="rsAminAmax", chr="CHROM", bp="GENPOS", suggestiveline=FALSE)
dev.off()

x = sign(data$Freq1 - 0.5) * sign(data$Effect)
table(x)
median(data$`P-value`[x > 0])
median(data$`P-value`[x < 0])
min(data$`P-value`[x > 0])
min(data$`P-value`[x < 0])
