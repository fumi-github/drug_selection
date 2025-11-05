library(tidyverse)

data = read.table("summary.rg.txt",header=T)
data$p1 = sub(".sumstats.gz", "", data$p1)
# data$p1 = sub("lipidaemiadrug.", "", data$p1)
data$p2 = sub(".sumstats.gz", "", data$p2)
# data$p2 = sub("lipidaemiadrug.", "", data$p2)

data = data[data$p1 != "p1", ]
data$rg = as.numeric(data$rg)
data$se = as.numeric(data$se)
data$z  = as.numeric(data$z)
data$p  = as.numeric(data$p)

target = c(
  "lipidaemiadrug_exclsupplement_plink.C10AA",   "C10AA",
  "lipidaemiadrug_exclsupplement_plink.C10AB",   "C10AB",
  "lipidaemiadrug_exclsupplement_plink.C10AX09", "C10AX09",
  "glgc-lipids2021_LDL",                         "LDL",
  "glgc-lipids2021_HDL",                         "HDL",
  "glgc-lipids2021_logTG",                       "logTG",
  "keaton2024_SBP",                              "SBP",
  "keaton2024_DBP",                              "DBP",
  "UKBB.GWAS1KG.EXOME.CAD",                      "CAD",
  "suzuki2024_T2D",                              "T2D"
)
target = c(
  "lipidaemiadrug_onlyaffectedexclsupplement_plink.C10AA",   "C10AA",
  "lipidaemiadrug_onlyaffectedexclsupplement_plink.C10AB",   "C10AB",
  "lipidaemiadrug_onlyaffectedexclsupplement_plink.C10AX09", "C10AX09",
  "glgc-lipids2021_LDL",                         "LDL",
  "glgc-lipids2021_HDL",                         "HDL",
  "glgc-lipids2021_logTG",                       "logTG",
  "keaton2024_SBP",                              "SBP",
  "keaton2024_DBP",                              "DBP",
  "UKBB.GWAS1KG.EXOME.CAD",                      "CAD",
  "suzuki2024_T2D",                              "T2D"
)
target = c(
  "hypertensiondrug_DBPge100SBPge160_plink.C03", "C03",
  "hypertensiondrug_DBPge100SBPge160_plink.C07", "C07",
  "hypertensiondrug_DBPge100SBPge160_plink.C08", "C08",
  "hypertensiondrug_DBPge100SBPge160_plink.C09", "C09",
  "glgc-lipids2021_LDL",                         "LDL",
  "glgc-lipids2021_HDL",                         "HDL",
  "glgc-lipids2021_logTG",                       "logTG",
  "keaton2024_SBP",                              "SBP",
  "keaton2024_DBP",                              "DBP",
  "UKBB.GWAS1KG.EXOME.CAD",                      "CAD",
  "nielsenAFib",                                 "AF"
)
target = c(
  "hypertensiondrug_onlyaffectedDBPge100SBPge160_plink.C03", "C03",
  "hypertensiondrug_onlyaffectedDBPge100SBPge160_plink.C07", "C07",
  "hypertensiondrug_onlyaffectedDBPge100SBPge160_plink.C08", "C08",
  "hypertensiondrug_onlyaffectedDBPge100SBPge160_plink.C09", "C09",
  "glgc-lipids2021_LDL",                         "LDL",
  "glgc-lipids2021_HDL",                         "HDL",
  "glgc-lipids2021_logTG",                       "logTG",
  "keaton2024_SBP",                              "SBP",
  "keaton2024_DBP",                              "DBP",
  "UKBB.GWAS1KG.EXOME.CAD",                      "CAD"
)
target = as.data.frame(matrix(target, nrow=length(target)/2, ncol=2, byrow=TRUE))
colnames(target) = c("filename", "Trait")

data = data[data$p1 %in% target$filename, ]
data = data[data$p2 %in% target$filename, ]
data$p1 = target$Trait[match(data$p1, target$filename)]
data$p2 = target$Trait[match(data$p2, target$filename)]

# x = sort(unique(c(data$p1,data$p2)))
x = target$Trait
data$p1 = factor(data$p1, levels=x)
data$p2 = factor(data$p2, levels=x)

###
library(corrplot)
x = data[, c("p2", "p1", "rg")]; colnames(x)[1:2] = c("p1", "p2")
x = rbind(x, data[, c("p1", "p2", "rg")])
x = x %>% pivot_wider(names_from=p2, values_from=rg)
M = as.matrix(x[, -1])
rownames(M) = x$p1
M = M[target$Trait, target$Trait]

x = data[, c("p2", "p1", "p")]; colnames(x)[1:2] = c("p1", "p2")
x = rbind(x, data[, c("p1", "p2", "p")])
x = x %>% pivot_wider(names_from=p2, values_from=p)
P = as.matrix(x[, -1])
rownames(P) = x$p1
P = P[target$Trait, target$Trait]

pl=c(colorRampPalette(c("blue", "white"))(60)[c(1, 21, 31, 41, 51)],
     colorRampPalette(c("white", "red"))(60)[c(10, 20, 30, 40, 60)])
p = corrplot(M, p.mat=P, method='color', type='upper', diag=FALSE,
             insig="n", col=pl)
x = p$corrPos
x$corr[x$p.value >= 0.01] = NA
text(p$corrPos$x, p$corrPos$y, round(x$corr, 2), cex=0.7)
# FIRST EXPORT AS EPS!!

###
library(ggcorrplot)
data <- data %>%
  mutate(Significance = ifelse(p < 0.05, "*", ""))

ggplot(data, aes(x = p1, y = p2, fill = rg)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation\nCoefficient") +
  geom_text(aes(label = paste0(round(rg, 2), Significance)), 
            color = "black", size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation Matrix Heatmap with Significance", x = "", y = "")

### https://github.com/mkanai/ldsc-corrplot-rg
library(corrplot)
library(reshape2)
corrplot_nsquare = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  x2 = dcast(rg, p1 ~ p2, value.var = "rg")
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1
  
  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0
  
  x2 = dcast(rg, p1 ~ p2, value.var = "q")
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1
  
  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = NA
  }
  
  trait1 = match(trait1, rownames(mat2))
  trait2 = match(trait2, colnames(mat2))
  mat2 = mat2[trait1, trait2]
  qmat2 = qmat2[trait1, trait2]
  
  corrplot(mat2, method = "square", order = order,
           p.mat = qmat2, sig.level = c(0, 0.05), insig = "label_sig",
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")
  return(list(mat = mat2, qmat = qmat2))
}

data$p1_category = 1
data$p2_category = 1
data$q = p.adjust(data$p, method="fdr")
rg = data

# duplicate lines
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
tmp$p1_category = rg$p2_category
tmp$p2_category = rg$p1_category
rg = rbind(rg, tmp)

corrplot_nsquare(rg, target$Trait, target$Trait, landscape=TRUE)
