# Choose hyperparameter glmnets that maximizes pgsZ obtained from validation data

library(data.table)
library(dplyr)
library(ggplot2)

# data = readRDS("aucR2ROC_lipidaemiadrug_exclsupplement/auc.C10AA.PC.wowinsorize.rds")
data = readRDS("aucR2ROC_lipidaemiadrug_onlyaffectedexclsupplement/auc.C10AA.PC.wowinsorize.rds")
# data = readRDS("aucR2ROC_hypertensiondrug_DBPge90SBPge140/auc.C02.PC.wowinsorize.pilot.rds")
# data = readRDS("aucR2ROC_hypertensiondrug_DBPge90SBPge140.withHF/auc.C02.PC.wowinsorize.rds")
# data = readRDS("aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.C03.PC.wowinsorize.rds")
data = readRDS("aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.C03.PC.wowinsorize.rds")

data$pgsauc = as.numeric(data$pgsauc)
plotdata = data %>%
  select(trait, glmnets, stest, strain, pgsZ, pgsZavg) %>%
  tidyr::pivot_longer(cols = c(pgsZ, pgsZavg)) %>%
  filter(name=="pgsZ" | strain==1)
# select(trait, glmnets, stest, strain, basalpgsZ, basalpgsZavg, basalpgsZadj) %>%
# tidyr::pivot_longer(cols = c(basalpgsZ, basalpgsZavg, basalpgsZadj)) %>%
# filter(name=="basalpgsZ" | strain==1)
# select(trait, glmnets, stest, strain, pgsauc, pgsaucavg) %>%
# tidyr::pivot_longer(cols = c(pgsauc, pgsaucavg)) %>%
# filter(name=="pgsauc" | strain==1)
# select(trait, glmnets, stest, strain, AUC.of.roc1, AUC.of.roc1avg, AUC.of.roc1adj) %>%
#   tidyr::pivot_longer(cols = c(AUC.of.roc1, AUC.of.roc1avg, AUC.of.roc1adj)) %>%
#   filter(name=="AUC.of.roc1" | strain==1)
# select(trait, glmnets, stest, strain, P, p, Pavg, pavg, Padj, padj) %>%
#   tidyr::pivot_longer(cols = c(P, p, Pavg, pavg, Padj, padj)) %>%
#   filter(name=="P" | strain==1)
x = plotdata %>% group_by(name, glmnets) %>%
  summarize(mean = mean(value, na.rm=TRUE), .groups="drop")

### aucR2ROC_lipidaemiadrug_exclsupplement
# C10AA.PC.winsorize     pgsZ peaks at glmnets=0.005 15.977149; basalpgsZadj 53.88812; AUC.of.roc1adj 0.7012322
# C10AA.PC.wowinsorize   pgsZ peaks at glmnets=0.006 15.914570; basalpgsZadj 53.87839; AUC.of.roc1adj 0.7009527
# C10AA.woPC.winsorize   pgsZ peaks at glmnets=0.005 15.966190; basalpgsZadj 53.88975; AUC.of.roc1adj 0.7012340
# C10AA.woPC.wowinsorize pgsZ peaks at glmnets=0.006 15.907954; basalpgsZadj 53.96538; AUC.of.roc1adj 0.7012292
# C10AB.PC.winsorize     pgsZ peaks at glmnets=0.02  7.488878; basalpgsZadj 9.295842;  AUC.of.roc1adj 0.6501708
# C10AB.PC.wowinsorize   pgsZ peaks at glmnets=0.08  7.509757; basalpgsZadj 10.280583; AUC.of.roc1adj 0.6643472
# C10AB.woPC.winsorize   pgsZ peaks at glmnets=0.03  7.456332; basalpgsZadj 9.376743;  AUC.of.roc1adj 0.6518369
# C10AB.woPC.wowinsorize pgsZ peaks at glmnets=0.08  7.502536; basalpgsZadj 10.271698; AUC.of.roc1adj 0.6639517
# C10AX09.PC.winsorize     pgsZ peaks at glmnets=0.02  4.241844; basalpgsZadj 10.843850; AUC.of.roc1adj 0.6286024
# C10AX09.PC.wowinsorize   pgsZ peaks at glmnets=0.03  4.171495; basalpgsZadj 10.796423; AUC.of.roc1adj 0.6281153
# C10AX09.woPC.winsorize   pgsZ peaks at glmnets=0.02  4.213890; basalpgsZadj 10.822441; AUC.of.roc1adj 0.6283636
# C10AX09.woPC.wowinsorize pgsZ peaks at glmnets=0.03  4.172436; basalpgsZadj 10.599086; AUC.of.roc1adj 0.6239227
# PC vs woPC almost identical
# Padj is slightly larget than padj
### aucR2ROC_hypertensiondrug_DBPge90SBPge140
# C02.PC.wowinsorize   pgsZ peaks at glmnets=0.01 2.09018413; basalpgsZadj ; AUC.of.roc1adj 
# C03.PC.wowinsorize   pgsZ peaks at glmnets=0.007 10.2453293; basalpgsZadj ; AUC.of.roc1adj 
# C07.PC.wowinsorize   pgsZ peaks at glmnets=0.004 8.517799; basalpgsZadj ; AUC.of.roc1adj 
# C08.PC.wowinsorize   pgsZ peaks at glmnets=0.004 12.268201; basalpgsZadj ; AUC.of.roc1adj 
# C09.PC.wowinsorize   pgsZ peaks at glmnets=0.003 16.239628; basalpgsZadj ; AUC.of.roc1adj 
### aucR2ROC_hypertensiondrug_DBPge90SBPge140.withHF
# C02.PC.wowinsorize   pgsZ peaks at glmnets=0.009(probably larger; not scanned) 2.057313; basalpgsZadj ; AUC.of.roc1adj 
# C03.PC.wowinsorize   pgsZ peaks at glmnets=0.007 10.221339; basalpgsZadj ; AUC.of.roc1adj 
# C07.PC.wowinsorize   pgsZ peaks at glmnets=0.004 8.615575; basalpgsZadj ; AUC.of.roc1adj 
# C08.PC.wowinsorize   pgsZ peaks at glmnets=0.005 12.15377; basalpgsZadj ; AUC.of.roc1adj 
# C09.PC.wowinsorize   pgsZ peaks at glmnets=0.003 16.25736; basalpgsZadj ; AUC.of.roc1adj 
### aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF
# C02.PC.wowinsorize   pgsZ peaks at glmnets=0.01 1.24518742; basalpgsZadj 18.71258; AUC.of.roc1adj 0.6600818
# C03.PC.wowinsorize   pgsZ peaks at glmnets=0.01 5.585021; basalpgsZadj 37.26141; AUC.of.roc1adj 0.6536437
# C07.PC.wowinsorize   pgsZ peaks at glmnets=0.005 6.171682; basalpgsZadj 20.04089; AUC.of.roc1adj 0.5885716
# C08.PC.wowinsorize   pgsZ peaks at glmnets=0.004 7.766874; basalpgsZadj 29.12031; AUC.of.roc1adj 0.6246002
# C09.PC.wowinsorize   pgsZ peaks at glmnets=0.002 9.445098; basalpgsZadj 36.03357; AUC.of.roc1adj 0.6252468
## aucR2ROC_hypertensiondrug.withHF
# C02.PC.wowinsorize   pgsZ peaks at glmnets=0.009(probably larger; not scanned) 1.2093045; basalpgsZadj ; AUC.of.roc1adj 
# C03.PC.wowinsorize   pgsZ peaks at glmnets=0.009(probably larger; not scanned) 4.428036; basalpgsZadj ; AUC.of.roc1adj 
# C07.PC.wowinsorize   pgsZ peaks at glmnets=0.007 5.408616; basalpgsZadj ; AUC.of.roc1adj 
# C08.PC.wowinsorize   pgsZ peaks at glmnets=0.004 6.779372; basalpgsZadj ; AUC.of.roc1adj 
# C09.PC.wowinsorize   pgsZ peaks at glmnets=0.003 7.508944; basalpgsZadj ; AUC.of.roc1adj 

ggplot(data = plotdata,
       aes(x = glmnets,
           y = value)) +
  geom_point(size = 0.1) +
  geom_smooth() +
  scale_x_log10() +
  facet_grid(cols = vars(name))
View(x)
