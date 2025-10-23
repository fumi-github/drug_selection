library(ggplot2)
library(dplyr)
library(data.table)

### If evaluating for original drug

train_trait = test_trait = "C10AA"
train_trait = test_trait = "C10AB"
train_trait = test_trait = "C10AX09"
x = readRDS(paste0(
  # "aucR2ROC_lipidaemiadrug_exclsupplement/auc.",
  "aucR2ROC_lipidaemiadrug_onlyaffectedexclsupplement/auc.",
  train_trait, ".PC.wowinsorize.rds"))

train_trait = test_trait = "C03"
train_trait = test_trait = "C07"
train_trait = test_trait = "C08"
train_trait = test_trait = "C09"
x = readRDS(paste0(
  # "aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.",
  "aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.",
  train_trait, ".PC.wowinsorize.rds"))

### If evaluating for different drug

train_trait = "C10AA";   test_trait = "C10AB"
train_trait = "C10AA";   test_trait = "C10AX09"
train_trait = "C10AB";   test_trait = "C10AA"
train_trait = "C10AB";   test_trait = "C10AX09"
train_trait = "C10AX09"; test_trait = "C10AA"
train_trait = "C10AX09"; test_trait = "C10AB"
x = readRDS(paste0(
  # "aucR2ROC_lipidaemiadrug_exclsupplement/auc.",
  "aucR2ROC_lipidaemiadrug_onlyaffectedexclsupplement/auc.",
  train_trait,
  ".test", test_trait, ".PC.wowinsorize.rds"))

train_trait = "C03"; test_trait = "C07"
train_trait = "C03"; test_trait = "C08"
train_trait = "C03"; test_trait = "C09"
train_trait = "C07"; test_trait = "C03"
train_trait = "C07"; test_trait = "C08"
train_trait = "C07"; test_trait = "C09"
train_trait = "C08"; test_trait = "C03"
train_trait = "C08"; test_trait = "C07"
train_trait = "C08"; test_trait = "C09"
train_trait = "C09"; test_trait = "C03"
train_trait = "C09"; test_trait = "C07"
train_trait = "C09"; test_trait = "C08"
x = readRDS(paste0(
  # "aucR2ROC_hypertensiondrug_DBPge100SBPge160.withHF/auc.",
  "aucR2ROC_hypertensiondrug_onlyaffectedDBPge100SBPge160/auc.",
  train_trait,
  ".test", test_trait, ".PC.wowinsorize.rds"))

# Under each stest independently, choose optimal hyperparameter glmnets

optglmnets = x %>%
  group_by(stest, glmnets) %>%
  summarize(meanpgsZ=mean(pgsZ, na.rm=TRUE), .groups="drop_last") %>%
  slice_max(meanpgsZ)
optglmnets
optdata = x %>%
  inner_join(optglmnets, join_by(stest, glmnets))

### Extract stats and plot

output = optdata %>% filter(strain==1) %>%
  summarize(
    pgsR2adjmean     =mean(pgsR2adj),      pgsR2adjsd     =sd(pgsR2adj),
    basalR2adjmean   =mean(basalR2adj),    basalR2adjsd   =sd(basalR2adj),
    basalpgsR2adjmean=mean(basalpgsR2adj), basalpgsR2adjsd=sd(basalpgsR2adj),
    encompassbasalpgsP=pnorm(-mean(encompassbasalpgsZadj))/2,
    pgsauc=mean(pgsaucadj),
    basalauc=mean(AUC.of.roc2adj),
    basalpgsauc=mean(AUC.of.roc1adj),
    aucp=median(padj))
output

m = optdata %>% filter(strain==1) %>%
  summarize(
    basal=   list(colMeans(do.call(rbind, basaldecileproportionadj))),
    pgs=     list(colMeans(do.call(rbind, pgsdecileproportionadj))),
    basalpgs=list(colMeans(do.call(rbind, basalpgsdecileproportionadj))))
s = optdata %>% filter(strain==1) %>%
  summarize(
    basal=   list(matrixStats::colSds(do.call(rbind, basaldecileproportionadj))),
    pgs=     list(matrixStats::colSds(do.call(rbind, pgsdecileproportionadj))),
    basalpgs=list(matrixStats::colSds(do.call(rbind, basalpgsdecileproportionadj))))
m = as.data.frame(do.call(cbind, lapply(m, unlist)))
m$decile = 1:10
m = m %>% tidyr::pivot_longer(cols = -c(decile), values_to="mean")
s = as.data.frame(do.call(cbind, lapply(s, unlist)))
s$decile = 1:10
s = s %>% tidyr::pivot_longer(cols = -c(decile), values_to="SD")
dataplot = left_join(m, s, join_by(decile, name))
dataplot$name = factor(dataplot$name,
                       levels=c("pgs", "basal", "basalpgs"),
                       labels=c("PGS", "Basic", "Combined"))
dataplot$train = train_trait
dataplot$test = test_trait
# dataplotcombine = rbind(dataplotcombine, dataplot)

ggplot(
  data = dataplot,
  # data = dataplotcombine[dataplotcombine$name=="PGS", ],
  aes(x = decile,
      y = mean)) +
  geom_col() +
  geom_linerange(
    aes(ymin = mean - SD,
        ymax = mean + SD)) +
  facet_grid(cols = vars(name)) +
  # facet_grid(cols = vars(train), rows=vars(test), scales="free_y") +
  xlab("Deciles") +
  ylab("Proportion taking drug") +
  scale_x_continuous(labels=NULL)
