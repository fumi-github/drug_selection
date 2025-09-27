library(glmnet)
library(ridge)
library(ggplot2)
library(dplyr)
library(data.table)
library(pROC)
library(R2ROC)

winsorize = function (x) {
  xmax = quantile(x, 0.99, na.rm = TRUE)
  xmin = quantile(x, 0.01, na.rm = TRUE)
  x[!is.na(x) & x > xmax] = xmax
  x[!is.na(x) & x < xmin] = xmin
  return(x)
}

data = read.table("ukb_imp_chr1_v3_qc.snps.raw", header=TRUE, sep="\t")
# data = read.table("ukb_imp_chr1_v3_qc.withHF.snps.raw", header=TRUE, sep="\t")
data = data[, -c(3:6)]
for (i in 2:22) {
  print(i)
  foo = read.table(paste0("ukb_imp_chr", i, "_v3_qc.snps.raw"), header=TRUE, sep="\t")
  # foo = read.table(paste0("ukb_imp_chr", i, "_v3_qc.withHF.snps.raw"), header=TRUE, sep="\t")
  print(identical(data$FID, foo$FID))
  data = cbind(data, foo[, -c(1:6)])
}

pheno = read.table(
  # "../gwas/regenie_lipidaemiadrug/ukb_lipidaemiadrug_BT.txt", #20240723
  # "../gwas/regenie_lipidaemiadrug_HDLlt40LDLge190/ukb_lipidaemiadrug_HDLlt40LDLge190_BT.txt", #20240819
  "../gwas/regenie_lipidaemiadrug_exclsupplement/ukb_lipidaemiadrug_exclsupplement_BT.txt", #20240827
  # "../gwas/regenie_hypertensiondrug_DBPge90SBPge140/ukb_hypertensiondrug_DBPge90SBPge140_BT.txt", #20241018
  # "../gwas/regenie_diabetesdrug_HbA1cge48/ukb_diabetesdrug_HbA1cge48_BT.txt", #20241018
  header=TRUE, sep=" ")
covariates = read.table(
  # "../gwas/regenie_lipidaemiadrug/ukb_lipidaemiadrug_covariates.txt", #20240723
  # "../gwas/regenie_lipidaemiadrug_HDLlt40LDLge190/ukb_lipidaemiadrug_HDLlt40LDLge190_covariates.txt", #20240819
  "../gwas/regenie_lipidaemiadrug_exclsupplement/ukb_lipidaemiadrug_exclsupplement_covariates.txt", #20240827
  # "../gwas/regenie_hypertensiondrug_DBPge90SBPge140/ukb_hypertensiondrug_DBPge90SBPge140_covariates.txt", #20241018
  # "../gwas/regenie_diabetesdrug_HbA1cge48/ukb_diabetesdrug_HbA1cge48_covariates.txt", #20241018
  header=TRUE, sep=" ")
identical(pheno$FID, covariates$FID)
data = data[match(pheno$FID, data$FID), ]
data = data[, -c(1:2)]
pheno = pheno[, -c(1:2)]
covariates = covariates[, -c(1:2)]


### split train test
splittrainvalidtest = function(n, ptest=0.2, stest=1, ptrain=0.6, strain=1) {
  set.seed(stest)
  x = sample(n)
  test = sort(head(x, round(n*ptest)))
  trainvalid = tail(x, n - round(n*ptest))
  set.seed(strain)
  trainvalid = sample(trainvalid)
  train = sort(head(trainvalid, round(n*ptrain)))
  valid = sort(tail(trainvalid, length(trainvalid) - round(n*ptrain)))
  return(list(train=train, valid=valid, test=test))
}

glm2steps = function(trait, covariates, data, glmnets=0.04) {
  a0 = glm(trait ~ .,
           data=covariates,
           family=binomial(),
           na.action=na.exclude)
  y = naresid(a0$na.action, a0$residual)
  # y = winsorize(y)
  
  keep = (rowSums(is.na(data)) == 0 & !is.na(y))
  table(keep)
  y2 = y[keep]
  data2 = data[keep, ]
  for (i in 1:ncol(data2)) { data2[, i] = data2[, i] - mean(data2[, i])}
  
  a2 = glmnet(
    x = as.matrix(data2),
    y = y2)
  coefSNP = as.numeric(coef(a2, s=glmnets))[-1]
  
  # # relaxed LASSO
  # a3 = lm(
  #   y ~ x,
  #   data=list(x = as.matrix(data2) %*% diag((coefSNP!=0)* 1),
  #   y = y2))
  # coefSNPb = coef(a3)[-1]
  # coefSNPb[is.na(coefSNPb)] = 0
  
  # a2 = linearRidge(
  #   y ~ .,
  #   data = cbind(data.frame(y = y2), data2))
  # coefSNP = coef(a2)[-1]
  
  # # 1 step; NG aborts
  # x = cbind(data.frame(y = phenotrain$trait), datatrain)
  # x = x[rowSums(is.na(x)) == 0, ]
  # a2 = logisticRidge(
  #   y ~ .,
  #   data = x,
  #   # scaling="corrForm")
  #   scaling="none")
  
  
  # # rescale
  # a3 = lm(y ~ x, data=data.frame(y = y2, x = as.matrix(data2) %*% coefSNP))
  # coefSNP = coefSNP * coef(a3)[2]
  
  return(list(coef1=coef(a0), coef2=coefSNP))
}

decilecut = function(x) {
  q = quantile(x, seq(0, 1, 0.1), na.rm=TRUE)
  if (length(unique(q)) < 11) {
    delta = min(diff(sort(unique(x)))) / 2
    set.seed(1)
    x = x + runif(length(x), max=delta)
    q = quantile(x, seq(0, 1, 0.1), na.rm=TRUE)
  }
  return(cut(x, q, include.lowest=TRUE))
}
decileproportion = function(x, trait) {
  x = data.frame(decile = decilecut(x), trait = trait) %>%
    na.omit() %>%
    group_by(decile) %>%
    summarize(proportion = sum(trait=="1") / n(), .groups="drop")
  return(x$proportion)
}

assess2steps = function(trait, covariates, data, coef1, coef2) {
  combined = data.frame(
    trait = trait,
    basal = as.matrix(cbind(1, covariates)) %*% coef1,
    pgs = as.matrix(data) %*% coef2)
  combined = combined[rowSums(is.na(combined)) == 0, ]
  combined$basalpgs = combined$basal + combined$pgs
  x = summary(a1 <- glm(factor(trait) ~ pgs, data=combined, family=binomial()))
  pgsZ = ifelse(nrow(x$coefficients)>=2, x$coefficients[2,3], NA)
  pgsR2 = anova(a1)[2,2] / anova(a1)[1,4]
  x = summary(a2 <- glm(factor(trait) ~ basal, data=combined, family=binomial()))
  basalZ = ifelse(nrow(x$coefficients)>=2, x$coefficients[2,3], NA)
  basalR2 = anova(a2)[2,2] / anova(a2)[1,4]
  x = summary(a3 <- glm(factor(trait) ~ basalpgs, data=combined, family=binomial()))
  basalpgsZ = ifelse(nrow(x$coefficients)>=2, x$coefficients[2,3], NA)
  basalpgsR2 = anova(a3)[2,2] / anova(a3)[1,4]
  x = summary(a4 <- glm(factor(trait) ~ basal + pgs, data=combined, family=binomial()))
  basalcoeff = ifelse(nrow(x$coefficients)>=2, x$coefficients[2,1], NA)
  pgscoeff = ifelse(nrow(x$coefficients)>=3, x$coefficients[3,1], NA)
  x = summary(a5 <- glm(factor(trait) ~ basal + basalpgs, data=combined, family=binomial()))
  encompassbasalpgsZ = ifelse(nrow(x$coefficients)>=3, x$coefficients[3,3], NA)
  res = data.table(
    pgsZ=pgsZ, pgsR2=pgsR2,
    basalZ=basalZ, basalR2=basalR2,
    basalpgsZ=basalpgsZ, basalpgsR2=basalpgsR2,
    basalcoeff=basalcoeff, pgscoeff=pgscoeff, encompassbasalpgsZ=encompassbasalpgsZ,
    basaldecileproportion    = list(decileproportion(combined$basal, combined$trait)),
    pgsdecileproportion      = list(decileproportion(combined$pgs, combined$trait)),
    basalpgsdecileproportion = list(decileproportion(combined$basalpgs, combined$trait)))
  
  a3 = pROC::roc.test(
    response=combined$trait, 
    predictor1=as.numeric(combined$basalpgs),
    predictor2=as.numeric(combined$basal))
  res2 = c(a3$estimate, a3$statistic, a3$p.value)
  names(res2)[4] = "P"
  res2 = data.frame(as.list(res2))
  res = cbind(res, res2)
  
  x = pROC::roc(response=combined$trait, predictor=as.numeric(combined$pgs))
  res$pgsauc = x$auc
  res$pgsaucL95 = ci.auc(x)[1]
  res$pgsaucU95 = ci.auc(x)[3]
  
  res2 = R2ROC::auc_diff(combined[, c("trait", "basalpgs", "basal")],1,2,nv=nrow(combined),kv=mean(combined$trait))
  res = cbind(res, data.frame(res2))
  
  return(res)
}

### main function
for (t in c("C10AA")) { # "C10AA", "C10AB", "C10AX09" "C02"
  print(t)
  ttest = t #"C10AB"
  results = data.table::data.table()
  # for (glmnets in c(seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), seq(0.1, 1, 0.1))) {
    for (glmnets in 0.006) {
    print(glmnets)
    res=
      lapply(1:20,
             function(stest) {
               print(stest)
               
               x = splittrainvalidtest(nrow(pheno), stest=stest)
               datatest  = data[x$test, ]
               phenotest  = pheno[x$test, ]
               covariatestest  = covariates[x$test, -c(4:5)]
               
               res=
                 lapply(1:20,
                        function (strain) {
                          print(strain)
                          
                          x = splittrainvalidtest(nrow(pheno), stest=stest, strain=strain)
                          datatrain = data[x$train, ]
                          datavalid = data[x$valid, ]
                          phenotrain = pheno[x$train, ]
                          phenovalid = pheno[x$valid, ]
                          # covariatestrain = covariates[x, ]
                          # covariatestest  = covariates[y, ]
                          # covariatestrain = covariates[x, -c(1:5)]
                          # covariatestest  = covariates[y, -c(1:5)]
                          covariatestrain = covariates[x$train, -c(4:5)] # with PC *1
                          covariatesvalid = covariates[x$valid, -c(4:5)] # with PC *1
                          # covariatestrain = covariates[x$train, 1:3] # wo PC (similar to *1)
                          # covariatesvalid = covariates[x$valid, 1:3]
                          # datatrain = cbind(datatrain, covariates[x, 1:3]) # sex age bmi SNP
                          # datatest  = cbind(datatest,  covariates[y, 1:3])
                          # datatrain = covariates[x, 1:3] # sex age bmi only
                          # datatest  = covariates[y, 1:3]
                          
                          # phenotrain$trait = phenotrain[, t]
                          
                          x = glm2steps(phenotrain[, t], covariatestrain, datatrain, glmnets=glmnets)
                          
                          
                          res = data.table::data.table(trait = t, glmnets=glmnets, stest = stest, strain = strain,
                                                       coef1=list(x$coef1[1:4]), coef2=list(x$coef2))
                          
                          r = assess2steps(phenovalid[, t], covariatesvalid[, c(1:3)], datavalid, x$coef1[1:4], x$coef2)
                          res = cbind(res, r)
                          # r = assess2steps(phenovalid[, t], covariatesvalid[, c(1:3)], datavalid, x$coef1[1:4], x$coef2b)
                          # res = cbind(res, r)
                          return(res)
                        })
               res = do.call(rbind, res)
               coef1 = colMeans(do.call(rbind, res$coef1))
               # coef1 = res$coef1[[1]]
               coef2 = res$coef2
               # for (i in 1:nrow(res)) {
               #   coef2[[i]] = coef2[[i]] * res$pgscoeff[i]
               # }
               coef2 = colMeans(do.call(rbind, coef2))
               if (all(is.na(res$pgscoeff))) {
                 coef2b = coef2 #zero
               } else {
                 # nrow(res) onwards compensates for the zero rows in res$coef2
                 coef2b = coef2 * mean(res$pgscoeff, na.rm=TRUE) * nrow(res) / (nrow(res) - sum(is.na(res$pgscoeff)))
               }
               
               x = assess2steps(phenotest[, ttest], covariatestest[, c(1:3)], datatest, coef1, coef2)
               colnames(x) = paste0(colnames(x), "avg")
               res=cbind(res, x)
               
               x = assess2steps(phenotest[, ttest], covariatestest[, c(1:3)], datatest, coef1, coef2b)
               colnames(x) = paste0(colnames(x), "adj")
               res=cbind(res, x)
               
               return(res)
             })
    res=do.call(rbind, res)
    results = rbind(results, res)
  }
  if (t == ttest) {
    saveRDS(results, file=paste0(t, ".rds"))
  } else {
    saveRDS(results, file=paste0(t, ".test", ttest, ".rds"))
  }
}
