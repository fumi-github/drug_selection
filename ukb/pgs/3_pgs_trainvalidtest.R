library(glmnet)
library(ridge)
library(ggplot2)
library(dplyr)
library(data.table)
library(pROC)
library(R2ROC)

# HELPER FUNCTIONS
# ---
# A function to Winsorize a numeric vector, capping extreme values at the 1st and 99th percentiles.
winsorize = function (x) {
  xmax = quantile(x, 0.99, na.rm = TRUE)
  xmin = quantile(x, 0.01, na.rm = TRUE)
  x[!is.na(x) & x > xmax] = xmax
  x[!is.na(x) & x < xmin] = xmin
  return(x)
}

# A function to split data into training, validation, and test sets.
splittrainvalidtest = function(n, ptest=0.2, stest=1, ptrain=0.6, strain=1) {
  # Set seed for reproducible test set split.
  set.seed(stest)
  shuffled_indices = sample(n)
  
  # Assign indices to the test set.
  test_indices = sort(head(shuffled_indices, round(n*ptest)))
  trainvalid_indices = tail(shuffled_indices, n - round(n*ptest))
  
  # Set seed for reproducible train/validation split.
  set.seed(strain)
  trainvalid_indices = sample(trainvalid_indices)
  
  # Assign indices to the training and validation sets.
  train_indices = sort(head(trainvalid_indices, round(n*ptrain)))
  valid_indices = sort(tail(trainvalid_indices, length(trainvalid_indices) - round(n*ptrain)))
  
  return(list(train=train_indices, valid=valid_indices, test=test_indices))
}

# A function to train a two-step model: 1) a base model with covariates, 2) a penalized model on the residuals.
train_two_step_model = function(outcome_variable, covariate_df, genotype_df, glmnets=0.04) {
  # Step 1: Fit a binomial GLM with covariates to get the base model.
  base_model = glm(outcome_variable ~ .,
                   data=covariate_df,
                   family=binomial(),
                   na.action=na.exclude)
  
  # Extract residuals from the base model, handling NAs appropriately.
  residuals_from_base_model = naresid(base_model$na.action, base_model$residual)

  # The user has the option to winsorize residuals, currently commented out.
  # residuals_from_base_model = winsorize(residuals_from_base_model)
  
  # Remove subjects with missing genotype data or residuals.
  complete_cases_mask = (rowSums(is.na(genotype_df)) == 0 & !is.na(residuals_from_base_model))
  table(complete_cases_mask)
  
  non_na_residuals = residuals_from_base_model[complete_cases_mask]
  non_na_genotypes = genotype_df[complete_cases_mask, ]
  
  # Mean-center the genotype data.
  for (i in 1:ncol(non_na_genotypes)) { 
    non_na_genotypes[, i] = non_na_genotypes[, i] - mean(non_na_genotypes[, i])
  }
  
  # Step 2: Fit a penalized (LASSO/elastic net) model using genotypes to predict the residuals.
  pgs_model_fit = glmnet(
    x = as.matrix(non_na_genotypes),
    y = non_na_residuals)
  
  # Extract SNP coefficients for the specified lambda (penalty parameter).
  snp_coefficients = as.numeric(coef(pgs_model_fit, s=glmnets))[-1]
  
  # Return coefficients from both steps.
  return(list(base_model_coeffs=coef(base_model), pgs_coeffs=snp_coefficients))
}

# A function to cut a continuous variable into deciles.
decilecut = function(score) {
  # Calculate decile quantiles.
  quantiles = quantile(score, seq(0, 1, 0.1), na.rm=TRUE)
  
  # If there are not enough unique quantile values (due to ties), add a small amount of noise to break ties.
  if (length(unique(quantiles)) < 11) {
    min_diff = min(diff(sort(unique(score)))) / 2
    set.seed(1)
    score_jittered = score + runif(length(score), max=min_diff)
    quantiles = quantile(score_jittered, seq(0, 1, 0.1), na.rm=TRUE)
  }
  return(cut(score, quantiles, include.lowest=TRUE))
}

# A function to calculate the proportion of cases per decile of a score.
decileproportion = function(score, outcome) {
  proportions_by_decile = data.frame(decile = decilecut(score), trait = outcome) %>%
    na.omit() %>%
    group_by(decile) %>%
    summarize(proportion = sum(trait=="1") / n(), .groups="drop")
  return(proportions_by_decile$proportion)
}

# A function to assess the performance of the two-step model on a dataset.
evaluate_model_performance = function(outcome_variable, covariate_df, genotype_df, base_coeffs, pgs_coeffs) {
  # Calculate the base model score and the polygenic score (PGS).
  scores_df = data.frame(
    trait = outcome_variable,
    basal_score = as.matrix(cbind(1, covariate_df)) %*% base_coeffs,
    pgs = as.matrix(genotype_df) %*% pgs_coeffs)
  
  # Remove rows with any missing values and calculate the combined score.
  scores_df = scores_df[rowSums(is.na(scores_df)) == 0, ]
  scores_df$combined_score = scores_df$basal_score + scores_df$pgs
  
  # Assess PGS alone.
  pgs_only_model_summary = summary(pgs_only_model <- glm(factor(trait) ~ pgs, data=scores_df, family=binomial()))
  pgs_Z = ifelse(nrow(pgs_only_model_summary$coefficients)>=2, pgs_only_model_summary$coefficients[2,3], NA)
  pgs_R2 = anova(pgs_only_model)[2,2] / anova(pgs_only_model)[1,4]
  
  # Assess base model score alone.
  base_only_model_summary = summary(base_only_model <- glm(factor(trait) ~ basal_score, data=scores_df, family=binomial()))
  basal_Z = ifelse(nrow(base_only_model_summary$coefficients)>=2, base_only_model_summary$coefficients[2,3], NA)
  basal_R2 = anova(base_only_model)[2,2] / anova(base_only_model)[1,4]
  
  # Assess combined score.
  combined_model_summary = summary(combined_model <- glm(factor(trait) ~ combined_score, data=scores_df, family=binomial()))
  combined_Z = ifelse(nrow(combined_model_summary$coefficients)>=2, combined_model_summary$coefficients[2,3], NA)
  combined_R2 = anova(combined_model)[2,2] / anova(combined_model)[1,4]
  
  # Assess base model and PGS as separate predictors.
  full_model_summary = summary(full_model <- glm(factor(trait) ~ basal_score + pgs, data=scores_df, family=binomial()))
  basal_coeff = ifelse(nrow(full_model_summary$coefficients)>=2, full_model_summary$coefficients[2,1], NA)
  pgs_coeff = ifelse(nrow(full_model_summary$coefficients)>=3, full_model_summary$coefficients[3,1], NA)
  
  # Encompassing test: Check if PGS adds information beyond the base model score.
  encompassing_model_summary = summary(encompassing_model <- glm(factor(trait) ~ basal_score + combined_score, data=scores_df, family=binomial()))
  encompass_combined_Z = ifelse(nrow(encompassing_model_summary$coefficients)>=3, encompassing_model_summary$coefficients[3,3], NA)
  
  # Compile results into a data.table.
  results_table = data.table(
    pgsZ=pgs_Z, pgsR2=pgs_R2,
    basalZ=basal_Z, basalR2=basal_R2,
    basalpgsZ=combined_Z, basalpgsR2=combined_R2,
    basalcoeff=basal_coeff, pgscoeff=pgs_coeff, encompassbasalpgsZ=encompass_combined_Z,
    basaldecileproportion    = list(decileproportion(scores_df$basal_score, scores_df$trait)),
    pgsdecileproportion      = list(decileproportion(scores_df$pgs, scores_df$trait)),
    basalpgsdecileproportion = list(decileproportion(scores_df$combined_score, scores_df$trait)))
  
  # Compare AUC of combined model vs. base model using pROC.
  roc_test_results = pROC::roc.test(
    response=scores_df$trait, 
    predictor1=as.numeric(scores_df$combined_score),
    predictor2=as.numeric(scores_df$basal_score))
  roc_test_df = c(roc_test_results$estimate, roc_test_results$statistic, roc_test_results$p.value)
  names(roc_test_df)[4] = "P"
  roc_test_df = data.frame(as.list(roc_test_df))
  results_table = cbind(results_table, roc_test_df)
  
  # Calculate AUC and CI for the PGS alone.
  pgs_roc = pROC::roc(response=scores_df$trait, predictor=as.numeric(scores_df$pgs))
  results_table$pgsauc = pgs_roc$auc
  results_table$pgsaucL95 = ci.auc(pgs_roc)[1]
  results_table$pgsaucU95 = ci.auc(pgs_roc)[3]
  
  # Compare AUC using R2ROC package.
  auc_diff_results = R2ROC::auc_diff(scores_df[, c("trait", "combined_score", "basal_score")],1,2,nv=nrow(scores_df),kv=mean(scores_df$trait))
  results_table = cbind(results_table, data.frame(auc_diff_results))
  
  return(results_table)
}


# DATA LOADING AND PREPARATION
# ---
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

# Load phenotype and covariate data.
phenotype_data = read.table(
# "../gwas/regenie_lipidaemiadrug/ukb_lipidaemiadrug_BT.txt", #20240723
# "../gwas/regenie_lipidaemiadrug_HDLlt40LDLge190/ukb_lipidaemiadrug_HDLlt40LDLge190_BT.txt", #20240819
  "../gwas/regenie_lipidaemiadrug_exclsupplement/ukb_lipidaemiadrug_exclsupplement_BT.txt", #20240827
# "../gwas/regenie_hypertensiondrug_DBPge90SBPge140/ukb_hypertensiondrug_DBPge90SBPge140_BT.txt", #20241018
# "../gwas/regenie_diabetesdrug_HbA1cge48/ukb_diabetesdrug_HbA1cge48_BT.txt", #20241018
    header=TRUE, sep=" ")
covariate_data = read.table(
# "../gwas/regenie_lipidaemiadrug/ukb_lipidaemiadrug_covariates.txt", #20240723
# "../gwas/regenie_lipidaemiadrug_HDLlt40LDLge190/ukb_lipidaemiadrug_HDLlt40LDLge190_covariates.txt", #20240819
  "../gwas/regenie_lipidaemiadrug_exclsupplement/ukb_lipidaemiadrug_exclsupplement_covariates.txt", #20240827
# "../gwas/regenie_hypertensiondrug_DBPge90SBPge140/ukb_hypertensiondrug_DBPge90SBPge140_covariates.txt", #20241018
# "../gwas/regenie_diabetesdrug_HbA1cge48/ukb_diabetesdrug_HbA1cge48_covariates.txt", #20241018
    header=TRUE, sep=" ")

# Align data frames and remove sample ID columns (FID, IID).
identical(phenotype_data$FID, covariate_data$FID)
genotype_data = genotype_data[match(phenotype_data$FID, genotype_data$FID), ]
genotype_data = genotype_data[, -c(1:2)]
phenotype_data = phenotype_data[, -c(1:2)]
covariate_data = covariate_data[, -c(1:2)]


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
               
               x = splittrainvalidtest(nrow(phenotype_data), stest=stest)
               datatest  = genotype_data[x$test, ]
               phenotest  = phenotype_data[x$test, ]
               covariatestest  = covariate_data[x$test, -c(4:5)]
               
               res=
                 lapply(1:20,
                        function (strain) {
                          print(strain)
                          
                          x = splittrainvalidtest(nrow(phenotype_data), stest=stest, strain=strain)
                          datatrain = genotype_data[x$train, ]
                          datavalid = genotype_data[x$valid, ]
                          phenotrain = phenotype_data[x$train, ]
                          phenovalid = phenotype_data[x$valid, ]
                          # covariatestrain = covariate_data[x, ]
                          # covariatestest  = covariate_data[y, ]
                          # covariatestrain = covariate_data[x, -c(1:5)]
                          # covariatestest  = covariate_data[y, -c(1:5)]
                          covariatestrain = covariate_data[x$train, -c(4:5)] # with PC *1
                          covariatesvalid = covariate_data[x$valid, -c(4:5)] # with PC *1
                          # covariatestrain = covariate_data[x$train, 1:3] # wo PC (similar to *1)
                          # covariatesvalid = covariate_data[x$valid, 1:3]
                          # datatrain = cbind(datatrain, covariate_data[x, 1:3]) # sex age bmi SNP
                          # datatest  = cbind(datatest,  covariate_data[y, 1:3])
                          # datatrain = covariate_data[x, 1:3] # sex age bmi only
                          # datatest  = covariate_data[y, 1:3]
                          
                          # phenotrain$trait = phenotrain[, t]
                          
                          x = train_two_step_model(phenotrain[, t], covariatestrain, datatrain, glmnets=glmnets)
                          
                          
                          res = data.table::data.table(trait = t, glmnets=glmnets, stest = stest, strain = strain,
                                                       coef1=list(x$base_model_coeffs[1:4]), coef2=list(x$pgs_coeffs))
                          
                          r = evaluate_model_performance(phenovalid[, t], covariatesvalid[, c(1:3)], datavalid, x$base_model_coeffs[1:4], x$pgs_coeffs)
                          res = cbind(res, r)
                          # r = evaluate_model_performance(phenovalid[, t], covariatesvalid[, c(1:3)], datavalid, x$base_model_coeffs[1:4], x$pgs_coeffsb)
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
               
               x = evaluate_model_performance(phenotest[, ttest], covariatestest[, c(1:3)], datatest, coef1, coef2)
               colnames(x) = paste0(colnames(x), "avg")
               res=cbind(res, x)
               
               x = evaluate_model_performance(phenotest[, ttest], covariatestest[, c(1:3)], datatest, coef1, coef2b)
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
