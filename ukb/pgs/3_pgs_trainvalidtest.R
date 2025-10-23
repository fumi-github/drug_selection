library(glmnet)
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
  non_na_genotypes <- scale(non_na_genotypes, center=TRUE, scale=FALSE)
  
  
  # Step 2: Fit a penalized (LASSO) model using genotypes to predict the residuals.
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
    score = score + runif(length(score), max=min_diff)
    quantiles = quantile(score, seq(0, 1, 0.1), na.rm=TRUE)
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
  pgs_roc = pROC::roc(response=scores_df$trait, predictor=as.numeric(scores_df$pgs),
                      direction="<")
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
  # "../gwas/regenie_lipidaemiadrug_exclsupplement/ukb_lipidaemiadrug_exclsupplement_BT.txt", #20240827
  "../gwas/regenie_lipidaemiadrug_onlyaffectedexclsupplement/ukb_lipidaemiadrug_onlyaffectedexclsupplement_BT.txt", #20250928
  # "../gwas/regenie_hypertensiondrug_DBPge90SBPge140/ukb_hypertensiondrug_DBPge90SBPge140_BT.txt", #20241018
  # "../gwas/regenie_hypertensiondrug_onlyaffectedDBPge100SBPge160/ukb_hypertensiondrug_onlyaffectedDBPge100SBPge160_BT.txt", #20250928
  # "../gwas/regenie_diabetesdrug_HbA1cge48/ukb_diabetesdrug_HbA1cge48_BT.txt", #20241018
  header=TRUE, sep=" ")
covariate_data = read.table(
  # "../gwas/regenie_lipidaemiadrug/ukb_lipidaemiadrug_covariates.txt", #20240723
  # "../gwas/regenie_lipidaemiadrug_HDLlt40LDLge190/ukb_lipidaemiadrug_HDLlt40LDLge190_covariates.txt", #20240819
  # "../gwas/regenie_lipidaemiadrug_exclsupplement/ukb_lipidaemiadrug_exclsupplement_covariates.txt", #20240827
  "../gwas/regenie_lipidaemiadrug_onlyaffectedexclsupplement/ukb_lipidaemiadrug_onlyaffectedexclsupplement_covariates.txt", #20250928
  # "../gwas/regenie_hypertensiondrug_DBPge90SBPge140/ukb_hypertensiondrug_DBPge90SBPge140_covariates.txt", #20241018
  # "../gwas/regenie_hypertensiondrug_onlyaffectedDBPge100SBPge160/ukb_hypertensiondrug_onlyaffectedDBPge100SBPge160_covariates.txt", #20250928
  # "../gwas/regenie_diabetesdrug_HbA1cge48/ukb_diabetesdrug_HbA1cge48_covariates.txt", #20241018
  header=TRUE, sep=" ")

# Align data frames and remove sample ID columns (FID, IID).
identical(phenotype_data$FID, covariate_data$FID)
genotype_data = genotype_data[match(phenotype_data$FID, genotype_data$FID), ]
genotype_data = genotype_data[, -c(1:2)]
phenotype_data = phenotype_data[, -c(1:2)]
covariate_data = covariate_data[, -c(1:2)]


# MAIN ANALYSIS
# ---
# Define the trait(s) to be analyzed and the specific trait to be used for testing.
for (current_trait in c("C10AA")) { # Example traits: "C10AA", "C10AB", "C10AX09", "C03"
  print(current_trait)
  for (test_trait in current_trait) {
  # for (test_trait in setdiff(c("C10AA", "C10AB", "C10AX09"), current_trait)) {
  # for (test_trait in setdiff(c("C03", "C07", "C08", "C09"), current_trait)) {
    all_results = data.table::data.table()
  
  # Loop over the glmnet penalty parameter (lambda). Currently set to a single value.
  # for (glmnets_param in c(seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), seq(0.1, 1, 0.1))) {
  for (glmnets_param in 0.006) {
    print(glmnets_param)
    
    # Outer loop for repeated hold-out validation, creating 20 different test sets.
    test_run_results_list =
      lapply(1:20,
             function(test_seed) {
               print(test_seed)
               
               # Split data based on the current test seed.
               split_indices = splittrainvalidtest(nrow(phenotype_data), stest=test_seed)
               genotype_test  = genotype_data[split_indices$test, ]
               phenotype_test  = phenotype_data[split_indices$test, ]
               covariate_test  = covariate_data[split_indices$test, -c(4:5)] # with PC *1
               
               # Inner loop to create 20 different train/validation splits for each test set.
               # This stabilizes model coefficient estimation.
               validation_run_results_list =
                 lapply(1:20,
                        function (train_seed) {
                          print(train_seed)
                          
                          # Split data into train and validation sets.
                          inner_split_indices = splittrainvalidtest(nrow(phenotype_data), stest=test_seed, strain=train_seed)
                          genotype_train = genotype_data[inner_split_indices$train, ]
                          genotype_valid = genotype_data[inner_split_indices$valid, ]
                          phenotype_train = phenotype_data[inner_split_indices$train, ]
                          phenotype_valid = phenotype_data[inner_split_indices$valid, ]
                          covariate_train = covariate_data[inner_split_indices$train, -c(4:5)] # with PC *1
                          covariate_valid = covariate_data[inner_split_indices$valid, -c(4:5)] # with PC *1
                          
                          # Train the two-step model on the training data.
                          trained_model = train_two_step_model(phenotype_train[, current_trait], covariate_train, genotype_train, glmnets=glmnets_param)
                          
                          # Store the results and coefficients.
                          single_run_result = data.table::data.table(trait = current_trait, glmnets=glmnets_param, stest = test_seed, strain = train_seed,
                                                                     coef1=list(trained_model$base_model_coeffs[1:4]), coef2=list(trained_model$pgs_coeffs))
                          
                          # Evaluate the trained model on the validation set.
                          validation_assessment = evaluate_model_performance(phenotype_valid[, current_trait], covariate_valid[, c(1:3)], genotype_valid, trained_model$base_model_coeffs[1:4], trained_model$pgs_coeffs)
                          
                          # Combine results.
                          single_run_result = cbind(single_run_result, validation_assessment)
                          return(single_run_result)
                        })
               
               # Combine results from all 20 validation runs.
               validation_results_df = do.call(rbind, validation_run_results_list)
               
               # Average the coefficients from the 20 validation runs to create a stable final model.
               avg_base_coeffs = colMeans(do.call(rbind, validation_results_df$coef1))
               avg_pgs_coeffs = colMeans(do.call(rbind, validation_results_df$coef2))
               
               # Create an adjusted version of the PGS coefficients, rescaled by their average effect size in the validation models.
               if (all(is.na(validation_results_df$pgscoeff))) {
                 adjusted_avg_pgs_coeffs = avg_pgs_coeffs # No adjustment if all are NA; actually is zero
               } else {
                 rescaling_factor = mean(validation_results_df$pgscoeff, na.rm=TRUE)
                 # Compensates for the zero rows in validation_results_df$coef2
                 rescaling_factor = rescaling_factor * nrow(validation_results_df) / (nrow(validation_results_df) - sum(is.na(validation_results_df$pgscoeff)))
                 adjusted_avg_pgs_coeffs = avg_pgs_coeffs * rescaling_factor
               }
               
               # Evaluate the averaged model on the held-out test set.
               test_assessment_avg = evaluate_model_performance(phenotype_test[, test_trait], covariate_test[, c(1:3)], genotype_test, avg_base_coeffs, avg_pgs_coeffs)
               colnames(test_assessment_avg) = paste0(colnames(test_assessment_avg), "avg")
               
               # Evaluate the adjusted model on the held-out test set.
               test_assessment_adj = evaluate_model_performance(phenotype_test[, test_trait], covariate_test[, c(1:3)], genotype_test, avg_base_coeffs, adjusted_avg_pgs_coeffs)
               colnames(test_assessment_adj) = paste0(colnames(test_assessment_adj), "adj")
               
               # Combine the validation results with the test assessments.
               final_results_for_test_seed = cbind(validation_results_df, test_assessment_avg, test_assessment_adj)
               
               return(final_results_for_test_seed)
             })
    
    # Combine results from all 20 test runs.
    combined_test_runs_df = do.call(rbind, test_run_results_list)
    all_results = rbind(all_results, combined_test_runs_df)
  }
  
  # Save the final results to an RDS file.
  if (current_trait == test_trait) {
    saveRDS(all_results, file=paste0(current_trait, ".rds"))
  } else {
    saveRDS(all_results, file=paste0(current_trait, ".test", test_trait, ".rds"))
  }
}
}
