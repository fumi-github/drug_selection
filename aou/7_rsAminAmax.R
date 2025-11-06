dbSNP = as.data.frame(data.table::fread("hg38_dbSNP155Common.tsv.gz", sep="\t", quote=""))
dbSNP$chrpos = paste0(dbSNP$`#chrom`, ":", dbSNP$chromEnd)

files = dir(
  # "regenie_hypolipidemics.LDL",
  # "regenie_hypolipidemics.LDLmyboxcox",
  # "regenie_hypolipidemics.LDLmyboxcox01",
  # "regenie_hypolipidemics.LDLmyboxcoxv2",
  # "regenie_hypolipidemics.LDLdelta",
  "regenie_hypolipidemics.LDLlogratio/C10AA",
  # "regenie_antihypertensive.SBP",
  # "regenie_antihypertensive.SBPmyboxcox",
  # "regenie_antihypertensive.SBPmyboxcox01",
  # "regenie_antihypertensive.SBPmyboxcoxv2",
  # "regenie_antidiabetic.HbA1c",
  # "regenie_antidiabetic.HbA1c/A10BB",
  # "regenie_antidiabetic.HbA1cmyboxcox",
  # "regenie_antidiabetic.HbA1cmyboxcox01",
  # "regenie_SBPnodrugmyboxcox-05nge4",
  # "regenie_MAPnodrugmyboxcox0nge4",
  # "regenie_PhysicalMeasurementsv2_MAPmyboxcoxnge2/MAPmyboxcox01s",
  pattern="*.regenie.gz",
  full.names=TRUE,
  recursive=TRUE)
#"regenie_antihypertensive.SBP/C08/aou_antihypertensive.SBP_step2_out_chr22_C08.regenie.gz"
for (f in files) {
  print(f)
  data = as.data.frame(data.table::fread(f, sep=" ", quote=""))
  data$chrpos = paste0("chr", data$CHROM, ":", data$GENPOS)
  data$rs = dbSNP$name[match(data$chrpos, dbSNP$chrpos)]
  
  output = data[!is.na(data$rs), ]
  output$rsAminAmax = paste0(output$rs,
                             ":", pmin(output$ALLELE0, output$ALLELE1),
                             ":", pmax(output$ALLELE0, output$ALLELE1))
  output = output[, ! colnames(output) %in% c("chrpos", "rs")]

  # # HbA1c ONLY !!!
  # # HbA1c [%]: convert to HbA1c [mmol/mol], IFCC = (10.93 x DCCT) - 23.5
  # output$BETA = 10.93 * output$BETA
  # output$SE   = 10.93 * output$SE

  # # LDL ONLY !!! only for raw not log !!!
  # # mg/dL: convert to mmol/L
  # output$BETA = 0.02586 * output$BETA
  # output$SE   = 0.02586 * output$SE
  
  write.table(output, file=sub("regenie.gz", "regenie.rs", f), quote=FALSE, row.names=FALSE)
}

files = dir(
  # "../ukb/gwas/regenie_antidiabetic.HbA1c.wpower1/",
  # "../ukb/gwas/regenie_antidiabetic.HbA1cmyboxcox.wpower1/",
  # "../ukb/gwas/regenie_antidiabetic.HbA1cmyboxcox01.wpower1/",
  # "../ukb/gwas/regenie_antihypertensive.SBP.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_antihypertensive.SBPmyboxcox.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_antihypertensive.SBPmyboxcox01.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_antihypertensive.SBPmyboxcoxv2.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_hypolipidemics.LDL.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox01.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_hypolipidemics.LDLmyboxcoxv2.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_hypolipidemics.LDLdelta.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_hypolipidemics.LDLlogratio.minafterstart28.wpower1/",
  # "../ukb/gwas/regenie_SBPnodrugpartialmyboxcox-04nge4/",
  # "../ukb/gwas/regenie_MAPnodrugpartialmyboxcox-02nge4/",
  "../ukb/gwas/regenie_MAPinstance0myboxcox01nge2/",
  pattern="*.regenie.gz",
  full.names=TRUE,
  recursive=TRUE)
# "../ukb/gwas/regenie_antidiabetic.HbA1c.wpower1/ukb_antidiabetic.HbA1c_step2_QT_chr10_biguanide.regenie.gz"
for (f in files) {
  print(f)
  data = as.data.frame(data.table::fread(f, sep=" ", quote=""))
  
  output = data[grep("^rs", data$ID), ]
  output$rsAminAmax = paste0(output$ID,
                             ":", pmin(output$ALLELE0, output$ALLELE1),
                             ":", pmax(output$ALLELE0, output$ALLELE1))
  write.table(output, file=sub("regenie.gz", "regenie.rs", f), quote=FALSE, row.names=FALSE)
}

files = dir(
  "../ukbaou/METAL/",
  # pattern="*.myboxcox01_1.tbl.P1e-6_HetDf1",
  # pattern="*SAMPLESIZE_1.tbl.P5e-8_HetDf1",
  # pattern="C10AA.delta_1.tbl.P5e-8_HetDf1",
  # pattern="C10AA.logratio_1.tbl.P5e-8_HetDf1",
  # pattern="SBPnodrugmyboxcox_s_1.tbl.P1e-6_HetDf1",
  pattern="MAPinstance0myboxcox_s.100x_1.tbl.P5e-8_HetDf1",
  full.names=TRUE,
  recursive=TRUE)
for (f in files) {
  data = read.table(f, header=TRUE)
  x = match(sub(":.*", "", data$MarkerName), dbSNP$name)
  data$Chr = as.numeric(sub("chr", "", dbSNP$`#chrom`[x]))
  data$Pos = dbSNP$chromEnd[x]
  data = data[order(data$Chr, data$Pos), ]
  write.table(data, file=paste0(f, ".ChrPos"), quote=FALSE, row.names=FALSE)
}

files = dir(
  "../ofh/plink_MAPclinicmyboxcox-01nge2",
  pattern="*.ADD.gz",
  full.names=TRUE,
  recursive=TRUE)
for (f in files) {
  print(f)
  data = as.data.frame(data.table::fread(f))
  data$chrpos = paste0("chr", data$`#CHROM`, ":", data$POS)
  data$rs = dbSNP$name[match(data$chrpos, dbSNP$chrpos)]
  
  output = data[!is.na(data$rs), ]
  output$rsAminAmax = paste0(output$rs,
                             ":", pmin(output$REF, output$ALT),
                             ":", pmax(output$REF, output$ALT))
  output = output[, ! colnames(output) %in% c("chrpos", "rs")]

  write.table(output, file=sub("ADD.gz", "ADD.rs", f), quote=FALSE, row.names=FALSE)
}