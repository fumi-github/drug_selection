#!/bin/bash

t=C10AA
# f1=../../ukb/gwas/regenie_hypolipidemics.LDL.minafterstart28.wpower1/ukb_hypolipidemics.LDL_step2_QT_chrall_statin.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDL/${t}/aou_hypolipidemics.LDL_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox_step2_QT_chrall_statin.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcox/${t}/aou_hypolipidemics.LDLmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox01.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox01_step2_QT_chrall_statin.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcox01/${t}/aou_hypolipidemics.LDLmyboxcox01_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLdelta.minafterstart28.wpower1/ukb_hypolipidemics.LDLdelta_step2_QT_chrall_statin.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDLdelta/${t}/aou_hypolipidemics.LDLdelta_step2_out_chrall_${t}.regenie.rs.geno01maf001mac0hwedev02
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLlogratio.minafterstart28.wpower1/ukb_hypolipidemics.LDLlogratio_step2_QT_chrall_statin.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDLlogratio/${t}/aou_hypolipidemics.LDLlogratio_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox01.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox01_step2_QT_chrall_statin.regenie.rs.100x
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcox01/${t}/aou_hypolipidemics.LDLmyboxcox01_step2_out_chrall_${t}.regenie.rs.100x
f1=../../ukb/gwas/regenie_hypolipidemics.LDLlogratio.minafterstart28.wpower1/ukb_hypolipidemics.LDLlogratio_step2_QT_chrall_statin.regenie.rs.100x
f2=../../aou/regenie_hypolipidemics.LDLlogratio/${t}/aou_hypolipidemics.LDLlogratio_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox01.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox01_step2_QT_chrall_statin.regenie.rs.P
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcox01/${t}/aou_hypolipidemics.LDLmyboxcox01_step2_out_chrall_${t}.regenie.rs.P
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLdelta.minafterstart28.wpower1/ukb_hypolipidemics.LDLdelta_step2_QT_chrall_statin.regenie.rs.P
# f2=../../aou/regenie_hypolipidemics.LDLdelta/${t}/aou_hypolipidemics.LDLdelta_step2_out_chrall_${t}.regenie.rs.P
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLlogratio.minafterstart28.wpower1/ukb_hypolipidemics.LDLlogratio_step2_QT_chrall_statin.regenie.rs.P
# f2=../../aou/regenie_hypolipidemics.LDLlogratio/${t}/aou_hypolipidemics.LDLlogratio_step2_out_chrall_${t}.regenie.rs.P
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcoxv2.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcoxv2_step2_QT_chrall_statin.regenie.rs.100x
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcoxv2/${t}/aou_hypolipidemics.LDLmyboxcoxv2_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02

# t=C10AX09
# f1=../../ukb/gwas/regenie_hypolipidemics.LDL.minafterstart28.wpower1/ukb_hypolipidemics.LDL_step2_QT_chrall_ezetimibe.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDL/${t}/aou_hypolipidemics.LDL_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox_step2_QT_chrall_ezetimibe.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcox/${t}/aou_hypolipidemics.LDLmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcox01.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcox01_step2_QT_chrall_ezetimibe.regenie.rs
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcox01/${t}/aou_hypolipidemics.LDLmyboxcox01_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_hypolipidemics.LDLmyboxcoxv2.minafterstart28.wpower1/ukb_hypolipidemics.LDLmyboxcoxv2_step2_QT_chrall_ezetimibe.regenie.rs.100x
# f2=../../aou/regenie_hypolipidemics.LDLmyboxcoxv2/${t}/aou_hypolipidemics.LDLmyboxcoxv2_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02

# t=A10BA
# f1=../../ukb/gwas/regenie_antidiabetic.HbA1c.wpower1/ukb_antidiabetic.HbA1c_step2_QT_chrall_biguanide.regenie.rs
# f2=../../aou/regenie_antidiabetic.HbA1c/${t}/aou_antidiabetic.HbA1c_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antidiabetic.HbA1cmyboxcox.wpower1/ukb_antidiabetic.HbA1cmyboxcox_step2_QT_chrall_biguanide.regenie.rs
# f2=../../aou/regenie_antidiabetic.HbA1cmyboxcox/${t}/aou_antidiabetic.HbA1cmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antidiabetic.HbA1cmyboxcox01.wpower1/ukb_antidiabetic.HbA1cmyboxcox01_step2_QT_chrall_biguanide.regenie.rs
# f2=../../aou/regenie_antidiabetic.HbA1cmyboxcox01/${t}/aou_antidiabetic.HbA1cmyboxcox01_step2_out_chrall_${t}.regenie.rs

# t=A10BB
# f1=../../ukb/gwas/regenie_antidiabetic.HbA1c.wpower1/ukb_antidiabetic.HbA1c_step2_QT_chrall_sulfonylurea.regenie.rs
# f2=../../aou/regenie_antidiabetic.HbA1c/${t}/aou_antidiabetic.HbA1c_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antidiabetic.HbA1cmyboxcox.wpower1/ukb_antidiabetic.HbA1cmyboxcox_step2_QT_chrall_sulfonylurea.regenie.rs
# f2=../../aou/regenie_antidiabetic.HbA1cmyboxcox/${t}/aou_antidiabetic.HbA1cmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antidiabetic.HbA1cmyboxcox01.wpower1/ukb_antidiabetic.HbA1cmyboxcox01_step2_QT_chrall_sulfonylurea.regenie.rs
# f2=../../aou/regenie_antidiabetic.HbA1cmyboxcox01/${t}/aou_antidiabetic.HbA1cmyboxcox01_step2_out_chrall_${t}.regenie.rs

# t=C03
# f1=../../ukb/gwas/regenie_antihypertensive.SBP.minafterstart28.wpower1/ukb_antihypertensive.SBP_step2_QT_chrall_thiazide.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBP/${t}/aou_antihypertensive.SBP_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox_step2_QT_chrall_thiazide.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox/${t}/aou_antihypertensive.SBPmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox01.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox01_step2_QT_chrall_thiazide.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox01/${t}/aou_antihypertensive.SBPmyboxcox01_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcoxv2.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcoxv2_step2_QT_chrall_thiazide.regenie.rs.100x
# f2=../../aou/regenie_antihypertensive.SBPmyboxcoxv2/${t}/aou_antihypertensive.SBPmyboxcoxv2_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02

# t=C07
# f1=../../ukb/gwas/regenie_antihypertensive.SBP.minafterstart28.wpower1/ukb_antihypertensive.SBP_step2_QT_chrall_b-blocker.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBP/${t}/aou_antihypertensive.SBP_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox_step2_QT_chrall_b-blocker.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox/${t}/aou_antihypertensive.SBPmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox01.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox01_step2_QT_chrall_b-blocker.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox01/${t}/aou_antihypertensive.SBPmyboxcox01_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcoxv2.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcoxv2_step2_QT_chrall_b-blocker.regenie.rs.100x
# f2=../../aou/regenie_antihypertensive.SBPmyboxcoxv2/${t}/aou_antihypertensive.SBPmyboxcoxv2_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02

# t=C08
# f1=../../ukb/gwas/regenie_antihypertensive.SBP.minafterstart28.wpower1/ukb_antihypertensive.SBP_step2_QT_chrall_CCB.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBP/${t}/aou_antihypertensive.SBP_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox_step2_QT_chrall_CCB.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox/${t}/aou_antihypertensive.SBPmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox01.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox01_step2_QT_chrall_CCB.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox01/${t}/aou_antihypertensive.SBPmyboxcox01_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcoxv2.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcoxv2_step2_QT_chrall_CCB.regenie.rs.100x
# f2=../../aou/regenie_antihypertensive.SBPmyboxcoxv2/${t}/aou_antihypertensive.SBPmyboxcoxv2_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02

# t=C09
# f1=../../ukb/gwas/regenie_antihypertensive.SBP.minafterstart28.wpower1/ukb_antihypertensive.SBP_step2_QT_chrall_ACEi_ARB.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBP/${t}/aou_antihypertensive.SBP_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox_step2_QT_chrall_ACEi_ARB.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox/${t}/aou_antihypertensive.SBPmyboxcox_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcox01.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcox01_step2_QT_chrall_ACEi_ARB.regenie.rs
# f2=../../aou/regenie_antihypertensive.SBPmyboxcox01/${t}/aou_antihypertensive.SBPmyboxcox01_step2_out_chrall_${t}.regenie.rs
# f1=../../ukb/gwas/regenie_antihypertensive.SBPmyboxcoxv2.minafterstart28.wpower1/ukb_antihypertensive.SBPmyboxcoxv2_step2_QT_chrall_ACEi_ARB.regenie.rs.100x
# f2=../../aou/regenie_antihypertensive.SBPmyboxcoxv2/${t}/aou_antihypertensive.SBPmyboxcoxv2_step2_out_chrall_${t}.regenie.rs.100x.geno01maf001mac0hwedev02

# t=SBPnodrugmyboxcox_s
# f1=../../ukb/gwas/regenie_SBPnodrugpartialmyboxcox-04nge4/ukb_SBPnodrugpartialmyboxcox-04nge4_step2_QT_chrall_s.regenie.rs
# f2=../../aou/regenie_SBPnodrugmyboxcox-05nge4/s/aou_SBPnodrugmyboxcox-05nge4_step2_out_chrall_s.regenie.rs

# t=MAPnodrugmyboxcox_s.100x
# f1=../../ukb/gwas/regenie_MAPnodrugpartialmyboxcox-02nge4/ukb_MAPnodrugpartialmyboxcox-02nge4_step2_QT_chrall_s.regenie.rs.100x
# f2=../../aou/regenie_MAPnodrugmyboxcox0nge4/s/aou_MAPnodrugmyboxcox0nge4_step2_out_chrall_s.regenie.rs.100x

# t=MAPinstance0myboxcox_s.100x
# f1=../../ukb/gwas/regenie_MAPinstance0myboxcox01nge2/ukb_MAPinstance0myboxcox01nge2_step2_QT_chrall_s.regenie.rs.100x
# f2=../../ofh/plink_MAPclinicmyboxcox-01nge2/ofh_snv.v8.chrall-ball.s.glm.linear.ADD.rs.100x

# outfile=${t}
# outfile=${t}.myboxcox
# outfile=${t}.myboxcox01
# outfile=${t}.delta
outfile=${t}.logratio
# outfile=${t}.myboxcox01.100x
#outfile=${t}.myboxcoxv2.100x
# outfile=${t}.logratio.100x
# outfile=${t}.myboxcox01.SAMPLESIZE
# outfile=${t}.delta.SAMPLESIZE
# outfile=${t}.logratio.SAMPLESIZE

~/final/opt/metal << EOF > ${outfile}.log

SCHEME STDERR
#SCHEME SAMPLESIZE
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize

#VERBOSE ON

###
MARKER rsAminAmax
WEIGHT N
LABEL TotalSampleSize as N
ALLELE ALLELE1 ALLELE0
FREQ A1FREQ
EFFECT BETA
STDERR SE
# PVALUE P

REMOVEFILTERS
ADDFILTER	INFO >= 0.3
ADDFILTER	A1FREQ >= 0.01
ADDFILTER	A1FREQ <= 0.99

PROCESS ${f1}

###
MARKER rsAminAmax
WEIGHT N
LABEL TotalSampleSize as N
ALLELE ALLELE1 ALLELE0
FREQ A1FREQ
EFFECT BETA
STDERR SE
# PVALUE P
# ### PLINK
# MARKER rsAminAmax
# WEIGHT OBS_CT
# LABEL TotalSampleSize as OBS_CT
# ALLELE A1 OMITTED
# FREQ A1_FREQ
# EFFECT BETA
# STDERR SE

REMOVEFILTERS
ADDFILTER	INFO >= 0.3
ADDFILTER	A1FREQ >= 0.01
ADDFILTER	A1FREQ <= 0.99
# ### PLINK
# REMOVEFILTERS
# ADDFILTER	A1_FREQ >= 0.01
# ADDFILTER	A1_FREQ <= 0.99

PROCESS ${f2}

OutFile ${outfile}_ .tbl
ANALYZE HETEROGENEITY

MINWEIGHT 200

EOF
