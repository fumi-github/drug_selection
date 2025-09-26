# keyrs keyhg38 are unique keys for bi-allele pairs
# N.B.
#  ukb$keyrs is on hg19 positive strand
#  dbsnp$keyrs (and aou$keyrs) is on hg38 positive strand
### Prepare ukb
ukb = 
  do.call(
    rbind,
    lapply(
      1:22,
      function (i) {
        as.data.frame(data.table::fread(
          paste0(
            "/baker/datasets/ukb55469/Genetic_data/Imputation/Pgen/Stringent_filter/ukb_imp_chr",
            i,
            "_v3_qc.pvar"))) }))
colnames(ukb)[1] = "CHROM"
all(! duplicated(ukb$ID))
ukb$keyrs = paste(ukb$ID,
                  pmin(ukb$REF, ukb$ALT),
                  pmax(ukb$REF, ukb$ALT),
                  sep=" ")
all(! duplicated(ukb$keyrs))

### Prepare dbsnp
dbsnp = as.data.frame(data.table::fread("../../aou/hg38_dbSNP155Common.tsv.gz"))
colnames(dbsnp)[1] = "chrom"
dbsnp = dbsnp[dbsnp$chrom %in% paste0("chr", 1:22), ]
dbsnp = tidyr::separate_rows(dbsnp, alts, sep=",")
dbsnp = dbsnp[! dbsnp$alts=="", ]
dbsnp$keyrs = paste(dbsnp$name,
                    pmin(dbsnp$ref, dbsnp$alts),
                    pmax(dbsnp$ref, dbsnp$alts),
                    sep=" ")
dbsnp$keyhg38 = paste(dbsnp$chrom,
                      dbsnp$chromStart + 1,
                      pmin(dbsnp$ref, dbsnp$alts),
                      pmax(dbsnp$ref, dbsnp$alts),
                      sep=" ")
all(! duplicated(dbsnp$keyrs))
all(! duplicated(dbsnp$keyhg38))

### Prepare aou
aou = read.table("../../aou/regenie_lipidaemiadrug/regenie_lipidaemiadrug.C10AA.ID",
                 header=TRUE)
aou = tidyr::separate(
  aou,
  ID,
  into=c("chr", "pos", "A1", "A2"),
  sep=":",
  remove=FALSE)
aou$keyhg38 = paste(aou$chr,
                    aou$pos,
                    pmin(aou$A1, aou$A2),
                    pmax(aou$A1, aou$A2),
                    sep=" ")
all(! duplicated(aou$keyhg38)) #FALSE
# There are wierd duplications
#                 ID  chr   pos A1 A2         keyhg38
# 13 chr1:31719:G:GA chr1 31719  G GA chr1 31719 G GA
# 14 chr1:31719:GA:G chr1 31719 GA  G chr1 31719 G GA
x = (aou$keyhg38 %in% dbsnp$keyhg38)
table(x) #2175024 6414424
aou = aou[x, ]
aou$keyrs = dbsnp$keyrs[match(aou$keyhg38, dbsnp$keyhg38)]
# Only keep overlap with ukb
x = (aou$keyrs %in% ukb$keyrs)
table(x) #291724 6122700
aou = aou[x, ]

### Get rs of genome-wide significant SNPs from literature.
# lipid not UK Biobank, AllofUs
foo = read.table("willer_glgc-lipids2013_157loci.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
foo[, 4] = toupper(foo[, 4])
foo[, 5] = toupper(foo[, 5])
colnames(foo)[2:4] = c("pos", "rs", "EA")
snps = foo[, 1:4]

# # lipid include UK Biobank; not AllofUs
# foo = read.table("glgc-lipids2021_ST5_1750indexvariants.txt", header=TRUE, sep="\t")
# dim(foo)
# colnames(foo)[2:3] = c("pos", "rs")
# snps = foo[, 1:4]

# BP not UK Biobank, AllofUs
foo = read.table("NatGenet_48_1151_30novel.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
foo$pos = NA
colnames(foo)[1:2] = c("chr", "rs")
snps = rbind(snps, foo[, c(1, 5, 2:3)])
#
foo = read.table("NatGenet_48_1151_50discovery.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
foo$EA = NA
colnames(foo)[1:3] = c("rs", "chr", "pos")
snps = rbind(snps, foo[, c(2:3, 1, 4)])
#
foo = read.table("NatGenet_48_1162_39known.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:4] = c("rs", "chr", "pos", "EA")
snps = rbind(snps, foo[, c(2:3, 1, 4)])
#
foo = read.table("NatGenet_48_1162_31novel.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:4] = c("rs", "chr", "pos", "EA")
snps = rbind(snps, foo[, c(2:3, 1, 4)])
#
foo = read.table("NatGenet_49_54_GERAICBPmeta.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:4] = c("rs", "chr", "pos", "EA")
snps = rbind(snps, foo[, c(2:3, 1, 4)])

# # BP include UK Biobank
# foo = read.table("NatGenet_48_1171_stage4.txt", header=TRUE, sep="\t")
# dim(foo); head(foo)
# colnames(foo)[1:4] = c("rs", "chr", "pos", "EA")
# snps = rbind(snps, foo[, c(2:3, 1, 4)])
#
# foo = read.table("GCST90310294_keaton2024_ST10_2103independentsignals.hg19.txt", header=TRUE, sep="\t")
# dim(foo)
# colnames(foo)[1:4] = c("rs", "chr", "pos", "EA")
# snps = rbind(snps, foo[, c(2:3, 1, 4)])

# T2D not UK Biobank, AllofUs
foo = read.table("Scott2017_DIAGRAM_128signals.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:3] = c("chr", "pos", "rs")
snps = rbind(snps, foo[, c(1:4)])

# # T2D include UK Biobank
# foo = read.table("Suzuki2024_ST4_1289independentsignals.txt", header=TRUE, sep="\t")
# dim(foo)
# colnames(foo) = c("chr", "rs", "pos", "EA")
# snps = rbind(snps, foo[, c(1, 3, 2, 4)])

# BMI not UK Biobank, AllofUs; 
foo = read.table("Locke2015_GIANTBMI_97loci.hg18.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:4] = c("rs", "chr", "pos", "EA")
foo$pos = NA #hg18!!!
snps = rbind(snps, foo[, c(2:3, 1, 4)])

# # BMI include UK Biobank; GRCh37
# foo = read.table("Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt", header=TRUE, sep="\t")
# dim(foo)
# colnames(foo) = c("rs", "chr", "pos", "EA")
# snps = rbind(snps, foo[, c(2, 3, 1, 4)])

# CAD not UK Biobank, AllofUs
foo = read.table("CARDIoGRAMplusC4D2015.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:4] = c("chr", "rs", "pos", "EA")
snps = rbind(snps, foo[, c(1, 3, 2, 4)])

# CAD; GRCh37; use UKB
# foo = read.table("Nelson2017_ST7_304independentvariants.txt", header=TRUE, sep="\t")
# dim(foo)
# colnames(foo)[1:3] = c("rs", "chr", "pos")
# snps = rbind(snps, foo[, c(2, 3, 1, 4)])

# AF not UK Biobank, AllofUs
# https://www.nature.com/articles/ng.3843
foo = read.table("NatGenet_49_946_23novelknown.txt", header=TRUE, sep="\t")
dim(foo); head(foo)
colnames(foo)[1:3] = c("rs", "chr", "EA")
foo$pos = NA
snps = rbind(snps, foo[, c(2, 5, 1, 3)])

snps = snps[order(snps$chr, snps$pos), ]
snps = snps[!duplicated(snps$rs), ]
dim(snps)
# We only care about rs; 610; withHF 633

### Extract from ukb by rs
table(snps$rs %in% ukb$ID) #583; withHF 606
snps[! snps$rs %in% ukb$ID, ]
# if allele frequency <0.01, OK to drop
# rs75305034 was merged into rs55688777
# rs115321690 was merged into rs3132535 (lowfreq)
# rs63418562 was merged into rs9508495
# rs12016871 was merged into rs9581854
# rs62023387 no proxy (not in 1000G). give up.
# rs115231027 was merged into rs4792830
# rs79349575 was merged into rs12941263
# rs141216986 chrX
# rs5945326 chrX
snpsukb = ukb[ukb$ID %in%
                c(snps$rs,
                  "rs55688777", "rs9508495", "rs9581854",
                  "rs4792830", "rs12941263"), ]
dim(snpsukb) #588; withHF 611

### Portable to aou
table(snpsukb$keyrs %in% aou$keyrs) #33 555; withHF 33 578

library(LDlinkR)

toadd = c()
for (snp in snpsukb$ID[! snpsukb$keyrs %in% aou$keyrs]) {
  print(snp)
  proxy = LDproxy(snp, genome_build="grch37", token="769f387ceaa7")
  if(ncol(proxy) == 1) { next }
  proxy = proxy[proxy$R2 >= 0.8, ]
  proxy = proxy[order(proxy$R2, decreasing=TRUE), ]
  proxy$Alleles = sub("\\)", "", sub("\\(", "", proxy$Alleles))
  proxy$A1 = sub("/.*", "", proxy$Alleles)
  proxy$A2 = sub(".*/", "", proxy$Alleles)
  proxy$keyrs = paste(proxy$RS_Number,
                      pmin(proxy$A1, proxy$A2),
                      pmax(proxy$A1, proxy$A2),
                      sep=" ")
  proxy = proxy[proxy$keyrs %in% aou$keyrs, ]
  if (nrow(proxy) > 0) {
    print(paste0("proxy: ", proxy$keyrs[1]))
    toadd = c(toadd, proxy$keyrs[1]) }
}
snpsukbaou = ukb[ukb$keyrs %in%
                   c(intersect(snpsukb$keyrs, aou$keyrs),
                     toadd), ]
dim(snpsukbaou) #570; withHF 593
snpsukbaou$IDaou = aou$ID[match(snpsukbaou$keyrs, aou$keyrs)]

write.table(snpsukb,
            # file="snps_ukb.txt",
            file="snps_ukb.withHF.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(snpsukbaou,
            # file="snps_ukb_aou.txt",
            file="snps_ukb_aou.withHF.txt",
            row.names=FALSE, quote=FALSE, sep="\t")

for (i in 1:22) {
  snpsinchr = snpsukbaou[snpsukbaou$CHROM == i, ]
  write.table(snpsinchr$ID,
              # file=paste0("snps_rs.chr", i, ".txt"),
              file=paste0("snps_rs.chr", i, ".withHF.txt"),
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(snpsinchr[, c("ID", "ALT")],
              # file=paste0("snps_rs_allele.chr", i, ".txt"),
              file=paste0("snps_rs_allele.chr", i, ".withHF.txt"),
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

### https://github.com/GraceSheng/triple-liftOver
### https://doi.org/10.1016/j.xhgg.2022.100159
# tail -n +2 snps_ukb_aou.txt | perl -ne '@a=split; print join(" ", @a[0, 2], 0, @a[1, 3, 4]), "\n"' > snps_ukb_aou.bim
# perl tripleliftover_v133.pl --bim ../snps_ukb_aou.bim --base hg19 --target hg38 --outprefix test
#
# snpid	chr	pos_[Hg38]	pos_[hg19]	category
# rs7948851	11	54600090	51519190	inverted
#
# Invert rs7948851 G/C !!!
