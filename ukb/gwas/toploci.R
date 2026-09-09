library(dplyr)

# r = "regenie_lipidaemiadrug_onlyaffectedexclsupplement/ukb_lipidaemiadrug_onlyaffectedexclsupplement_step2_BT_chr"
# t = "C10AA"
# t = "C10AB"
# t = "C10AX09"
r = "regenie_hypertensiondrug_onlyaffectedDBPge100SBPge160/ukb_hypertensiondrug_onlyaffectedDBPge100SBPge160_step2_BT_chr"
# t = "C03"
# t = "C07"
# t = "C08"
t = "C09"

data = data.table::fread(
  cmd = paste0("xzcat ", r,
               1, "_",
               t, ".regenie.xz"),
  sep = " ") %>%
  filter(LOG10P > 2)
for (i in 2:22) {
  data = rbind(
    data,
    data.table::fread(
      cmd = paste0("xzcat ", r,
                   i, "_",
                   t, ".regenie.xz"),
      sep = " ") %>%
      filter(LOG10P > 2)
  )
}

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

data = data %>% filter(LOG10P > 3)
x = snpstowindow(data$CHROM, data$GENPOS, data$LOG10P)
x = x %>%
  dplyr::arrange(- pow) %>%
  dplyr::top_n(20, pow)

output = data[match(paste0(x$chr, ":", x$peak),
                    paste0(data$CHROM, ":", data$GENPOS)), ]
output = output %>%
  mutate(A1isminor = A1FREQ < 0.5,
         EA    = ifelse(A1isminor, ALLELE1, ALLELE0),
         OA    = ifelse(A1isminor, ALLELE0, ALLELE1),
         EAF   = ifelse(A1isminor, A1FREQ, 1 - A1FREQ),
         logOR = ifelse(A1isminor, BETA, - BETA),
         P = 10^(- LOG10P),
         drug = t) %>%
  dplyr::select(drug, CHROM, GENPOS, ID, EA, OA, EAF, logOR, P)
write.csv(output, file="foo.csv")