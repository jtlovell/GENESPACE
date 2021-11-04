# 1. map the positions of var
nodels <- proc_garbageStatsFormat(
  "proc_new_full_set.condensedGenes.GT2_reads.LE90_inds.noINDELS.stats.dat",
  splitInd2list = T)
d <- nodels[,list(lib = unlist(localIndividuals)),
            by = c("chr", "start", "end", "nInd","localVariants","uniqID","type","geneName")]
maxAFthresh <- with(d, (.9*uniqueN(lib)))
ulibs <- unique(d$lib)
# -- Pull all sites at < 95% of individuals
tp <- with(subset(nodels, nInd <= maxAFthresh & nInd > 0 & (end-start) > 1), data.frame(
  Chr = chr, Start = start, End = end,
  Value = ifelse(nInd <= 2, 1, ifelse(nInd <= 6, 3, ifelse(nInd < 20, 4, 5)))))

# -- plot full karyotype
ideogram(
  karyotype = human_karyotype,
  overlaid = tp, colorset1 = c("red","gold","lightblue","dodgerblue"),
  output = "final_noINDEL_fullAF.svg")

# -- sites at 10% or fewer
ideogram(
  karyotype = human_karyotype,
  overlaid = subset(tp, Value < 4), colorset1 = c("red"),
  output = "final_noINDEL_0.1AF.svg")

# -- plot karyotype with just private sites
ideogram(
  karyotype = human_karyotype,
  overlaid = subset(tp, Value == 1), colorset1 = c("red"),
  output = "final_noINDEL_privAF.svg")

# -- for each library pull the number of private sites with n total libraries
# -- get mean for sample of 100 combinations for each
i <- "NA19030"

k <- 10
nsamp = 10
d <- nodels[,list(lib = unlist(localIndividuals)),
            by = c("chr", "start", "end", "nInd","localVariants","uniqID","type","geneName")]
ind <- d[,c("lib", "uniqID")]
spl <- split(ind, by = "lib")
testat <- c(
  1:20,
  unique(GENESPACE::round_toInteger(22:38,2)),
  42, 46, 50, 55, 61, 68, 76, 85, 95, 106)
mx <- max(mx) + 1
nsamp <- ceiling(ifelse(testat < 20, (mx - testat),
                        ifelse(testat < 40, (mx - testat)/2, (mx - testat)/8))/2)
setDTthreads(1)
sens <- rbindlist(mclapply(ulibs, mc.cores = 6, mc.preschedule = F, function(i){
  js <- ulibs[ulibs != i]
  kd <- rbindlist(lapply(1:length(testat), function(k){
    nsampk <- nsamp[k]
    nk = testat[k]
    sampd <- rbindlist(lapply(1:nsampk, function(samp){
      x <- spl[[i]]
      y <- rbindlist(spl[sample(js, k, replace = F)])
      out <- data.table(
        lib = i,
        k = nk,
        nsamp = samp,
        npriv = sum(!x$uniqID %in% y$uniqID),
        nshared = sum(x$uniqID %in% y$uniqID))
      return(out)
    }))
    return(sampd)
  }))
  ko <- kd[,list(prop = sum(npriv)/sum(nshared + npriv),
                 meanPriv = mean(npriv),
                 meanShared = mean(nshared)),
           by = c("lib","k")]
  return(ko)
}))
fwrite(sens, file = "privatePAV_sens.csv.gz")

libs <- c(
  "CSER-011-C", "CSER-105-C", "CSER-110-C", "CSER-221-C", "CSER-225-C", "CSER-238-C",
  "CSER-242-C", "CSER-A-00002-C", "CSER-A-00004-C", "CSER-A-00005-C")
ggplot(subset(sens, lib %in% libs), aes(x = k, y = meanPriv, col = lib))+
  geom_line()+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "grey5",
                                        color = "white"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 2, size = .15),
        strip.background = element_rect(fill="white",
                                        colour = "white"))+
  labs(x = "sequenced population size (n genotypes)",
       y = "Number of 'private' SVs",
       title = "Sensitivity of private SV to sequenced pop. size")
library(ggplot2)
pdf("privateSV_sens3.pdf", height = 3, width = 5)
# ggplot(sens, aes(x = k, y = meanPriv, group = lib, col = lib == "CSER-242-C"))+
#   geom_line(alpha = .2)+
#   scale_color_manual(values = c("white","green"), guide = F)+
#   theme(panel.background = element_rect(fill = "grey5",
#                                         color = "white"),
#         panel.border = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 2, size = .15),
#         strip.background = element_rect(fill="white",
#                                         colour = "white"))+
#   labs(x = "sequenced population size (n genotypes)",
#        y = "Number of 'private' SVs",
#        title = "Sensitivity of private SV to sequenced pop. size")
ggplot(sens, aes(x = k, y = meanPriv, group = lib, col = lib == "CSER-242-C"))+
  geom_line(alpha = .2)+
  scale_color_manual(values = c("white","green"), guide = F)+
  theme(panel.background = element_rect(fill = "grey5",
                                        color = "white"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 2, size = .15),
        strip.background = element_rect(fill="white",
                                        colour = "white"))+
  scale_y_log10()+
  labs(x = "sequenced population size (n genotypes)",
       y = "Number of 'private' SVs",
       title = "Sensitivity of private SV to sequenced pop. size")
dev.off()

maxAFthresh <- 1e3
ds <- subset(d, (end-start) > 1 & lib == "CSER-242-C" & nInd > 0)
tp <- with(subset(ds, nInd <= maxAFthresh), data.frame(
  Chr = chr, Start = start, End = end, Value = 1))

ideogram(
  karyotype = human_karyotype,
  overlaid = tp, colorset1 = rgb(1,0,0,.2),
  output = "final_noINDEL_fullAF_CSER-242-C.svg")

# -- sites at 10% or fewer
maxAFthresh <- 35
tp <- with(subset(ds, nInd <= maxAFthresh), data.frame(
  Chr = chr, Start = start, End = end, Value = 1))

ideogram(
  karyotype = human_karyotype,
  overlaid = tp, colorset1 = rgb(1,0,0,.2),
  output = "final_noINDEL_35lib_CSER-242-C.svg")


# -- plot karyotype with just private sites
maxAFthresh <- 1
tp <- with(subset(ds, nInd <= maxAFthresh), data.frame(
  Chr = chr, Start = start, End = end, Value = 1))

ideogram(
  karyotype = human_karyotype,
  overlaid = tp, colorset1 = "red",
  output = "final_noINDEL_1lib_CSER-242-C_test.svg")



ds <-
tp <- with(subset(d, (end-start) > 1 & lib %in% libs & nInd == 1), data.frame(
  Chr = chr, Start = start, End = end, Value = as.numeric(as.factor(lib))))
cols <- c(
  "#8B1D00", "#D05100", "#ED9004", "#F9C70E", "#EAE075", "#BAE0DB",
  "#8BEDF9", "#74B8FC", "#4871F9", "#040DC9","#0E004C","#5E09A3",
  "#C054F9","#E6BDFC")[c(1,2,3,4,7,8,10,12,13,14)]
ideogram(
  karyotype = human_karyotype,
  overlaid = tp, colorset1 = cols,
  output = "final_noINDEL_allPriv.svg", width = 400)


# -- plot curves

# -- pull karyotype of chr19
maxAFthresh <- 1
ds <- subset(d, (end-start) > 1 & lib %in% libs & nInd == 1 & chr == 19)
c19 <- human_karyotype[19,]
for(i in libs){
  print(i)
  mi <- with(subset(ds, lib == i), data.frame(
    Chr = chr, Start = start, End = end, Value = 1))
  col <- cols[i]
  if(nrow(mi) > 1){
    ideogram(
      karyotype = c19,
      overlaid = mi,
      output = sprintf("%s_nodels19.svg", i),
      width = 10)
  }
}
tp <- with(subset(ds, nInd <= maxAFthresh), data.frame(
  Chr = chr, Start = start, End = end, Value = 1))

# -- plot positions of unique variants for 10 libraries.

