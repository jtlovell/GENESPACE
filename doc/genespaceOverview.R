## ----require------------------------------------------------------------------
library(GENESPACE)
runwd <- file.path("~/Desktop/testGenespace")

## ----rmExistDir, echo = F-----------------------------------------------------
if(dir.exists(runwd))
  unlink(runwd, recursive = T)

## ----mkdat--------------------------------------------------------------------
make_exampleDataDir(writeDir = runwd)

## ----showFiles----------------------------------------------------------------
list.files(runwd, recursive = T, full.names = F)

## ----init---------------------------------------------------------------------
gpar <- init_genespace(
  genomeIDs = c("human","chimp","rhesus"),
  speciesIDs = c("human","chimp","rhesus"),
  versionIDs = c("human","chimp","rhesus"),
  outgroup = NULL,
  ploidy = rep(1,3),
  diamondMode = "fast",
  orthofinderMethod = "fast",
  wd = runwd,
  orthofinderInBlk = FALSE, 
  overwrite = F, 
  verbose = T,
  nCores = 4,
  minPepLen = 50,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = "~/Documents/comparative_genomics/programs/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

## ----qInit, eval = FALSE------------------------------------------------------
#  spec <- c("human","chimp","rhesus")
#  gpar <- init_genespace(
#    genomeIDs = spec,  speciesIDs = spec,  versionIDs = spec, ploidy = rep(1,3),
#    diamondMode = "fast", orthofinderMethod = "fast", wd = runwd,
#    path2mcscanx = "~/Documents/comparative_genomics/programs/MCScanX",
#    rawGenomeDir = file.path(runwd, "rawGenomes"))

## ----parseAnnotations---------------------------------------------------------
parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=")

## ----parseNCBI, eval = FALSE--------------------------------------------------
#  parse_ncbi(gsParam = gpar, overwrite = T)

## ----orthofinder--------------------------------------------------------------
gpar <- run_orthofinder(
  gsParam = gpar)

## ----setParam-----------------------------------------------------------------
gpar <- set_syntenyParams(gsParam = gpar)

## ----prntParams, echo = FALSE-------------------------------------------------
knitr::kable(gpar$params$synteny)

## ----synteny------------------------------------------------------------------
gpar <- synteny(gsParam = gpar)

## ----riparian, fig.width = 10-------------------------------------------------
plot_riparianHits(gpar)

## ----riparian2, fig.width = 10------------------------------------------------
ripSouceDat <- plot_riparianHits(
  gpar, 
  refGenome = "chimp",
  invertTheseChrs = data.frame(genome = "rhesus", chr = 2),
  genomeIDs = c("chimp", "human", "rhesus"),
  labelTheseGenomes = c("chimp", "rhesus"),
  gapProp = .001,
  refChrCols = c("#BC4F43", "#F67243"),
  blackBg = FALSE, 
  returnSourceData = T, 
  verbose = F)

## ----ripReg, fig.width = 10---------------------------------------------------
regs <- data.frame(
  genome = c("human", "human", "chimp", "rhesus"),
  chr = c(3, 3, 4, 5),
  start = c(0, 50e6, 0, 60e6),
  end = c(10e6, 70e6, 50e6, 90e6),
  cols = c("pink", "gold", "cyan", "dodgerblue"))

plot_riparianHits(
  gpar, 
  onlyTheseRegions = regs)

## ----ripReg2, fig.width = 10--------------------------------------------------
regs2 <- data.frame(
  genome = c("human", "human"),
  chr = c(3, 3),
  start = c(0, 50e6),
  end = c(10e6, 70e6),
  cols = c("pink", "gold"))

plot_riparianHits(
  gpar, 
  onlyTheseRegions = regs2, 
  excludeNoRegChr = T)

## ----ripOvly, fig.width = 10--------------------------------------------------
plot_riparianHits(
  gpar, 
  refChrCols = "grey80", 
  annotatePlot = F)

plot_riparianHits(
  gpar, 
  onlyTheseRegions = regs,
  add2plot = T)

## ----buildPangenome-----------------------------------------------------------
pg <- pangenome(gpar)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

