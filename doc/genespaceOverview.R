## ----setup, include=FALSE-----------------------------------------------------
if(dir.exists("~/Desktop/testGenespace"))
  unlink("~/Desktop/testGenespace", recursive = T)
if(!dir.exists("~/Desktop/testGenespace"))
  dir.create("~/Desktop/testGenespace")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "~/Desktop/testGenespace")

## ---- eval = FALSE------------------------------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  if (!requireNamespace("GENESPACE", quietly = TRUE))
#      devtools::install_github("jtlovell/GENESPACE", upgrade = F)

## ---- eval = FALSE------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  if (!requireNamespace("Biostrings", quietly = TRUE))
#    BiocManager::install("Biostrings")
#  if (!requireNamespace("rtracklayer", quietly = TRUE))
#     BiocManager::install("rtracklayer")

## -----------------------------------------------------------------------------
library(GENESPACE)
runwd <- file.path("~/Desktop/testGenespace")

## -----------------------------------------------------------------------------
make_exampleDataDir(writeDir = runwd)

## -----------------------------------------------------------------------------
list.files(runwd, recursive = T, full.names = F)

## -----------------------------------------------------------------------------
gpar <- init_genespace(
  genomeIDs = c("human","chimp","rhesus"),
  speciesIDs = c("human","chimp","rhesus"),
  versionIDs = c("human","chimp","rhesus"),
  ploidy = rep(1,3),
  diamondMode = "fast",
  orthofinderMethod = "fast",
  wd = runwd,
  nCores = 4,
  minPepLen = 50,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2mcscanx = "~/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

## ----parse annotations--------------------------------------------------------
parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=")

## ----orthofinder--------------------------------------------------------------
gpar <- run_orthofinder(
  gsParam = gpar)

## ----synteny------------------------------------------------------------------
gpar <- synteny(gsParam = gpar)

## ----riparian, fig.width=5----------------------------------------------------
ripdat <- plot_riparian(
  gpar,
  colByChrs = c("#BC4F43", "#F67243"))

## ----build pangenome----------------------------------------------------------
pg <- pangenome(gpar)

## -----------------------------------------------------------------------------
query_pangenome(
  pg = pg, 
  refChrom = "3", 
  startOrder = 1, 
  endOrder = 5)

## -----------------------------------------------------------------------------
query_pangenome(
  pg = pg, 
  refChrom = "3", 
  startOrder = 38, 
  endOrder = 42)

