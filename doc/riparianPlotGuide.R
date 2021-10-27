## ----setup, include=FALSE-----------------------------------------------------
if(dir.exists("~/Desktop/testGenespace"))
  unlink("~/Desktop/testGenespace", recursive = T)
if(!dir.exists("~/Desktop/testGenespace"))
  dir.create("~/Desktop/testGenespace")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "~/Desktop/testGenespace")

## ----run, echo = TRUE, error = FALSE, warning = FALSE, message = FALSE, include = FALSE----
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("GENESPACE", quietly = TRUE))
    devtools::install_github("jtlovell/GENESPACE", upgrade = F)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
if (!requireNamespace("rtracklayer", quietly = TRUE))
   BiocManager::install("rtracklayer")


library(GENESPACE)
runwd <- file.path("~/Desktop/testGenespace")

make_exampleDataDir(writeDir = runwd)

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

parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=")

gpar <- run_orthofinder(
  gsParam = gpar)

gpar <- synteny(gsParam = gpar)

## ----riparian, fig.width = 5, fig.height = 3, fig.align = 'center'------------
ripSourceData <- plot_riparian(
  gpar, 
  returnSourceData = T, 
  highlightRef = "cyan")

## ----bp, fig.width = 5, fig.height = 3, fig.align = 'center'------------------
plot_riparian(
  gpar, 
  plotRegions = F, 
  useOrder = F, 
  colByChrs = c("gold", "cyan"),
  braidAlpha = .25,
  verbose = F)

## ----invert, fig.width = 5, fig.height = 3, fig.align = 'center'--------------
invertThisGenomeChr <- data.table(genome = "rhesus", chr = "2")
plot_riparian(
  gpar, 
  invertTheseChrs = invertThisGenomeChr,
  verbose = F)

## ----hchr, fig.width = 5, fig.height = 3, fig.align = 'center'----------------
plot_riparian(
  gpar, 
  onlyTheseChrs = "4")

## ----ichr, fig.width = 5, fig.height = 3, fig.align = 'center'----------------
regs <- data.table(
  genome = c("human", "rhesus"),
  chr = c(3, 5),
  start = c(0, 0),
  end = c(1e7, 1e7))
plot_riparian(
  gpar, 
  useOrder = F, 
  onlyTheseRegions = regs,
  verbose = F)

## ----ic, fig.width = 5, fig.height = 3, fig.align = 'center'------------------
regs <- data.table(
  genome = c("human"),
  chr = c(3),
  start = c(0),
  end = c(5e7))
plot_riparian(
  gpar, 
  useOrder = F, 
  onlyTheseRegions = regs, 
  excludeChrOutOfRegion = F, 
  colByChrs = "white",
  verbose = F)

## ----orange, fig.width = 5, fig.height = 3, fig.align = 'center'--------------
plot_riparian(
  gpar,
  blackBg = FALSE,
  chrFill = "orange",
  chrBorder = "grey",
  braidAlpha = 1,
  verbose = F)

## ----bp1, fig.width = 5, fig.height = 3, fig.align = 'center'-----------------
plot_riparian(
  gpar,
  chrLabCex = 1,
  chrRectBuffer = 1.1,
  verbose = F)

## ----lab, fig.width = 5, fig.height = 3, fig.align = 'center'-----------------
plot_riparian(
  gpar,
  labelTheseGenomes = "human",
  verbose = F)

