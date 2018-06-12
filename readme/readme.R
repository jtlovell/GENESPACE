# Pipeline:

#############################################
#############################################
# 1. Set up all variables you will need:
rm(list = ls())

library(data.table)
library(dbscan)
library(Biostrings)
library(wordspace)
library(TSP)

library(compareGeneSpace)

peptide.dir = "/Users/jlovell/Documents/comparative_genomics/peptide"
transcript.dir = "/Users/jlovell/Documents/comparative_genomics/transcript/"
mappings.dir = "/Users/jlovell/Documents/comparative_genomics/mappings"
gff.dir = "/Users/jlovell/Documents/comparative_genomics/gff"
mcscan.dir = "/Users/jlovell/Documents/comparative_genomics/mcs_mapping"
Concorde.path = "/Users/jlovell/Documents/comparative_genomics/programs/concorde/TSP/"
MCScanX.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX"
# Parse the names of the peptide and transcript fastas so that they match geneIDs the gffs

# parse_fastaHeader(fasta.dir = peptide.dir, is.peptide = T,
#                   pattern = "pep.fa", verbose = T)


abbrev.dict = list(Phallii = "Pa",
                   PhalliiHAL = "Ph",
                   Sviridis = "Sv",
                   Sbicolor = "Sb",
                   Pvirgatum = "Pv")
ploidy.dict = list(Phallii = 2,
                   PhalliiHAL = 2,
                   Sviridis = 2,
                   Sbicolor = 2,
                   Pvirgatum = 4)
id.mat = make_inputFileMatrix(
  peptide.dir = "/Users/jlovell/Documents/comparative_genomics/peptide",
  mappings.dir = "/Users/jlovell/Documents/comparative_genomics/mappings",
  gff.dir = "/Users/jlovell/Documents/comparative_genomics/gff",
  mcscan.dir = "/Users/jlovell/Documents/comparative_genomics/mcs_mapping",
  ploidy.dict = ploidy.dict,
  abbrev.dict = abbrev.dict,
  ref.id = "PhalliiHAL"
)

id.mat = id.mat[id.mat$id1 != "Pvirgatum" & id.mat$id2 != "Pvirgatum",]

MCScanX.params = "-a -s 5 -m 5"
buffer = 1

tsp.method = "Concorde"
Concorde.path = NULL
max.jump = 5
ref.id = "PhalliiHAL"
min.block.size = 5

plotit = T
chr1toplot = "Chr03"
chr2toplot = "Chr_03"
altGenome2plot = "Sviridis"

verbose = T

mappings = pipe_Diamond2MCScanX(inputFileMatrix = id.mat,
                                mcscan.dir = mcscan.dir,
                                nthreads = 6, onlyParseMap = T,
                                topPerc = 75,
                                dbs_radii = c(100, 40, 20),
                                dbs_mappingsInRadius = c(10, 10, 10))

mcs.out = run_MCScanX(ref.id = ref.id,
                      MCScanX.path = MCScanX.path,
                      mcs_mapping.dir = mcscan.dir,
                      MCScanX.params = MCScanX.params,
                      buffer = buffer)

plot_blocksAndMapping(map = mcs.out$map,
                      blk= mcs.out$block,
                      ref.id = ref.id,
                      altGenome2plot = altGenome2plot,
                      chr1toplot = chr1toplot,
                      chr2toplot = chr2toplot,
                      main = "Raw MCScanX blocks")

merged = merge_overlappingBlocks(map = mcs.out$map,
                                 blk = mcs.out$block,
                                 verbose = verbose)
plot_blocksAndMapping(map = merged$map,
                      blk= merged$blk,
                      ref.id = ref.id,
                      altGenome2plot = altGenome2plot,
                      chr1toplot = chr1toplot,
                      chr2toplot = chr2toplot,
                      main = "Merged MCScanX blocks")


tsped <- split_blocksByTSP(map = merged$map,
                           Concorde.path = Concorde.path)

if(plotit){
  plot_blocksAndMapping(map = tsped$map,
                        blk= tsped$block,
                        ref.id = ref.id,
                        altGenome2plot = altGenome2plot,
                        chr1toplot = chr1toplot,
                        chr2toplot = chr2toplot,
                        main = "Merged MCScanX blocks")
}
