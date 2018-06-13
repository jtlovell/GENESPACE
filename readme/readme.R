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

## get into py 2.7
# conda create -n py27 python=2.7 anaconda
# conda config --add channels defaults
# conda config --add channels conda-forge
# conda config --add channels bioconda
# conda install orthofinder
# conda install diamond

## make sure to add mcscanx to your path
# export PATH=$PATH:/Users/jlovell/Documents/comparative_genomics/programs/MCScanX/


# library(compareGeneSpace)

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

MCScanX.params = "-a -s 5 -m 25"
buffer = 1

tsp.method = "Concorde"
max.jump = 5
ref.id = "PhalliiHAL"
min.block.size = 5

plotit = T
chr1toplot = c("Chr03","Chr07","Chr08")
chr2toplot = c("Chr_03","Chr_07","Chr_08")
chr1toplot = c("Chr07","Chr08")
chr2toplot = c("Chr_07","Chr_08")
altGenome2plot = "Sviridis"

verbose = T

mappings = pipe_Diamond2MCScanX(inputFileMatrix = id.mat,
                                mcscan.dir = mcscan.dir,
                                nthreads = 6, onlyParseMap = T,
                                topPerc = 50,
                                dbs_radii = c(100, 40, 15),
                                dbs_mappingsInRadius = c(10, 5, 5))

mcs.out = run_MCScanX(ref.id = ref.id,
                      MCScanX.path = MCScanX.path,
                      mcs_mapping.dir = mcscan.dir,
                      MCScanX.params = MCScanX.params,
                      buffer = buffer)

plot_blocksAndMapping(map = mcs.out$map,
                      blk= mcs.out$block,
                      ref.id = ref.id,
                      altGenome2plot = altGenome2plot,
                      chr1toplot = paste0("Chr0",9),
                      chr2toplot =  paste0("Chr_0",9),
                      main = "Raw MCScanX blocks")

merged = merge_overlappingBlocks(map = mcs.out$map,
                                 blk = mcs.out$block,buffer = 0,
                                 verbose = verbose)
plot_blocksAndMapping(map = merged$map,
                      blk= merged$block,
                      ref.id = ref.id,
                      altGenome2plot = altGenome2plot,
                      chr1toplot = paste0("Chr0",1:9),
                      chr2toplot =  paste0("Chr_0",1:9),
                      main = "Merged MCScanX blocks")



spled <-split_byDensity(map = merged$map,quantile = .999)

spled2 <-split_byDensity(map = spled$map,quantile = .999,rerank = F)
if(plotit){
  plot_blocksAndMapping(map = spled$map,
                        blk= spled$block,
                        ref.id = ref.id,
                        altGenome2plot = altGenome2plot,
                        chr1toplot = paste0("Chr0",1:9),
                        chr2toplot =  paste0("Chr_0",1:9),
                        main = "Split MCScanX blocks Final")
}

final_blocks = chop_genome(gff, fasta, block, map, buffer, ref.id){
  block = spled$b
}
