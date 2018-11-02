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




# Pipeline
1. get representative peptide for each gene model, rename to gene ID, place in single directory
2. run orthofinder w/ diamond
3. read in orthofinder blast
4. Parse orthofinder orthoblock file, calculate n / genome for each block
5. drop orthofinder "single-copy" genes from blast
6. drop orthofinder "multi-copy" genes above expected threshold (2x ploidy)
7. Parse orthofinder blast files, concatenate, reformat for MCScanX
8. Run MCScanX, merging and splitting blocks to ensure no overlapping.
9. Parse genes in orthogroups by blocks, dropping genes that are not in the same blocks
10. Subset to orthogroups that
10. Parse orthogroups, extracting























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

remap.dir = "/Users/jlovell/Documents/comparative_genomics/re_mapping/"
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
  remap.dir = "/Users/jlovell/Documents/comparative_genomics/re_mapping",
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
                      MCScanX.params = "-a -s 5 -m 25",
                      buffer = 1)

plot_blocksAndMapping(map = mcs.out$map,
                      blk= mcs.out$block,
                      ref.id = "PhalliiHAL",
                      altGenome2plot = "Sviridis",
                      chr1toplot = paste0("Chr0",1:9),
                      chr2toplot =  paste0("Chr_0",1:9),
                      main = "Raw MCScanX blocks")

merged = merge_overlappingBlocks(map = mcs.out$map,
                                 blk = mcs.out$block,buffer = 1,
                                 verbose = TRUE)
plot_blocksAndMapping(map = merged$map,
                      blk= merged$block,
                      ref.id = "PhalliiHAL",
                      altGenome2plot = "Sviridis",
                      chr1toplot = paste0("Chr0",7),
                      chr2toplot =  paste0("Chr_0",7),
                      main = "Raw MCScanX blocks")



spled <-split_byDensity(map = merged$map,
                        quantile = .999)

if(plotit){
  plot_blocksAndMapping(map = spled$map,
                        blk= spled$block,
                        ref.id = "PhalliiHAL",
                        altGenome2plot = "Sviridis",
                        chr1toplot = paste0("Chr0",7),
                        chr2toplot =  paste0("Chr_0",7),
                        main = "Raw MCScanX blocks")

}


tmp = make_blocks(spled$map)
map = tmp$map
blk = tmp$block


gr.ref = with(blk[blk$genome1 == ref.id,],
              GRanges(seqnames = chr1,
                      ranges = IRanges(start = start1,
                                       end = end1),
                      names = block.id))
gr.q = gr.ref[1,]
test = findOverlaps(gr.q, gr.ref)

gr.alt = with(blk[blk$genome1 == ref.id,],
              GRanges(seqnames = chr2,
                      ranges = IRanges(start = start2,
                                       end = end2,
                                       names = block.id)))

cat_blocks = function(blk, map){

}
final_blocks = chop_genome(id.mat,
                           block,
                           map,
                           buffer=2,
                           ref.id){

  if(verbose)
    cat("Reading in annotations\n")

  gff.u = c(id.mat$gff1, id.mat$gff2)[!duplicated(c(id.mat$gff1, id.mat$gff2))]
  names(gff.u)<-c(id.mat$id1, id.mat$id2)[!duplicated(c(id.mat$gff1, id.mat$gff2))]
  gff.list = sapply(gff.u, simplify = F, USE.NAMES = T, parse_quickGff)

  all.genes = unique(unlist(lapply(gff.list, function(x) as.character(data.frame(x)[,1]))))
  orphan = all.genes[!all.genes %in% unique(c(map$geneID1,map$geneID2))]

  gff.orphan = rbindlist(lapply(names(gff.list), function(i){
    x = gff.list[[i]]
    x$genome1 = i
    x = x[x$id %in% orphan, ]
    return(x)
  }))

  bp.buffer = 10000

  blk$block.id = as.character(blk$block.id)
  mapdata = rbindlist(lapply(1:nrow(blk), function(i){
    if(verbose & i %% 100 == 0)
      cat("completed", i, "/", nrow(blk),"\n")

    block.id = blk$block.id[i]
    genome1 = blk$genome1[i]
    genome2 = blk$genome2[i]
    start1 = blk$start1[i]-bp.buffer
    start2 = blk$start2[i]-bp.buffer
    end1 = blk$end1[i]+bp.buffer
    end2 = blk$end2[i]+bp.buffer
    chr1 = blk$chr1[i]
    chr2 = blk$chr2[i]

    wh1 = which(gff.orphan$genome1 == genome1 &
      gff.orphan$chr == chr1 &
      gff.orphan$end>= start1 &
      gff.orphan$start<=end1)

    if(length(wh1)>0){
      coords1 = gff.orphan$id[wh1]
    }else{
      coords1 = NA
    }

    wh2 = which(gff.orphan$genome1 == genome2 &
      gff.orphan$chr == chr2 &
      gff.orphan$end>= start2 &
      gff.orphan$start<=end2)

    if(length(wh2)>0){
      coords2 = gff.orphan$id[wh2]
    }else{
      coords2 = NA
    }

    o1 = data.frame(block.id = block.id,
                    genome.id = genome1,
                    gene.id = coords1,
                    is.subject = F,
                    class = "orphan.inblock",
                    stringsAsFactors = F)
    o2 = data.frame(block.id = block.id,
                    genome.id = genome2,
                    gene.id = coords2,
                    class = "orphan.inblock",
                    is.subject = T,
                    stringsAsFactors = F)
    i1 = data.frame(block.id = block.id,
                    genome.id = genome1,
                    gene.id = map$geneID1[map$block.id == block.id],
                    class = "mapped",
                    is.subject = T,
                    stringsAsFactors = F)
    i2 = data.frame(block.id = block.id,
                    genome.id = genome2,
                    gene.id = map$geneID2[map$block.id == block.id],
                    class = "mapped",
                    is.subject = T,
                    stringsAsFactors = F)


    out = data.frame(rbind(i1,o1,i2,o2))
    out = out[complete.cases(out),]
    return(out)
  }))

  upep = unique(c(id.mat$remap.pep1, id.mat$remap.pep2))
  library(Rsamtools)
  indexed = lapply(upep, indexFa)
  fal = lapply(upep, FaFile)
  seq = lapply(fal, scanFa, as = "AAStringSet")
  seq = do.call(c, seq)

  spl = split(mapdata, mapdata$block.id)

  odir = "/Users/jlovell/Documents/comparative_genomics/re_mapping"
  spl = spl[sapply(spl, function(y) length(unique(y$genome.id))>1)]
  n = names(spl)[order(-sapply(spl, nrow))]

  com = paste("orthofinder -f",peptide.dir,"-S diamond -t 32")

  orthos = sapply(n, simplify = F, USE.NAMES = T, function(i){


    x = spl[[i]]
    if(verbose)
      cat(which(names(spl)==i),"/", length(spl),"... n.genes = ",nrow(x),"\n")
    rid = ref.id
    if(!rid %in% x$gene.id){
      rid = x$genome.id[1]
    }
    alt.id = unique(x$genome.id[x$genome.id!=rid])
    x.ref = x[x$genome.id == rid,]
    x.alt = x[x$genome.id == alt.id,]
    s1= seq[as.character(x.ref$gene.id)]
    s2= seq[as.character(x.alt$gene.id)]
    wdir = file.path(odir,i)
    system(paste("mkdir",wdir))
    writeXStringSet(s1, file.path(wdir,paste0(rid,".fa")))
    writeXStringSet(s2, file.path(wdir,paste0(alt.id,".fa")))

    com = paste("orthofinder -f",wdir,"-og -S diamond -t 1 > /dev/null")
    system(com)
    of = list.files(wdir,
                    pattern = "^Orthogroups.txt$",
                    recursive = T,
                    full.names = T)
    of = readLines(of)
    of2 = lapply(of, function(y) strsplit(y,": ")[[1]])
    ortholist = lapply(of2, function(y) strsplit(y[[2]], " ")[[1]])
    names(ortholist)<-sapply(of2, function(y) y[1])
    system(paste("rm -r",wdir))
    return(ortholist)
  })
  system()

  extend_blocks = function(blk, map, gff.list){


  }
  find_blocks4orphans = function(genome, gff.list, start, end, rank.buffer = )

  block = spled$block
  map = spled$map
}
