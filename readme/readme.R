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
peptide.dir = "/Users/jlovell/Documents/comparative_genomics/peptide"
transcript.dir = "/Users/jlovell/Documents/comparative_genomics/transcript/"
mappings.dir = "/Users/jlovell/Documents/comparative_genomics/mappings"
gff.dir = "/Users/jlovell/Documents/comparative_genomics/gff"
mcscan.dir = "/Users/jlovell/Documents/comparative_genomics/mcs_mapping"
Concorde.path = "/Users/jlovell/Documents/comparative_genomics/programs/concorde/TSP/"
MCScanX.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX"
# Parse the names of the peptide and transcript fastas so that they match geneIDs the gffs
# in.pep = list.files(peptide.dir, pattern = "pep.fa", full.names = T)
# ss = lapply(in.pep, function(i){
#   print(i)
#   x = Biostrings::readAAStringSet(i)
#   names(x)<-sapply(gsub(".*locus=","",names(x)), function(y) strsplit(y," ")[[1]][1])
#   print(head(names(x)))
#   writeXStringSet(x, filepath = i)
# })
#
# in.trs = list.files(transcript.dir, pattern = "transcript.fa", full.names = T)
# ss = lapply(in.trs, function(i){
#   print(i)
#   x = Biostrings::readDNAStringSet(i)
#   names(x)<-sapply(gsub(".*locus=","",names(x)), function(y) strsplit(y," ")[[1]][1])
#   print(head(names(x)))
#   writeXStringSet(x, filepath = i)
# })


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
mappings = pipe_Diamond2MCScanX(inputFileMatrix = id.mat,
                                mcscan.dir = mcscan.dir,
                                nthreads = 6, onlyParseMap = T)

ref.id = "PhalliiHAL"

init.mcscan2 = run_MCScanX(ref.id = "PhalliiHAL",
                           mapping.obj = mappings,
                           MCScanX.path = MCScanX.path,
                           mcs_mapping.dir = mcscan.dir,
                           MCScanX.params = "-a -s 5 -m 10", plotit = T,
                           altGenome2plot = "Sviridis",
                           chr1toplot = c("Chr08"),
                           chr2toplot = c("Chr_08"),
                           buffer = 1)

parsed_mcscan = split_blocksByTSP(map = init.mcscan2$map,
                                  Concorde.path = Concorde.path)

merge_overlappingBlocks = function(blk,map,
                                   verbose = T,
                                   buffer = 1){




init.mcscan2 = data.table(init.mcscan2)

map = data.frame(init.mcscan2$map)
blk = data.frame(init.mcscan2$block)

test = merge_overlappingBlocks(map = map, blk = blk)




    within.it = with(yt, which(rankend1<=x$rankend1 &
                      rankend2<=x$rankend2 &
                      rankstart1>=x$rankstart1 &
                      rankstart2>=x$rankstart2))
    over.right = with(yt, which(rankstart1<=x$rankend1 &
                                rankstart2<=x$rankend2 &
                                rankend1>=x$rankstart1 &
                                rankend2>=x$rankstart2))
    yt[yt$rankend1>=x$rankstart1 &
         yt$rankend2>=x$rankstart2 ,]
  })
}



pep.list = lapply(in.pep, readAAStringSet)
names(pep.list)<-gsub("/","",gsub(".pep.fa","",gsub(".*peptide/","",in.pep)), fixed = T)

trs.list = lapply(in.trs, readDNAStringSet)
names(trs.list)<-gsub("/","",gsub(".transcript.fa","",gsub(".*transcript/","",in.trs)), fixed = T)

block = data.frame(init.mcscan2$block, stringsAsFactors = F)
map = data.frame(init.mcscan2$map, stringsAsFactors = F)
tpb = block[block$genome1 == "PhalliiHAL" & block$genome2 == "Sviridis" &
              grepl("3|7|8",block$chr1) & grepl("3|7|8", block$chr2),]
tp = map[map$genome1 == "PhalliiHAL" & map$genome2 == "Sviridis" &
           grepl("3|7|8",map$chr_genome1) & grepl("3|7|8", map$chr_genome2),]

setnames(tp, gsub("_genome","", colnames(tp)))
setnames(tp, "block", "block.id")
ggplot(tp, aes(x = start1, y = start2, col = as.factor(block.id)))+
  geom_point(pch = ".")+
  geom_segment(data = tpb, aes(x = start1, xend = end1, y = start2, yend = end2))+
  scale_color_manual(values = rep_len(1:8, length.out = length(unique(tp$block.id))), guide = F)+
  facet_grid(chr1 ~ chr2, as.table = F)

lapply(1:nrow(block), function(i){
  x = block[i,]
  gff1 = gff.list[[x$genome1]]
  gff1 = gff1[with(gff1, chr == x$chr1 & start <= x$end1 & end >= x$start1),]

  gff2 = gff.list[[x$genome2]]
  gff2 = gff2[with(gff2, chr == x$chr2 & start <= x$end2 & end >= x$start2),]

  pep1 = pep.list[[x$genome1]]
  pep1 = pep1[gff1$id]

  pep2 = pep.list[[x$genome2]]
  pep2 = pep2[gff2$id]

  trs1 = trs.list[[x$genome1]]
  trs1 = trs1[gff1$id]

  trs2 = trs.list[[x$genome2]]
  trs2 = trs2[gff2$id]



})


make_blockTrs = function()

  allmap = rbindlist(lapply(mappings, function(x){
    out1 = x[,1:5]
    out2 = x[,6:10]
    setnames(out1,c("geneID","chr","start","end","strand"))
    setnames(out2,c("geneID","chr","start","end","strand"))
    return(rbind(out1,out2))
  }))
allmap = allmap[!duplicated(allmap),]
allmap1 = allmap
allmap2 = allmap

