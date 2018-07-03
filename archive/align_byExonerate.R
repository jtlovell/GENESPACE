####
# - Part 0: Format peptide fastas, place into `peptide.dir`
# - Part 1: Run BLASTs with ortho.finder


peptide.dir = "/Users/jlovell/Documents/comparative_genomics/peptide"
tmp.dir = "/Users/jlovell/Documents/comparative_genomics/tmp"
blast.dir = "/Users/jlovell/Documents/comparative_genomics/all_blast"
gff.dir = "/Users/jlovell/Documents/comparative_genomics/gff"

results.dir = "/Users/jlovell/Documents/comparative_genomics/results"


ploidy = c(2,2,2,2,4)
names(ploidy) = c("PhalliiHAL","Phallii","Sviridis","Sbicolor","Pvirgatum")

abbrevs = c("Ph","Pa","Sv","Sb","Pv")
names(abbrevs) = c("PhalliiHAL","Phallii","Sviridis","Sbicolor","Pvirgatum")

mcscanx.input.dir = "/Users/jlovell/Documents/comparative_genomics/mcscanx_input/"



library(data.table)
blast.results = runParse_orthofinder(runOF = FALSE,
                            peptide.dir= peptide.dir,
                            tmp.dir = tmp.dir,
                            blast.dir = blast.dir,
                            gff.dir = gff.dir,
                            ploidy = ploidy)
mcscan.raw = run_MCScanX(blast.results = blast.results,
                         abbrevs = abbrevs,
                         mcscanx.input.dir = mcscanx.input.dir,
                         MCScanX.path = MCScanX.path,
                         MCScanX.params = MCScanX.params)
mcscan.blocks = make_blocks(mcscan.raw)
merge1 = merge_overlappingBlocks(map = mcscan.blocks$map,
                                 blk = mcscan.blocks$block,
                                 buffer = 1.5,
                                 verbose = F)
merged = merge_overlappingBlocks(map = merge1$map,
                                 blk = merge1$block,
                                 buffer = 1.5,
                                 verbose = T)


tmp.dir = "/Users/jlovell/Documents/comparative_genomics/within_blocks"
output.dir = file.path(results.dir,"orthogroups_byblock")
blk = merged$block
buffer = 50000
ogs.bygroup = run_orthoFinderInBlock(blk = merged$block,
                                     tmp.dir = "/Users/jlovell/Documents/comparative_genomics/within_blocks",
                                     output.dir = file.path(results.dir,"orthogroups_byblock"),
                                     buffer = 50000,
                                     peptide.dir = peptide.dir,
                                     gff.dir = gff.dir,
                                     ncores = 8, onlyParseOF = T)

gff.files = list.files(gff.dir, full.names = T)
names(gff.files) = gsub(".gff3$","",basename(gff.files))

parse_gff = function(gff){
  g = suppressWarnings(
    data.table::fread(gff,showProgress = F, verbose = F))
  g = g[g$V3 == "gene",c(9,1,4,5,7)]
  g$V9 = sapply(g$V9, function(x) gsub("Name=","",strsplit(x,";")[[1]][2]))
  data.table::setnames(g, c("id","chr","start","end","strand"))
  return(g)
}

gff = sapply(names(gff.files), USE.NAMES = T, simplify = F, function(i){
  tmp = parse_gff(gff.files[[i]])
  tmp$genome = i
  tmp$order = frank(tmp[,c("chr","start")], ties.method = "random")
  return(tmp)
})

single.genes = lapply(ogs.bygroup$orthogroups, function(x) x[["1"]])


b1 = blk
b1 = b1[,c("block.id","genome1","genome2","chr1","chr2","start1","start2","end1","end2")]
b2 = blk[,c("block.id","genome2","genome1","chr2","chr1","start2","start1","end2","end1")]
colnames(b2)<-c("block.id","genome1","genome2","chr1","chr2","start1","start2","end1","end2")
blko = rbind(b1,b2)

bed.dir = file.path(results.dir,"block.beds")
block.fasta.dir = file.path(results.dir,"block.fasta")
assembly.dir = "/Users/jlovell/Documents/comparative_genomics/genomes"
system(paste("mkdir", bed.dir))
system(paste("mkdir", block.fasta.dir))

spl = split(blko, blko$genome1)
trsh = lapply(names(spl),function(i){
  x = spl[[i]][,c("chr1","start1","end1","block.id")]
  bedf = file.path(bed.dir,paste0(i,".bed"))
  faf = file.path(assembly.dir,paste0(i,".fa"))
  fafo = file.path(block.fasta.dir,paste0(i,".fa"))
  write.table(x, file=bedf,
              quote=F, sep="\t", row.names=F, col.names=F)
  system(paste("bedtools getfasta -fi",
               faf, "-bed", bedf,"-name -fo",fafo))
})

system(paste("mkdir",file.path(tmp.dir,"single.gene.fasta")))

peptide.fastas = lapply(list.files(peptide.dir, full.names = T), function(i){
  Biostrings::readAAStringSet(i)
})
names(peptide.fastas)<-gsub(".fa", "", list.files(peptide.dir))

cds.fastas = lapply(list.files(cds.dir, full.names = T), function(i){
  Biostrings::readAAStringSet(i)
})
names(cds.fastas)<-gsub(".fa", "", list.files(cds.dir))

blk.fastas = lapply(list.files(block.fasta.dir, full.names = T), function(i){
  Biostrings::readDNAStringSet(i)
})
names(blk.fastas)<-gsub(".fa", "", list.files(block.fasta.dir))

trsh = lapply(names(single.genes), function(i){
  x = as.character(single.genes[[i]])
  tmp.fa = file.path(tmp.dir,"single.gene.fasta",i)
  system(paste("mkdir", tmp.fa))

  blk.info = blk[blk$block.id == i,]

  g1 = blk.info$genome1[1]
  g2 = blk.info$genome2[1]

  gblk.info = blk[blk$block.id == i,]
  map.info = map[map$block.id == i,]

  s1 = x[x %in% gff[[g1]]$id]
  s2 = x[x %in% gff[[g2]]$id]

  tc1 = cds.fastas[[g1]][s1]
  tc2 = cds.fastas[[g2]][s2]
  fc1 = file.path(tmp.fa,paste0(g1,".cds.fa"))
  fc2 = file.path(tmp.fa,paste0(g2,".cds.fa"))

  Biostrings::writeXStringSet(tc1, filepath = fc1)
  Biostrings::writeXStringSet(tc2, filepath = fc2)

  t1 = peptide.fastas[[g1]][s1]
  t2 = peptide.fastas[[g2]][s2]
  f1 = file.path(tmp.fa,paste0(g1,".pep.fa"))
  f2 = file.path(tmp.fa,paste0(g2,".pep.fa"))

  Biostrings::writeXStringSet(t1, filepath = f1)
  Biostrings::writeXStringSet(t2, filepath = f2)

  wg1 = blk.fastas[[g1]][i]
  wg2 = blk.fastas[[g2]][i]
  fg1 = file.path(tmp.fa,paste0(g1,".blk.fa"))
  fg2 = file.path(tmp.fa,paste0(g2,".blk.fa"))

  Biostrings::writeXStringSet(wg1, filepath = fg1)
  Biostrings::writeXStringSet(wg2, filepath = fg2)

  blat1 = file.path(tmp.fa,paste0(g1,"_",g2,".blast"))
  bed1 = file.path(tmp.fa,paste0(g1,"_",g2,".bed"))
  blfa1 = file.path(tmp.fa,paste0(g1,"_",g2,".blat.fa"))
  paste("blat", fg1, fc2, "-t=dna -q=dna -out=blast8 -maxIntron=20000 -minScore=100", blat1)
  d = fread(blat1)
  d[, prop := V12/max(V12),
    by = list(V1)]
  d <- d[d$prop>=.5,]
  b.start = tapply(apply(d[,9:10],1,min),d$V1,min)-blk.buffer
  b.end = tapply(apply(d[,9:10],1,max),d$V1,min)+blk.buffer
  b.start[b.start<0]<-0
  b.end[b.end>width(wg1)]<-width(wg1)
  b.strand = tapply(ifelse(d$V10 > d$V9,"+","-"),d$V1, function(x) names(sort(table(x),decreasing=TRUE))[1])
  b.strand =
  bo1 = data.frame(chr = i,
                   start = b.start
                   end = b.end,
                   id = paste0(g1,"_",names(b.start)),
                   score = "1",
                   strand = b.strand)
  write.table(bo1, file=bed1, quote=F, sep="\t", row.names=F, col.names=F)

  system(paste("bedtools getfasta -fi", ))


  blen = length()

  tapply(d$start, d$V1, min)- tapply(d$end, d$V1, max)
  --query q1.fa q2.fa t2.fa --target t1.fa t2.fa t3.fa
  system(paste("exonerate --model protein2genome",
               "--forcescan target --softmasktarget",
               "--gappedextension --refine region",
               "--dpmemory 2000 --maxintron 10000",
               "--showalignment no --showtargetgff yes -n 1",
               "--query",f2, "--target", fg1))
})


cat ~/Downloads/test.gff | sed '/^Command/ d' | sed '/^Hostname/ d'
exonerate --model protein2genome --forcescan target --softmasktarget --gappedextension --refine region --dpmemory 2000 --maxintron 10000 --showvulgar no --showalignment no --showtargetgff yes -n 1 --query /Users/jlovell/Documents/test_compareGenespace/tmp/single.gene.fasta/1/Phallii.fa --target /Users/jlovell/Documents/test_compareGenespace/tmp/single.gene.fasta/1/PhalliiHAL.blk.fa | sed '/^#/ d' > ~/Downloads/test.gff

plot_mapping(blk = merged$block, map = merged$map)
good.blocks = ogs.bygroup$proptype
good.blocks = good.blocks$block.id[good.blocks$p.2>=.1]
good.map = merged$map[merged$map$block.id %in% good.blocks,]
good.block = merged$block[merged$block$block.id %in% good.blocks,]
plot_mapping(blk = good.block, map = good.map, genomes = c("PhalliiHAL","Sviridis"))





