
nothing = chop_assemblyByBlock(blk = blk,
                               out.dir = tmp.dir,
                               assembly.dir = assembly.dir,
                               verbose = T)

# set up run:
blk$unique = with(blk, paste(genome1, genome2, block.id, sep = "_"))
x = blk[i,]

cds.fastas = do.call(c,lapply(genomeIDs, function(j)
  Biostrings::readAAStringSet(file.path(dirs$cds, paste0(j,".fa")))))
test = mclapply(c(1,2), mc.cores = 1, mc.preschedule = F, function(i){
  print(i)
  return(pipe_blast2exonerate(x=blk[i,]))
})

pipe_blast2exonerate<-function(x){
  wd = file.path(tmp.dir,x$unique)
  if(!dir.exists(wd)){
    dir.create(wd)
  }

  chop_assemblyByBlock(out.dir = wd,
                       x = x,
                       assembly.dir = assembly.dir)

  orphans = find_unMappedGenes(gff = gff,
                               map = map,
                               blk = x)[[1]]
  orphans[[1]]<-orphans[[1]][orphans[[1]] %in% hiqual.genes]
  orphans[[2]]<-orphans[[2]][orphans[[2]] %in% hiqual.genes]

  chop_unMapPeptideByBlock(orphans = orphans,
                           peptide.dir = peptide.dir,
                           tmp.dir = wd)

  pep1.fa = file.path(wd,"pep1.fa")
  pep2.fa = file.path(wd,"pep2.fa")
  fa1.fa = file.path(wd,"seq1.fa")
  fa2.fa = file.path(wd,"seq2.fa")

  blast12 = file.path(wd,"x1v2.blast")
  blast21 = file.path(wd, "x2v1.blast")

  dmnd1 = file.path(wd,"db1")
  dmnd2 = file.path(wd,"db2")

  system(paste("diamond makedb --quiet --in", pep1.fa, "-d",dmnd1))
  system(paste("diamond blastx --quiet --sensitive",
               "--threads",1,
               other.blast.param,
               "-d", dmnd1,
               "-q", fa2.fa,
               "-o",blast12))

  system(paste("diamond makedb --quiet --in", pep2.fa, "-d",dmnd2))
  system(paste("diamond blastx --quiet --sensitive",
               "--threads",1,
               other.blast.param,
               "-d", dmnd2,
               "-q", fa1.fa,
               "-o",blast21))

  bl12 = fread(blast12)
  bl21 = fread(blast21)

  hit.reg = get_hitRegion(bl = rbind(bl12,bl21),
                          assembly.dir = assembly.dir)
  hit.reg$reg.name = with(hit.reg, paste(genome,start,end,id,sep = "_"))

  exon.out = pipe_exonerate(hit.reg = hit.reg,
                            assembly.dir = dirs$assembly,
                            cds.fastas = cds.fastas,
                            genomeIDs = genomeIDs,
                            tmp.dir = wd,
                            verbose = F)
  return(exon.out)
}


get_hitRegion = function(bl,
                         buffer = 1e3,
                         assembly.dir,
                         max.dist2besthit = 5e3){
  genome.info = lapply(bl$V1, function(x) strsplit(x,".", fixed = T)[[1]])
  bl$genome = sapply(genome.info, function(x) x[1])
  bl$chr = sapply(genome.info, function(x) x[2])
  bl$genome.start = as.numeric(sapply(genome.info, function(x) x[3]))
  bl$genome.end = bl$genome.start+bl$V8
  bl$genome.start = bl$genome.start+bl$V7
  bl$length = bl$genome.end - bl$genome.start
  bl$start = with(bl, ifelse(length < 0, genome.end, genome.start))
  bl$end = with(bl, ifelse(length > 0, genome.end, genome.start))
  bl$length = bl$end - bl$start
  setnames(bl, "V2", "id")
  bl$neg.score = bl$V12 * (-1)
  bl$neg.len = bl$length * (-1)
  setkey(bl, id, neg.score, neg.len)

  bl[,rank := frank(neg.score, ties.method = "random"),
     by = list(id)]
  bs = bl[bl$rank == 1,c("id","start","end")]
  setnames(bs, c("id","best.start","best.end"))
  setkey(bl,id)
  setkey(bs,id)
  bo = merge(bl, bs)
  bo$dist2best = with(bo, abs(best.start - start))
  bo = bo[bo$dist2best <= max.dist2besthit,]



  bi = bo[,list(start = min(start),
                end = max(end)),
          by = list(genome, chr, id)]
  bi$length = with(bi, end-start)

  if(buffer != 0){
    fais <- rbindlist(lapply(unique(bi$genome), function(i){
      tmp = read.delim(file.path(assembly.dir, paste0(i,".fa.fai")),
                       header = F, stringsAsFactors = F,
                       col.names = c("chr",  "chr.length", "v1", "v2", "v3"))[, 1:2]
      tmp$genome = i
      return(tmp)
    }))
    setkey(fais, genome, chr)
    bo = merge(fais, bi)
    bo$start <- bo$start - buffer
    bo$end <- bo$end + buffer
    bo$end[bo$end>bo$chr.length]<-bo$chr.length[bo$end>bo$chr.length]
    bo$start[bo$start<=0]<-1
    bi <- bo
  }

  return(bi)
}

chop_unMapPeptideByBlock = function(orphans,
                                    peptide.dir,
                                    tmp.dir){

  fastas = do.call(c,lapply(list.files(peptide.dir, full.names = T), function(j)
    Biostrings::readAAStringSet(j)))

  fas1 = fastas[orphans[[1]]]
  fas2 = fastas[orphans[[2]]]

  p1 = file.path(tmp.dir, "pep1.fa")
  p2 = file.path(tmp.dir,"pep2.fa")
  Biostrings::writeXStringSet(fas1, filepath = p1)
  Biostrings::writeXStringSet(fas2, filepath = p2)

}


run_blastxFromGeneList = function(tmp.dir,
                                  verbose = T,
                                  max.threads = 1,
                                  other.blast.param = "--matrix BLOSUM62 --frameshift 15"){

  pep1.fa = file.path(tmp.dir,"pep1.fa")
  pep2.fa = file.path(tmp.dir,"pep2.fa")
  fa1.fa = file.path(tmp.dir,"seq1.fa")
  fa2.fa = file.path(tmp.dir,"seq2.fa")

  blast12 = file.path(tmp.dir,paste0(i,".1v2.blast"))
  blast21 = file.path(tmp.dir,paste0(i,".2v1.blast"))

  dmnd1 = file.path(tmp.dir,paste0(i,".db1"))
  dmnd2 = file.path(tmp.dir,paste0(i,".db2"))

  system(paste("diamond makedb --quiet --in", pep1.fa, "-d",dmnd1))
  system(paste("diamond blastx --quiet --sensitive",
               "--threads",max.threads,
               other.blast.param,
               "-d", dmnd1,
               "-q", fa2.fa,
               "-o",blast12))

  system(paste("diamond makedb --quiet --in", pep2.fa, "-d",dmnd2))
  system(paste("diamond blastx --quiet --sensitive",
               "--threads",max.threads,
               other.blast.param,
               "-d", dmnd2,
               "-q", fa1.fa,
               "-o",blast21))

  bl12 = fread(blast12)
  bl21 = fread(blast21)
  if(verbose)
    cat("Done!\n\t")
  out[[i]]$blast_1p2g = bl12
  out[[i]]$blast_2p1g = bl21
  return(out)
}





make_blockDB = function(tmp.dir,
                        assembly.dir,
                        peptide.dir,
                        only.these.genes = NULL,
                        blk,
                        map,
                        gff,
                        buffer = 1000,
                        max.blast.threads = 1,
                        other.blast.param = "--sensitive --matrix BLOSUM62 --frameshift 15",
                        verbose = T){

  nothing = chop_assemblyByBlock(blk = blk,
                                 out.dir = tmp.dir,
                                 assembly.dir = assembly.dir,
                                 verbose = T)

  if(verbose)
    cat("Done!\nFinding unmapped genes ...")
  gene.list = find_unMappedGenes(gff = gff,
                                 only.these.genes = only.these.genes,
                                 map = map,
                                 blk = blk)

  if(verbose)
    cat("Done!\nChopping peptide fastas by genes in blocks ...")
  nothing = chop_unMapPeptideByBlock(gene.list = gene.list,
                                     peptide.dir = peptide.dir,
                                     tmp.dir = tmp.dir)

  if(verbose)
    cat("Done!\n")

  blastx.out2 <- run_blastx(blk,
                            gene.list = gene.list[2],
                                        verbose = T,
                                        max.threads = 8,
                                        other.blast.param = "--min-score 10")

  if(verbose)
    cat("Found",sum(sapply(unlist(blastx.out, recursive = F),nrow)),"hits\n")

  if(verbose)
    cat("Pulling regions with best BLASTX hit ...")
  hit.reg = get_hitRegion(blast.output.list = blastx.out[2],
                          assembly.dir = dirs$assembly,
                          buffer = buffer)
  if(verbose)
    cat("Done!\n")
  exon.out = pipe_exonerate(hit.reg = hit.reg,
                            assembly.dir = dirs$assembly,
                            cds.dir = dirs$cds, genomeIDs = genomeIDs,
                            tmp.dir = dirs$tmp.dir,
                            verbose = T)
  if(verbose)
    cat("Completed!\n")

  return(list(exonerate.out = exon.out,
              hit.regions = hit.reg,
              blast.out = blastx.out,
              unmapped.genes = gene.list))

}
blk = res.all$rerun.results$block
map = res.all$rerun.results$map
gff = res.all$initial.results$gff
assembly.dir = dirs$assembly
peptide.dir = dirs$peptide
gff.dir = dirs$gff
cds.dir = dirs$cds
tmp.dir = dirs$tmp.dir
nothing = chop_assemblyByBlock(blk = blk,
                            out.dir = dirs$tmp.dir,
                            assembly.dir = dirs$assembly)







tmp.dir = dirs$tmp.dir

gene.qual = import_ovlps(ovlp.dir = file.path(directory,"ovlps"),
                         genomeIDs = genomeIDs)

hiqual.genes = gene.qual$id[gene.qual$confidence == "high"]
modqual.genes = gene.qual$id[gene.qual$confidence == "moderate"]
loqual.genes = gene.qual$id[gene.qual$confidence == "low"]

test = unlist(gene.list)
table(hiqual.genes %in% test) / length(hiqual.genes)
table(modqual.genes %in% test) / length(modqual.genes)
table(loqual.genes %in% test) / length(loqual.genes)

test2 = hit.reg$id

table(test2 %in% hiqual.genes) / length(test2)
table(test2 %in% modqual.genes) / length(test2)
table(test2 %in% loqual.genes) / length(test2)

blastx.out = run_blastxFromGeneList(gene.list = gene.list)
hit.reg = get_hitRegion(blast.output.list = blastx.out, assembly.dir = dirs$assembly,
                        buffer = 1000)
exon.out = pipe_exonerate(hit.reg = hit.reg,
                          assembly.dir = dirs$assembly,
                          cds.dir = dirs$cds, genomeIDs = genomeIDs,
                          tmp.dir = dirs$tmp.dir)
save(exon.out, file = "exon.tmp.rda")
pipe_exonerate <- function(hit.reg,
                           assembly.dir,
                           cds.fastas,
                           genomeIDs,
                           tmp.dir,
                           verbose = T){

  if(verbose)
    cat("Importing CDS fastas\n")

  if(verbose)
    cat("Running",nrow(hit.reg),"exonerate searches\n\tCompleted: ")

  exon.out = lapply(1:nrow(hit.reg), function(i){
    if(verbose)
      if(i %% 100 == 0)
        cat(paste0(i,", "))
    target.gene = hit.reg$id[i]
    genome = hit.reg$genome[i]
    bed = hit.reg[i,c("chr","start","end","reg.name"),with = F]
    colnames(bed)[4]<-"id"
    tmp <- run_exonerate(assembly.dir = dirs$assembly,
                        cds.fastas = cds.fastas,
                        target.gene = target.gene,
                        genome = genome,
                        cds.genome = cds.genome,
                        bed = bed,
                        tmp.dir = tmp.dir)
    return(tmp)
  })

  exon_cds <- do.call(c, lapply(exon.out, function(x) x$cds.seq))
  gff.gene <- rbindlist(lapply(exon.out, function(x) x$gff))
  gff.gene$ex.start = with(gff.gene, as.numeric(start)+as.numeric(exon.start))
  gff.gene$ex.end = with(gff.gene, as.numeric(start)+as.numeric(exon.end))
  gff.out <- gff.gene[,list(start = min(ex.start),
                            end = max(ex.end)),
                      by = list(id, chr, strand)]

  if(verbose)
    cat("\nDone!\n")

  return(list(gff = gff.out,
              exonerate_cds = exon_cds))

}


i = 1






  if(any(bi$length > max.hit.length)){
    hilen = bi$id[bi$length > max.hit.length]
    hitmp = bl[bl$id %in% hilen,]
    hitmp = hitmp[!duplicated(hitmp$id),names(bi),with = F]
    bi = rbind(bi[!bi$id %in% hilen,], hitmp)
  }
  setkey(bi, genome, chr)

  if(buffer != 0){
    fais <- rbindlist(lapply(unique(bi$genome), function(i){
      tmp = read.delim(file.path(assembly.dir, paste0(i,".fa.fai")),
                 header = F, stringsAsFactors = F,
                 col.names = c("chr",  "chr.length", "v1", "v2", "v3"))[, 1:2]
      tmp$genome = i
      return(tmp)
    }))
    setkey(fais, genome, chr)
    bo = merge(fais, bi)
    bo$start <- bo$start - buffer
    bo$end <- bo$end + buffer
    bo$end[bo$end>bo$chr.length]<-bo$chr.length[bo$end>bo$chr.length]
    bo$start[bo$start<=0]<-1
    bi <- bo
  }
  bi$pep.genome = sapply(bi$block.id, function(x) strsplit(x,"_")[[1]][2])
  bi$reg.name = with(bi, paste(genome,id,sep = "_"))

  return(bi)
}


chop_assemblyByBlock = function(x,
                                out.dir,
                                assembly.dir,
                                verbose = T){
  f = out.dir

  g1 = x$genome1
  g2 = x$genome2
  id = paste(g1,g2,x$block.id,sep = "_")
  bed1 = x[,c("chr1","start1","end1"),with = F]
  bed1$unique = with(bed1, paste(g1,chr1,start1,end1, sep = "."))
  bed2 = x[,c("chr2","start2","end2"), with = F]
  bed2$unique = with(bed2, paste(g2,chr2,start2,end2, sep = "."))

  bed1f = file.path(f,"bed1.bed")
  bed2f = file.path(f,"bed2.bed")
  fa1 = file.path(assembly.dir, paste0(g1,".fa"))
  fa2 = file.path(assembly.dir, paste0(g2,".fa"))
  fafo1 = file.path(f,"seq1.fa")
  fafo2 = file.path(f,"seq2.fa")


  write.table(bed1, file=bed1f,
              quote=F, sep="\t", row.names=F, col.names=F)
  write.table(bed2, file=bed2f,
              quote=F, sep="\t", row.names=F, col.names=F)
  system(paste("bedtools getfasta -fi",
               fa1, "-bed", bed1f,"-name -fo",fafo1))
  system(paste("bedtools getfasta -fi",
               fa2, "-bed", bed2f,"-name -fo",fafo2))
}


import_ovlps = function(ovlp.dir, genomeIDs){

  ret = rbindlist(lapply(genomeIDs, function(i){
    o = fread(file.path(ovlp.dir, paste0(i,".ovlps")))
    if(ncol(o)==18){
      o$V12<-NULL
    }
    setnames(o,c("trs.id", "EST.list","intron.valid.cdsratio", "cscore",
                 "homology.cov", "gene.completescore", "cdslength", "cdscount",
                 "cdsratio", "EST.CDS.overall.overlap.ratio",
                 "hitName", "score", "gapRatio", "peptideMotif",
                 "TEdomainCount", "domainCount", "expressionLevel"))
    o$id = substr(o$trs.id,1, nchar(o$trs.id)-2)
    o$genome = i
    o$final_score = with(o, (expressionLevel > 1) + homology.cov + cscore)
    o$ignoreExpr_score = with(o, homology.cov + cscore)

    o$confidence = with(o, ifelse(final_score <=1, "low",
                                  ifelse(final_score <= 2, "moderate",
                                         ifelse(final_score > 2, "high", NA))))

    o$ignoreExpr_confidence = with(o, ifelse(ignoreExpr_score <.5, "low",
                                  ifelse(ignoreExpr_score < 2, "moderate",
                                         ifelse(ignoreExpr_score >= 2, "high", NA))))

    out = o[,c("genome","id","expressionLevel","homology.cov","intron.valid.cdsratio","cscore","final_score","confidence","ignoreExpr_confidence")]
    return(out)
  }))

  return(ret)
}


find_unMappedGenes = function(gff,
                              map,
                              blk,
                              only.these.genes = NULL,
                              ovlp.dir = NULL){
  blk$unique = with(blk, paste(genome1,genome2,block.id,sep = "_"))
  map$unique = with(map, paste(genome1,genome2,block.id,sep = "_"))
  ol = sapply(blk$unique, USE.NAMES = T, simplify = F, function(i){
    tm = map[map$unique == i,]
    tb = blk[blk$unique == i,]

    g1 = gff[with(gff, genome == tb$genome1 &
                    chr == tb$chr1 &
                    start <= tb$end1 &
                    end >= tb$start1),]

    g2 = gff[with(gff, genome == tb$genome2 &
                    chr == tb$chr2 &
                    start <= tb$end2 &
                    end >= tb$start2),]


    go1 = g1$id[!g1$id %in% tm$id1]
    go2 = g2$id[!g2$id %in% tm$id2]

    if(!is.null(only.these.genes)){
      go1 = go1[go1 %in% only.these.genes]
      go2 = go2[go2 %in% only.these.genes]
    }

    return(list(genome1 = go1, genome2 = go2))
  })

  return(ol)
}


chop_unMapPeptideByBlock = function(gene.list,
                                    peptide.dir,
                                    tmp.dir){

  fastas = do.call(c,lapply(list.files(peptide.dir, full.names = T), function(j)
    Biostrings::readAAStringSet(j)))

  lapply(names(gene.list), function(i){
    fas1 = fastas[gene.list[[i]][[1]]]
    fas2 = fastas[gene.list[[i]][[2]]]

    p1 = file.path(tmp.dir,paste0(i,".pep1.fa"))
    p2 = file.path(tmp.dir,paste0(i,".pep2.fa"))
    Biostrings::writeXStringSet(fas1, filepath = p1)
    Biostrings::writeXStringSet(fas2, filepath = p2)
  })
}

run_exonerate <- function(assembly.dir,
                          cds.fastas,
                          target.gene,
                          genome,
                          cds.genome,
                          bed,
                          tmp.dir){
  faf = file.path(assembly.dir,paste0(genome,".fa"))
  bedf = file.path(tmp.dir, "tmp.bed")
  fab = file.path(tmp.dir,"tmp.fa")
  cdsf = file.path(tmp.dir,"tmp.cds.fa")
  write.table(bed, file=bedf,
              quote=F, sep="\t", row.names=F, col.names=F)
  system(paste("bedtools getfasta -fi",
               faf, "-bed", bedf,"-name -s -fo",fab))

  cds = cds.fastas[target.gene]
  Biostrings::writeXStringSet(cds, filepath = cdsf)

  out1 = system(paste("exonerate --model est2genome",
                      # "--forcescan target --softmasktarget",
                      # "--gappedextension --refine region",
                      "--bestn 1",
                      "--dpmemory 2000 --maxintron 10000",
                      "--showalignment no --showtargetgff yes --showvulgar no --showcigar no",
                      '--ryo ">\n%tas\n"',
                      "--query",cdsf, "--target", fab), intern = T)

  if(length(grep("^>",out1))==0){
    out1 = DNAStringSet("")
    names(out1)<-bed$id
    gff1 = data.frame(exon.start = NA,
                      exon.end = NA,
                      strand = NA,
                      align.info = NA)
    gffo = data.frame(bed,gff1)
  }else{
    end = which(out1 == "")[1]+1
    out1 = out1[1:end]
    gff1 = do.call(rbind,lapply(out1[grepl(bed$id, out1) & grepl("exon\t", out1)],
                                function(x) strsplit(x,"\t")[[1]][c(4:5,7,9)]))
    gff1 = data.frame(gff1,
                      stringsAsFactors = F)
    colnames(gff1)<-c("exon.start","exon.end","strand","align.info")
    gffo = data.frame(bed,gff1)

    wh = grep("^>",out1)+1
    out1 = DNAStringSet(paste(out1[wh:(length(out1)-2)], collapse = ""))
    names(out1)<-bed$id
  }
  return(list(cds.seq = out1, gff = gffo))
}




