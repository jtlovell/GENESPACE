nomap.list = lapply(names(of.inblk), function(i){
  cat(i,"\n")
  out = find_missingHitRegion(gff = gff, inblk.of.out = of.inblk[[i]])
  return(out)
})
nomap<-rbindlist(nomap.list)

cds.fastas = sapply(genomeIDs, simplify = F, USE.NAMES = T,function(i){
  Biostrings::readDNAStringSet(file.path(cds.dir,paste0(i,".fa")))
})

exonerate.list = lapply(1:nrow(nomap), function(i){
  cat("seq",paste0(i,"/",nrow(nomap)),"... target length = ")
  bed = nomap[i,c("chr","start","end","id"),with = F]
  genome = nomap$query.genome[i]
  cds.genome = nomap$target.genome[i]
  target.gene = nomap$gene[i]
  cat(width(cds.fastas[[cds.genome]][target.gene]))
  out = run_exonerate(assembly.dir = assembly.dir,
                      target.gene = target.gene,
                      cds.fastas = cds.fastas,
                      genome = genome,
                      cds.genome = cds.genome,
                      bed = bed,
                      tmp.dir = tmp.dir)
  cat(" align width =",width(out$cds.seq),"\n")
  return(out)
})
nomap$align.width = sapply(exonerate.list, function(x)  width(x$cds.seq))
nomap$cds.width = apply(nomap,1, function(x)  width(cds.fastas[[x[6]]][x[7]]))
nomap$align.QC = with(nomap, align.width <= (cds.width*1.5) &
                        align.width >= (cds.width*.05))





test  =find_missingHitRegion(gff = gff, inblk.of.out = of.inblk[[1]])



}

blast.inblock = lapply(list.files(base_dir,full.names = T), function(x){
  blast.dir = file.path(x,"blast")

})

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

  cds = cds.fastas[[cds.genome]][target.gene]
  Biostrings::writeXStringSet(cds, filepath = cdsf)

  out1 = system(paste("exonerate --model est2genome",
                      "--forcescan target --softmasktarget",
                      "--gappedextension --refine region",
                      "--dpmemory 2000 --maxintron 10000 -n 1 ",
                      "--showalignment no --showtargetgff yes --showvulgar no --showcigar no",
                      '--ryo ">\n%tas\n"',
                      "--query",cdsf, "--target", fab), intern = T)

  if(length(grep("^>",out1))==0){
    out1 = DNAStringSet("")
    names(out1)<-paste(bed$id, target.gene, sep = "-")
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
    names(out1)<-paste(bed$id, target.gene, sep = "-")
  }
  return(list(cds.seq = out1, gff = gffo))
}



find_missingHitRegion = function(inblk.of.out,
                                 gff,
                                 no.bound.buffer = 1e4,
                                 bound.buffer = 1e3){
  x = inblk.of.out

  only1 = x$orthogroups[x$orthogroups$og %in% x$no.ortho,]
  if(nrow(only1)==0){
    return(NULL)
  }else{
    cols = c("genome","id","chr","start","end","order")
    cols1 = paste0(cols, 1)
    cols2 = paste0(cols,2)
    both1 = data.frame(x$cull.blast[,c(cols1,"score"),with =F])
    both2 = data.frame(x$cull.blast[,c(cols2,"score"),with =F])
    colnames(both1)<-c(cols,"score")
    colnames(both2)<-c(cols,"score")
    both1$pair = 1:nrow(both1)
    both2$pair = 1:nrow(both2)
    both = rbind(both1, both2)

    out = rbindlist(lapply(1:nrow(only1), function(i){
      x1 = only1[i,]

      gff.both = gff[gff$genome == x1$genome & gff$id %in% both$id,]
      gff.x = gff[gff$genome == x1$genome & gff$id == x1$id,]
      gff.both$diff.x = gff.both$order - gff.x$order
      gff.both$left = gff.both$order < gff.x$order
      gff.both$right = gff.both$order > gff.x$order
      gff.left = gff.both[gff.both$left,]
      gene.left = gff.left$id[which.min(abs(gff.left$diff.x))]
      gff.right = gff.both[gff.both$right,]
      gene.right = gff.right$id[which.min(abs(gff.right$diff.x))]

      if(length(gene.left)==1 & length(gene.right)==1){
        pairs = unique(both$pair[both$id %in% c(gene.right, gene.left) & both$genome == x1$genome])
        tmp = both[both$pair %in% pairs,]
        good.pairs = sapply(split(tmp[tmp$genome ==  x1$genome,],
                                  tmp$id[tmp$genome ==  x1$genome]), function(y){
                                    y$pair[which.max(y$score)]
                                  })
        tmp = tmp[tmp$pair %in% good.pairs,]
        coords = tmp[tmp$genome != x1$genome,]
        coords = data.frame(chr = coords$chr[1],
                            start = min(unlist(coords[,c("start","end")]))-bound.buffer,
                            end = max(unlist(coords[,c("start","end")]))+bound.buffer,
                            id = paste(coords$id, collapse = "_"),
                            query.genome = coords$genome[1],
                            target.genome = x1$genome,
                            gene = x1$id,
                            stringsAsFactors = F)
      }else{
        if(length(gene.left)==0 & length(gene.right)==0){
          coords = both[both$genome != x1$genome,]
          coords = data.frame(chr = coords$chr[1],
                              start = min(unlist(coords[,c("start","end")]))-no.bound.buffer,
                              end = max(unlist(coords[,c("start","end")]))+no.bound.buffer,
                              id = paste(coords$id, collapse = "_"),
                              query.genome = coords$genome[1],
                              target.genome = x1$genome,
                              gene = x1$id,
                              stringsAsFactors = F)
        }else{
          if(length(gene.left)==0){
            gene.left  = gff.both$id[gff.both$genome == x1$genome & gff.both$start == min(gff.both$start[gff.both$genome == x1$genome])]
            pairs = unique(both$pair[both$id %in% c(gene.right, gene.left) & both$genome == x1$genome])
            tmp = both[both$pair %in% pairs,]
            good.pairs = sapply(split(tmp[tmp$genome ==  x1$genome,],
                                      tmp$id[tmp$genome ==  x1$genome]), function(y){
                                        y$pair[which.max(y$score)]
                                      })
            tmp = tmp[tmp$pair %in% good.pairs,]
            coords = tmp[tmp$genome != x1$genome,]
            coords = data.frame(chr = coords$chr[1],
                                start = min(unlist(coords[,c("start","end")]))-no.bound.buffer,
                                end = max(unlist(coords[,c("start","end")]))+bound.buffer,
                                id = paste(coords$id, collapse = "_"),
                                query.genome = coords$genome[1],
                                target.genome = x1$genome,
                                gene = x1$id,
                                stringsAsFactors = F)
          }else{
            gene.right = gff.both$id[gff.both$genome == x1$genome & gff.both$end == max(gff.both$end[gff.both$genome == x1$genome])]
            pairs = unique(both$pair[both$id %in% c(gene.right, gene.left) & both$genome == x1$genome])
            tmp = both[both$pair %in% pairs,]
            good.pairs = sapply(split(tmp[tmp$genome ==  x1$genome,],
                                      tmp$id[tmp$genome ==  x1$genome]), function(y){
                                        y$pair[which.max(y$score)]
                                      })
            tmp = tmp[tmp$pair %in% good.pairs,]
            coords = tmp[tmp$genome != x1$genome,]
            coords = data.frame(chr = coords$chr[1],
                                start = min(unlist(coords[,c("start","end")]))-bound.buffer,
                                end = max(unlist(coords[,c("start","end")]))+no.bound.buffer,
                                id = paste(coords$id, collapse = "_"),
                                query.genome = coords$genome[1],
                                target.genome = x1$genome,
                                gene = x1$id,
                                stringsAsFactors = F)
          }
        }
      }
      return(coords)
    }))
  }
  return(out)
}
