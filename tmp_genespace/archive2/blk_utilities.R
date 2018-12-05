blk = make_blockDir(blk = merged.close$block, block.dir = block.dir)
sgff = split_gffByBlock(gff = gff, blk = blk)
# chop_assemblyByBlock(blk = blk, assembly.dir = assembly.dir)
chop_annotationByBlock(sgff = sgff, blk = blk)
of_inblock(blk =  blk)

od = og.inblk$orthogroup.datatable
od2 = od[od$og.n.genomes > 1,]
odsplit = split.data.table(od2, by = "unique")
odg = rbindlist(lapply(odsplit, function(x) expand.grid(x$id, x$id)))
g = graph_from_data_frame(data.frame(odg))
cl = clusters(g)
mem = cl$membership


og.inblk = read_orthoGroups(block.dir = block.dir,
                            gff = gff,
                            ncores = 6)
concatenate_orthoGroups = function(og.inblk, n.iter = 1, verbose){
  og = og.inblk$orthogroup.datatable
  og$unique = paste0(og$block.id, "_", og$og)
  spl.gen = split(og, og$genome)
  if(verbose)
    cat("Building species-specific databases ... completed:")
  spl.out = rbindlist(lapply(names(spl.gen), function(i){
    if(verbose)
      cat(paste0(i,", "))
    print(i)
    x = spl.gen[[i]]
    spl.id = split(x, x$id)
    og.uniq = lapply(spl.id, function(y) unique(y$unique))
    og.uout = data.table(id = rep(names(spl.id),sapply(og.uniq,length)), unique = unlist(og.uniq))
    og.uout$genome = i
    return(og.uout)
  }))

  spl.gen = split(og, og$block.id)
  spl.out2 = rbindlist(mclapply(names(spl.gen), mc.cores = 6, function(i){
    if(verbose)
      cat(paste0(i,", "))
    x = spl.gen[[i]]
    spl.id = split(x, x$og)
    gene.comb = rbindlist(lapply(spl.id, function(y) expand.grid(y$id, y$id)))
    return(gene.comb)
  }))
  splt = split(spl.out, spl.out$unique)
  outl = list()
  outu = vector()
  ns = 0
  u = unique(spl.out$id)
  for(i in u){
    if(which(u == i) %% 50 == 0) cat("\tCompleted",which(u == i),"/", length(u),"\n")
    if(!i %in% outu){
      ns = ns+1
      t1 = spl.out[spl.out$id == i,]
      t2 = spl.out[spl.out$unique %in% t1$unique,]
      t3 = spl.out[spl.out$id %in% t2$id,]

      t4 = spl.out[spl.out$unique %in% t3$unique,]
      t5 = spl.out[spl.out$id %in% t4$id,]
      out = unique(t5$id)
      outl[[i]] = out
      outu = c(outu, out)
    }
  }


  graph_from_data_frame(d, directed = TRUE, vertices = NULL)

  if(verbose)
    cat("Connecting orthogroups across pairwise comparisons ... completed\n")

}
make_blockDir = function(blk,
                         block.dir,
                         verbose = T){s
  if(file.exists(block.dir)){
    system(paste("rm -rf", block.dir))
  }
  system(paste("mkdir", block.dir))

  if(verbose)
    cat("Building directories for", nrow(blk), "blocks\n")
  blk$unique = with(blk, paste(genome1, genome2, block.id, sep = "_"))
  blk$dir = sapply(1:nrow(blk), function(i){
    x = file.path(block.dir,blk$unique[i])
    system(paste("mkdir",x))
    system(paste("mkdir",file.path(x,"genome")))
    for(j in c("assembly","cds", "gff", "peptide", "transcript")){
      system(paste("mkdir",file.path(x,"genome",j)))
    }
    for(j in c("tmp","blast")){
      system(paste("mkdir",file.path(x,j)))
    }
    return(x)
  })
  if(verbose)
    cat("Done!\n")
  return(blk)
}

split_gffByBlock = function(blk,
                          gff,
                          verbose = T){
  gff$unique = with(gff, paste(genome, chr, sep = "_"))
  blk$unique1 = with(blk, paste(genome1, chr1, sep = "_"))
  blk$unique2 = with(blk, paste(genome2, chr2, sep = "_"))

  if(verbose)
    cat("Splitting gff into", nrow(blk),"blocks ... ")
  sgff = split.data.table(gff, "unique")
  gffo <- lapply(1:nrow(blk), function(i){
    btmp = blk[i,]
    g1 = sgff[[btmp$unique1]]
    g2 = sgff[[btmp$unique2]]
    g1 = g1[g1$start<=btmp$end1 & g1$end>=btmp$start1,]
    g2 = g2[g2$start<=btmp$end2 & g2$end>=btmp$start2,]
    go = rbind(g1, g2)
    # go$blk = btmp$unique
    return(go)
  })

  names(gffo)<-with(blk, paste(genome1, genome2, block.id))
  if(verbose)
    cat("Done!\n")
  return(gffo)
}

chop_assemblyByBlock = function(blk,
                                assembly.dir,
                                verbose = T){

  if(verbose)
    cat("Chopping up assembly fasta by block breakpoints ... \n")
  for(i in 1:nrow(blk)){
    if(verbose)
      if(i %% 100 == 0) cat("\tCompleted:",i,"/", nrow(blk),"\n")
    x = blk[i,]
    g1 = x$genome1
    g2 = x$genome2
    f = x$dir
    bed1 = x[,c("chr1","start1","end1"),with = F]
    bed1$unique = with(bed1, paste(g1,chr1,start1,end1, sep = "."))
    bed2 = x[,c("chr2","start2","end2"), with = F]
    bed2$unique = with(bed2, paste(g2,chr2,start2,end2, sep = "."))

    bed1f = file.path(f,"genome","bed1.bed")
    bed2f = file.path(f,"genome","bed2.bed")
    fa1 = file.path(assembly.dir,paste0(g1,".fa"))
    fa2 = file.path(assembly.dir,paste0(g2,".fa"))
    fafo1 = file.path(f,"genome","assembly","seq1.fa")
    fafo2 = file.path(f,"genome","assembly","seq2.fa")


    write.table(bed1, file=bed1f,
                quote=F, sep="\t", row.names=F, col.names=F)
    write.table(bed2, file=bed2f,
                quote=F, sep="\t", row.names=F, col.names=F)
    system(paste("bedtools getfasta -fi",
                 fa1, "-bed", bed1f,"-name -fo",fafo1))
    system(paste("bedtools getfasta -fi",
                 fa2, "-bed", bed2f,"-name -fo",fafo2))
  }
  if(verbose)
    cat("Done!\n")
}

chop_annotationByBlock = function(sgff,
                                  blk){


  flist = list(peptide = peptide.dir,
               cds = cds.dir,
               transcript = transcript.dir)

  for(i in names(flist)){
    if(verbose)
      cat(i,"fastas ... Importing ...")

    fastas = lapply(list.files(flist[[i]], full.names = T), function(j){
      if(i == "peptide"){
        Biostrings::readAAStringSet(j)
      }else{
        Biostrings::readDNAStringSet(j)
      }
    })
    names(fastas)<-gsub(".fa","",list.files(flist[[i]], full.names = F),fixed = F)

    if(verbose)
      cat("Parsing and writing ...")
    test = lapply(1:nrow(blk), function(j){
      x = blk[j,]
      g = sgff[[x$unique]]
      sg = split(g$id, g$genome)
      t1 = fastas[[names(sg)[1]]][sg[[1]]]
      t2 = fastas[[names(sg)[2]]][sg[[2]]]
      p1 = file.path(x$dir,"genome",i,"seq1.fa")
      p2 = file.path(x$dir,"genome",i,"seq2.fa")
      Biostrings::writeXStringSet(t1, filepath = p1)
      Biostrings::writeXStringSet(t2, filepath = p2)
    })

    if(verbose)
      cat("Done\n")
  }
}


of_inblock = function(blk, ncores = 6){
  if(verbose)
    cat("Running orthofinder in blocks\n")
  test = mclapply(1:nrow(blk), mc.cores = ncores, mc.preschedule = F, function(i){
    x = blk[i,]
    run_orthofinder(
      peptide.dir = file.path(x$dir, "genome", "peptide"),
      og.threads = 1,
      blast.threads = 1,
      og.silent = T,
      verbose = F,
      tmp.dir = file.path(x$dir, "tmp"),
      blast.dir = file.path(x$dir, "blast"))
  })
  if(verbose)
    cat("Done!\n")
}

read_orthoGroups = function(block.dir,
                            gff,
                            ncores = 6,
                            verbose = T){
  if(verbose)
    cat("Looking for Orthogroup text files\n")

  ogf = list.files(block.dir, pattern = "Orthogroups.txt", full.names = T, recursive = T)
  bn <- basename(ogf)
  ogf<- ogf[bn == "Orthogroups.txt" & grepl("blast/Orthogroups.txt", ogf, fixed = T)]
  dirs = gsub("/","",gsub("/blast","",gsub(block.dir,"",dirname(ogf), fixed = T),fixed = T), fixed = T)
  names(ogf)<-dirs

  if(verbose)
    cat("Reading and parsing", length(ogf),"files\n")

  ogl1 <- rbindlist(mclapply(names(ogf), mc.cores = ncores, mc.preschedule = T, function(i){
    print(i)
    og <- readLines(ogf[[i]])
    og <- lapply(og, function(x) strsplit(x, " ")[[1]])
    og.name <- sapply(og, function(x) x[1])
    og.length <- sapply(og, length)-1
    og.ids <- sapply(og, function(x) x[-1])
    og2 = data.table(block.id = i,
                     og = rep(og.name, og.length),
                     id = unlist(og.ids),
                     stringsAsFactors = F)
    return(og2)
  }))

  if(verbose)
    cat("Compiling metadata\n")

  gffi = gff[,c("id","genome")]
  setkey(gffi, id)
  setkey(ogl1,id)
  ogo = merge(gffi, ogl1)
  setkey(ogo, block.id, genome, og, id)
  ogo[, og.n.genes:=length(unique(id)), by = list(block.id, og)]
  ogo[, og.n.genomes:=length(unique(genome)), by = list(block.id, og)]

  ogo.meta = ogo[!duplicated(ogo[,-c(1:2),with = F]),-c(1:2), with = F]
  ogo.sum = ogo.meta[,list(n.og.singleGene = sum(og.n.genes == 1),
                           n.og.singleGenome = sum(og.n.genomes == 1 & og.n.genes > 1),
                           n.og.RBHortho = sum(og.n.genes == 2 & og.n.genomes == 2),
                           n.og.multiOrtho = sum(og.n.genes > 2 & og.n.genomes == 2)),
                     by = list(block.id)]

  if(verbose)
    cat("Constructing graphs\n")
  ogn <- rbindlist(lapply(names(ogf), function(i){
    print(i)
    og <- readLines(ogf[[i]])
    og <- lapply(og, function(x) strsplit(x, " ")[[1]])
    og.name <- sapply(og, function(x) x[1])
    og.length <- sapply(og, length)-1
    og.ids <- sapply(og, function(x) x[-1])
    ogn = sapply(og.ids, length)
    ognet1 = unlist(og.ids[ogn == 1])
    ognet1 = data.table(Var1 = ognet1, Var2 = ognet1, stringsAsFactors = F)
    ognet2 = rbindlist(lapply(og.ids[ogn > 1], function(x) expand.grid(x,x)))
    ognet = rbind(ognet1, ognet2)
    return(ognet)
  }))
  if(verbose)
    cat("\tDone!\n")
  return(list(metadata = ogo.sum, orthogroup.datatable = ogo))
}

make_ofInputInBlk <- function(initial.blast.dir,
                              out.dir,
                              ogff,
                              species.mappings,
                              of.speciesIDs,
                              verbose = T){
  if(verbose)
    cat("Copying orthofinder output to", out.dir,"\n")
  tmp.blast.dir = file.path(out.dir,"tmp")
  tmp.blast.dir.files = file.path(out.dir,"tmp","blasts")
  if(file.exists(tmp.blast.dir)){
    system(paste("rm -rf", tmp.blast.dir))
  }
  system(paste("mkdir", tmp.blast.dir))
  system(paste("cp -r", initial.blast.dir, out.dir))
  ogn1 = sapply(names(ogff), function(x) strsplit(x, " ")[[1]][1])
  ogn2 = sapply(names(ogff), function(x) strsplit(x, " ")[[1]][2])
  if(verbose)
    cat("Culling blast hits ...\n\t")
  for(i in 1:nrow(species.mappings)){
    x = species.mappings[i,]
    g = ogff[ogn1 == x$ref & ogn2 == x$alt]
    xfile = basename(x$filename)
    if(verbose)
      cat(xfile)
    f = fread(file.path(out.dir,xfile), header = F, stringsAsFactors = F, check.names = F)
    if(verbose)
      cat(paste0(" (", nrow(f)," hits) ... "))
    setkey(f, V1, V2)
    fl = rbindlist(lapply(1:length(g), function(j){
      y = g[[j]]
      genes = y$gene.num
      return(f[f$V1  %in% genes & f$V2 %in% genes,])
    }))
    fl = fl[!duplicated(fl[,1:2,with = F]),]
    if(verbose)
      cat("retained", nrow(fl),"hits\n\t")
    write.table(fl,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = file.path(out.dir,xfile))
  }

  for(i in genomeIDs){
    ogl = ogff[ogn1 == i | ogn2 == i]
    ogs = lapply(ogl, function(x) x$gene.num[x$genome == i])
    ognum = of.speciesIDs$genome.num[of.speciesIDs$genome == i]
    xfile = paste0("Blast",ognum,"_",ognum,".txt")
    if(verbose)
      cat(xfile)
    f = fread(file.path(out.dir,xfile), header = F, stringsAsFactors = F, check.names = F)
    if(verbose)
      cat(paste0(" (", nrow(f)," hits) ... "))
    setkey(f, V1, V2)
    fl = rbindlist(lapply(1:length(ogs), function(j){
      genes = ogs[[j]]
      return(f[f$V1  %in% genes & f$V2 %in% genes,])
    }))
    fl = fl[!duplicated(fl[,1:2,with = F]),]
    if(verbose)
      cat("retained", nrow(fl),"hits\n\t")
    write.table(fl,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = file.path(out.dir,xfile))
  }

  if(verbose)
    cat("Done!\n\t")
}

get_ofIDs = function(sgff, of.geneIDs, of.speciesIDs, verbose = T){
  if(verbose)
    cat("Adding orthofinder IDs to gffs ... ")
  of.geneIDs <- data.table(of.geneIDs)
  setkey(of.geneIDs, "id")
  of.speciesIDs <- data.table(of.speciesIDs)
  setkey(of.speciesIDs, "genome")

  ofg.out = lapply(sgff, function(x){
    setkey(x,"genome")
    x1 = merge(of.speciesIDs, x)
    setkey(x1, "id")
    x2 = merge(of.geneIDs, x1)
  })
  names(ofg.out) <- names(sgff)
  if(verbose)
    cat("Done")
  return(ofg.out)
}
