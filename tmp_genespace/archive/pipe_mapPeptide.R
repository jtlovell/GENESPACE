pipe_mapping <- function(genomeIDs,
                         reciprocal = TRUE,
                         mapping.program = "diamond",
                         calculate_diamond_db = TRUE,
                         n.cores = 1,
                         verbose = TRUE){


  if(mapping.program == "blat"){
    com.fun = function(db,
                       pep,
                       out){
      paste("blat", pep, db,
            "-prot -out=blast8 -mask=lower -qMask=lower",
            out, "&>/dev/null")
    }
  }else{
    if(mapping.program == "diamond"){
      com.fun = function(pep,
                         db,
                         out,
                         n.cores = 1,
                         diamond.opts = "",
                         switch.names = F){

        if(!switch.names){
          cols = "6 sseqid qseqid pident length mismatch gapopen sstart send qstart qend evalue bitscore"
        }else{
          cols = "6"
        }
        return(
          paste("diamond blastp",
              "--query",pep,
              "--db", db,
              "--min-score",10,
              "--max-target-seqs",12,
              "--threads",n.cores,
              "--quiet",
              diamond.opts,
              "--outfmt",cols,
              "--out",out))
      }

    }else{
      stop("Only blat and diamond mapping.programs are currently implemented")
    }
  }


  if(mapping.program == "diamond" & calculate_diamond_db){
    if(verbose)
      cat("Generating diamond databases ...\n")

    for(i in genomeIDs){
      if(verbose)
        cat("\t ...",i,"\n")
      system(paste("diamond makedb",
                   "--in",file.path(peptide.dir, paste0(i,".fa")),
                   "--db", file.path(blast.dir, i),
                   "--threads",n.cores,
                   "--quiet"))
    }
  }

  map.index = combn(genomeIDs,2, simplify = F)
  if(verbose)
    cat("Running", mapping.program, "searches :\n")

  if(mapping.program == "diamond"){
    results = lapply(1:length(map.index), function(i){
      x = map.index[[i]]
      if(reciprocal){
        if(verbose)
          cat("\t",x[1],"vs.",x[2],paste0("... (",i,"/",length(map.index),")"),"\n")
        system(
          com.fun(pep = file.path(peptide.dir, paste0(x[2],".fa")),
                  db = file.path(blast.dir,paste0(x[1],".dmnd")),
                  out = file.path(blast.dir,paste0(x[1],"_",x[2],".blast")),
                  n.cores = n.cores))
        if(verbose)
          cat("\t",x[2],"vs.",x[1],"\n")
        system(
          com.fun(pep = file.path(peptide.dir, paste0(x[1],".fa")),
                  db = file.path(blast.dir,paste0(x[2],".dmnd")),
                  out = file.path(blast.dir,paste0(x[2],"_",x[1],".blast")),
                  n.cores = n.cores, switch.names = T))
      }else{
        if(verbose)
          cat("\t",x[1],"vs.",x[2],paste0("... (",i,"/",length(map.index),")"),"\n")
        system(
          com.fun(pep = file.path(peptide.dir, paste0(x[2],".fa")),
                  db = file.path(blast.dir,paste0(x[1],".dmnd")),
                  out = file.path(blast.dir,paste0(x[1],"_",x[2],".blast")),
                  n.cores = n.cores))
      }
    })
  }else{
    results = parallel::mclapply(map.index,
                                 mc.cores = n.cores,
                                 mc.preschedule = F, function(x){
      if(reciprocal){
        if(verbose)
          cat("\t",x[1],"vs.",x[2],"\n")
        system(
          com.fun(pep = file.path(peptide.dir, paste0(x[2],".fa")),
                  db = file.path(peptide.dir,paste0(x[1],".fa")),
                  out = file.path(blast.dir,paste0(x[1],"_",x[2],".blast"))))
        if(verbose)
          cat("\t",x[2],"vs.",x[1],"\n")
        system(
          com.fun(pep = file.path(peptide.dir, paste0(x[1],".fa")),
                  db = file.path(peptide.dir,paste0(x[2],".fa")),
                  out = file.path(blast.dir,paste0(x[2],"_",x[1],".blast"))))

      }else{
        if(verbose)
          cat("\t",x[1],"vs.",x[2],"\n")
        system(
          com.fun(pep = file.path(peptide.dir, paste0(x[2],".fa")),
                  db = file.path(peptide.dir,paste0(x[1],".dmnd")),
                  out = file.path(blast.dir,paste0(x[1],"_",x[2],".blast"))))
      }
    })
  }

  if(verbose)
    cat("All BLAST runs, completed!")

  info = data.frame(t(combn(genomeIDs,2, simplify = T)),stringsAsFactors = F)
  colnames(info)<-c("id1","id2")
  info$blast.file = file.path(blast.dir,paste0(info$id1,"_",info$id2,".blast"))
  if(reciprocal){
    info$recip.file = file.path(blast.dir,paste0(info$id2,"_",info$id1,".blast"))
  }else{
    info$recip.file = NA
  }
  return(info)
}

test = parse_blastResults(genomeIDs = genomeIDs,
                          gff.dir = gff.dir,
                          info = blast.info,
                          ploidy = ploidy)
parse_blastResults = function(genomeIDs,
                              gff.dir,
                              info,
                              verbose = TRUE,
                              ploidy,
                              min.propMax = .5,
                              min.score = 30,
                              nmapsPerHaplotype = 1,
                              eps.radius = c(100,50,20),
                              n.mappingWithinRadius = c(10,10,10),
                              plotit = T){
  require(data.table)

  if(verbose)
    cat("Parsing gff files\n")
  gff.files <- list.files(gff.dir,
                          full.names = T)
  names(gff.files) <- gsub(".gff3$", "",
                           basename(gff.files))

  parse_gff <- function(gff){
    g <- suppressWarnings(
      data.table::fread(gff,
                        showProgress = F,
                        verbose = F))
    g <- g[g$V3 == "gene",c(9,1,4,5,7)]
    g$V9 <- sapply(g$V9, function(x) gsub("Name=", "",
                                          strsplit(x,";")[[1]][2]))
    data.table::setnames(g, c("id","chr","start","end","strand"))
    return(g)
  }

  gff <- rbindlist(lapply(names(gff.files), function(i){
    tmp <- parse_gff(gff.files[[i]])
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr","start")],
                       ties.method = "random")
    return(tmp)
  }))
  gff1 <- data.table(gff)
  gff2 <- data.table(gff)
  setnames(gff1, paste0(colnames(gff),"1"))
  setkey(gff1,id1,genome1)
  setnames(gff2, paste0(colnames(gff),"2"))
  setkey(gff2,id2,genome2)

  blast.in = lapply(1:nrow(info), function(i){

    if(verbose)
      cat("Parsing",info$id1[i],"mapped against", info$id2[i],"\n")
    f1 = fread(info$blast.file[i], showProgress = F,
               col.names = c("id1", "id2", "pident", "length", "mismatch", "gapopen",
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    if(!is.na(info$recip.file[i])){
      f2 = fread(info$recip.file[i], showProgress = F,
                 col.names = c("id1", "id2", "pident", "length", "mismatch", "gapopen",
                               "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
      raw = data.table(rbind(f1, f2))
    }else{
      raw = data.table(f1)
    }
    raw$genome1 = info$id1[i]
    raw$genome2 = info$id2[i]


    raw$negbitscore = raw$bitscore*(-1)
    setkey(raw,id1,id2,negbitscore)
    uniq <- raw[!duplicated(raw[,c("id1","id2"), with = F]),]
    if(verbose)
      cat("\tRetained", nrow(uniq),"unique mappings\n")

    d <- uniq[uniq$bitscore >= min.score,]

    d[, rank1 := frank(bitscore, ties.method = "dense"),
      by = list(id1)]
    d[, rank2 := frank(bitscore, ties.method = "dense"),
      by = list(id2)]

    ploidy1 <- ploidy[info$id1[i]]
    ploidy2 <- ploidy[info$id2[i]]

    cull <- data.table(d[(d$rank1 <= ploidy1*nmapsPerHaplotype |
                           d$rank2 <= ploidy2*nmapsPerHaplotype) &
                           d$rank1 <= ploidy2*nmapsPerHaplotype*2 &
                              d$rank2 <= ploidy1*nmapsPerHaplotype*2,])

    cull[, prop1 := bitscore/max(bitscore),
      by = list(id1)]
    cull[, prop2:= bitscore/max(bitscore),
         by = list(id2)]
    cull2 <- data.table(cull[cull$prop1>=min.propMax |
                               cull$prop2>=min.propMax,])

    if(verbose)
      cat("\tRetained", nrow(cull2),"mappings that pass score thresholds\n")

    setkey(cull2,id2,genome2)
    merged = merge(gff2, cull2)
    setkey(merged,id1,genome1)
    merged = merge(gff1, merged)


    merged$rank1 <- frank(merged[,c("chr1","start1")],
                   ties.method = "dense")
    merged$rank2<- frank(merged[,c("chr2","start2")],
                   ties.method = "dense")

    x = data.table(merged)

    nn <- frNN(x[,c("rank1","rank2")],
               eps = eps.radius[1])
    dbs <- dbscan(nn,
                  minPts = n.mappingWithinRadius[1])

    cull.dbs1 = data.table(x[dbs$cluster != 0,])
    if(verbose)
      cat("\tRetained", nrow(cull.dbs1),"mappings that pass initial density thresholds (" ,n.mappingWithinRadius[1],"hits within",eps.radius[1],") \n")

    x = data.table(cull.dbs1)
    x$rank1 <- frank(x[,c("chr1","start1")],
                     ties.method = "dense")
    x$rank2<- frank(x[,c("chr2","start2")],
                    ties.method = "dense")

    nn <- frNN(x[,c("rank1","rank2")],
               eps = eps.radius[2])
    dbs <- dbscan(nn,
                  minPts = n.mappingWithinRadius[2])

    cull.dbs2 = data.table(x[dbs$cluster != 0,])
    if(verbose)
      cat("\tRetained", nrow(cull.dbs2),"mappings that pass secondary density thresholds (" ,n.mappingWithinRadius[2],"hits within",eps.radius[2],") \n")

    x = data.table(cull.dbs2)
    x$rank1 <- frank(x[,c("chr1","start1")],
                     ties.method = "dense")
    x$rank2<- frank(x[,c("chr2","start2")],
                    ties.method = "dense")

    nn <- frNN(x[,c("rank1","rank2")],
               eps = eps.radius[3])
    dbs <- dbscan(nn,
                  minPts = n.mappingWithinRadius[3])

    cull.dbs3 = data.table(x[dbs$cluster != 0,])
    if(verbose)
      cat("\tRetained", nrow(cull.dbs3),"mappings that pass final density thresholds (" ,n.mappingWithinRadius[3],"hits within",eps.radius[3],") \n")

    if(plotit){
      par(mfrow = c(2,2))
      with(merged, plot(rank1, rank2, pch = ".", col = rgb(0,0,0,.1),
                        xlab = info$id1[i],
                        ylab = info$id2[i],
                        main = "Passed score culling"))
      with(cull.dbs1, plot(rank1, rank2, pch = ".", col = rgb(0,0,0,.1),
                        xlab = info$id1[i],
                        ylab = info$id2[i],
                        main = "Passed initial dbs culling"))
      with(cull.dbs2, plot(rank1, rank2, pch = ".", col = rgb(0,0,0,.1),
                           xlab = info$id1[i],
                           ylab = info$id2[i],
                           main = "Passed second dbs culling"))
      with(cull.dbs1, plot(rank1, rank2, pch = ".", col = rgb(0,0,0,.1),
                           xlab = info$id1[i],
                           ylab = info$id2[i],
                           main = "Passed final dbs culling"))
    }
  })
  return(rbindlist(blast.in))
}

