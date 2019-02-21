#' @title Simple orthogroup-constrained tandem array inference.
#'
#' @description
#' \code{find_arrayClusters} Simple orthogroup-constrained tandem array inference
#' @param map the blast dataset to screen for syntenic hits
#' @param rerank Logical, should ranks be remade before each step?
#' @param clean.radius Passed on to clean_blocks
#' @param clean.mappings Passed on to clean_blocks
#' @param merge.buffer Passed on to merge_blocks
#' @param min.hits1 The minimum amount of hits in either genome
#' to be counted as an array
#' @param min.hits2 The minimum amount of hits in the over-represented
#' genome to be counted as an array
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
find_tandemArrays <- function(map,
                              min.hits1 = 1,
                              min.hits2 = 3,
                              cull.gff2map = F,
                              method = "cluster",
                              gff = NULL,
                              mcscan.dir = NULL,
                              DCTA.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX/downstream_analyses/detect_collinear_tandem_arrays",
                              verbose = T,
                              MCScanX.param = "-m 50 -w 10 -s 10 -b 0 -k",
                              ...){

  #######################################################
  #######################################################
  cluster_map2arrays <- function(map,
                                   min.hits1 = 2,
                                   min.hits2 = 3,
                                   ...){
    map$array.id <- with(map, paste0(block.id, "|", og1))
    map$unique.genome <- with(map, paste(genome1, genome2))
    spl = split(map, "unique.genome")

    pw.cta <- rbindlist(lapply(spl, function(x){
      if (verbose)
        cat("\t",x$genome1[1], "<-->", x$genome2[1], "... ")

      x <- data.table(x)
      x[,n.hits1 := length(unique(id1)),
          by = list(array.id)]

      x[,n.hits2 := length(unique(id2)),
          by = list(array.id)]

      x$array.id[with(x,
                        (n.hits1 <= min.hits1 &
                           n.hits2 <= min.hits2) |
                          (n.hits1 <= min.hits2 &
                             n.hits2 <= min.hits1))]<-NA
      if (verbose)
        cat("found",
            length(unique(x$og1[!is.na(x$array.id)])),
            "tandem arrays across",
            sum(!is.na(x$array.id)), "blast hits\n")
      return(x)
    }))

    return(pw.cta)
  }
  #######################################################
  #######################################################
  prep_gff4mcscan <- function(gff,
                              genomeIDs){
    gff <- gff[gff$genome %in% genomeIDs,]
    gff[,chr.num := frank(chr, ties.method = "dense"),
        by = genome]
    gff$genome <- factor(gff$genome,
                         levels = genomeIDs)
    gff$rank.start = frank(gff, genome, chr, start,
                           ties.method = "dense")
    gff$rank.end = frank(gff, genome, chr, end,
                         ties.method = "dense")
    gff[,genome.num := frank(genome,
                             ties.method = "dense")]
    lets <- paste0(letters,letters)[1:length(unique(gff$genome))]
    gff$genome.abbrev <- paste0(lets[gff$genome.num], gff$chr.num)
    gff <- gff[,c("genome.abbrev", "id", "start", "end")]
    return(gff)
  }
  #######################################################
  #######################################################
  prep_blast4mcscan <- function(map,
                                genomeIDs){
    map <- map[with(map,
                    genome2 %in% genomeIDs &
                      genome1 %in% genomeIDs), ]
    map <- map[,c("id1", "id2",
                  "perc.iden", "align.length",
                  "n.mismatch", "n.gapOpen",
                  "q.start", "q.end",
                  "s.start", "s.end",
                  "eval", "score")]
    return(map)
  }
  #######################################################
  #######################################################
  propagate_connectedDCTAs <- function(dcta.dt){
    aids2check <- unique(dcta.dt$array.id[duplicated(dcta.dt$og) |
                                            duplicated(dcta.dt$og, fromLast = T)])
    while (length(aids2check) > 0) {
      x <- aids2check[1]
      tid <- dcta.dt$array.id[dcta.dt$id %in% dcta.dt$id[dcta.dt$array.id %in% x]]
      tod <- vector()
      while (length(tod) < length(tid)) {
        tod <- tid
        tid <- dcta.dt$array.id[dcta.dt$id %in% dcta.dt$id[dcta.dt$array.id %in% tid]]
      }
      dcta.dt$array.id[dcta.dt$array.id %in% unique(tid)] <- tid[1]
      aids2check <- aids2check[!aids2check %in% tid]
    }
    dcta.dt <- dcta.dt[!duplicated(dcta.dt),]
    return(dcta.dt)
  }
  #######################################################
  #######################################################
  propagate_connectedOGs <- function(dcta.dt){
    ud <- dcta.dt[!duplicated(dcta.dt[,c("og","array.id")]),]
    dup <- duplicated(ud$og)

    aids2check <- ud$array.id[which(dup)]

    while (length(aids2check) > 0) {
      x = aids2check[1]
      tid = dcta.dt$array.id[dcta.dt$og %in% dcta.dt$og[dcta.dt$array.id %in% x]]
      tod = vector()
      while (length(tod) < length(tid)) {
        tod = tid
        tid = dcta.dt$array.id[dcta.dt$og %in% dcta.dt$og[dcta.dt$array.id %in% tid]]
      }
      dcta.dt$array.id[dcta.dt$array.id %in% unique(tid)] <- tid[1]
      aids2check <- aids2check[!aids2check %in% tid]
    }
    dcta.dt <- dcta.dt[!duplicated(dcta.dt),]
    return(dcta.dt)
  }
  #######################################################
  #######################################################
  load_DCTAs <- function(dtca.output,
                         map){
    og.dt <- with(map,
                  data.table(id = c(id1, id2),
                             og = c(og1, og2)))
    og.dt <- og.dt[!duplicated(og.dt),]
    raw <- fread(dtca.output, header = T)
    r2 <- apply(raw, 1, function(x)
      strsplit(
        paste(x[1], x[3],
              x[2], x[4],
              sep = ","), ",")[[1]])

    for (i in 1:length(r2))
      r2[[i]] <- data.table(id = r2[[i]],
                            array.id = i)
    r2 <- rbindlist(r2)
    r2 <- r2[!duplicated(r2),]
    r3 <- propagate_connectedDCTAs(r2)
    r3 <- r3[!duplicated(r3),]
    r3 <- merge(og.dt, r3, by = "id")
    r4 <- propagate_connectedOGs(r3)
    r4 <- r4[!duplicated(r4),]
    r4$array.id <- frank(r4$array.id, ties.method = "dense")
    return(r4)
  }
  #######################################################
  #######################################################
  rerank_map <- function(map,
                         gff,
                         cull.gff2map){

    if(cull.gff2map)
      gff <- gff[gff$id %in% unique(c(map$id1, map$id2))]

    rr = with(map,
              rerank_fromIDs(id1 = id1,
                             id2 = id2,
                             gff = gff))
    map$rank1 <- NULL
    map$rank2 <- NULL
    map <- merge(rr[,c("id1", "id2", "rank1", "rank2")],
                 map,
                 by = c("id1","id2"))
  }
  #######################################################
  #######################################################
  prep_gff4mcscan <- function(gff,
                              genomeIDs){
    gff <- gff[gff$genome %in% genomeIDs,]
    gff[,chr.num := frank(chr, ties.method = "dense"),
        by = genome]
    gff$genome <- factor(gff$genome,
                         levels = genomeIDs)
    gff$rank.start = frank(gff, genome, chr, start,
                           ties.method = "dense")
    gff$rank.end = frank(gff, genome, chr, end,
                         ties.method = "dense")
    gff[,genome.num := frank(genome,
                             ties.method = "dense")]
    lets <- paste0(letters,letters)[1:length(unique(gff$genome))]
    gff$genome.abbrev <- paste0(lets[gff$genome.num], gff$chr.num)
    gff <- gff[,c("genome.abbrev", "id", "start", "end")]
    return(gff)
  }
  #######################################################
  #######################################################
  prep_blast4mcscan <- function(map,
                                genomeIDs){
    map <- map[with(map,
                    genome2 %in% genomeIDs &
                      genome1 %in% genomeIDs), ]
    map <- map[,c("id1", "id2",
                  "perc.iden", "align.length",
                  "n.mismatch", "n.gapOpen",
                  "q.start", "q.end",
                  "s.start", "s.end",
                  "eval", "score")]
    return(map)
  }
  #######################################################
  #######################################################
  propagate_connectedDCTAs <- function(dcta.dt){
    aids2check <- unique(dcta.dt$array.id[duplicated(dcta.dt$og) |
                                            duplicated(dcta.dt$og, fromLast = T)])
    while (length(aids2check) > 0) {
      x <- aids2check[1]
      tid <- dcta.dt$array.id[dcta.dt$id %in% dcta.dt$id[dcta.dt$array.id %in% x]]
      tod <- vector()
      while (length(tod) < length(tid)) {
        tod <- tid
        tid <- dcta.dt$array.id[dcta.dt$id %in% dcta.dt$id[dcta.dt$array.id %in% tid]]
      }
      dcta.dt$array.id[dcta.dt$array.id %in% unique(tid)] <- tid[1]
      aids2check <- aids2check[!aids2check %in% tid]
    }
    dcta.dt <- dcta.dt[!duplicated(dcta.dt),]
    return(dcta.dt)
  }
  #######################################################
  #######################################################
  propagate_connectedOGs <- function(dcta.dt){
    ud <- dcta.dt[!duplicated(dcta.dt[,c("og","array.id")]),]
    dup <- duplicated(ud$og)

    aids2check <- ud$array.id[which(dup)]

    while (length(aids2check) > 0) {
      x = aids2check[1]
      tid = dcta.dt$array.id[dcta.dt$og %in% dcta.dt$og[dcta.dt$array.id %in% x]]
      tod = vector()
      while (length(tod) < length(tid)) {
        tod = tid
        tid = dcta.dt$array.id[dcta.dt$og %in% dcta.dt$og[dcta.dt$array.id %in% tid]]
      }
      dcta.dt$array.id[dcta.dt$array.id %in% unique(tid)] <- tid[1]
      aids2check <- aids2check[!aids2check %in% tid]
    }
    dcta.dt <- dcta.dt[!duplicated(dcta.dt),]
    return(dcta.dt)
  }
  #######################################################
  #######################################################
  load_DCTAs <- function(dtca.output,
                         map){
    og.dt <- with(map,
                  data.table(id = c(id1, id2),
                             og = c(og1, og2)))
    og.dt <- og.dt[!duplicated(og.dt),]
    raw <- fread(dtca.output, header = T)
    r2 <- apply(raw, 1, function(x)
      strsplit(
        paste(x[1], x[3],
              x[2], x[4],
              sep = ","), ",")[[1]])

    for (i in 1:length(r2))
      r2[[i]] <- data.table(id = r2[[i]],
                            array.id = i)
    r2 <- rbindlist(r2)
    r2 <- r2[!duplicated(r2),]
    r3 <- propagate_connectedDCTAs(r2)
    r3 <- r3[!duplicated(r3),]
    r3 <- merge(og.dt, r3, by = "id")
    r4 <- propagate_connectedOGs(r3)
    r4 <- r4[!duplicated(r4),]
    r4$array.id <- frank(r4$array.id, ties.method = "dense")
    return(r4)
  }
  #######################################################
  #######################################################
  find_collinearArrays<- function(gff,
                                  map,
                                  mcscan.dir,
                                  rerank = T,
                                  DCTA.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX/downstream_analyses/detect_collinear_tandem_arrays",
                                  verbose = T,
                                  MCScanX.param = "-s 2 -m 250 -w 1",
                                  ...){



    if (rerank){
      if (verbose)
        cat("Using the gff gene ranks as positions\n")
      rr = with(map,
                rerank_fromIDs(id1 = id1,
                               id2 = id2,
                               gff = gff))
      map$rank1 <- NULL
      map$rank2 <- NULL
      map <- merge(rr[,c("id1", "id2", "rank1", "rank2")],
                   map,
                   by = c("id1","id2"))

    }
    if (verbose)
      cat("Running tandem array searches for each pair of genomes ... \n")
    genome.cmbns <- map[, c("genome1", "genome2")]
    genome.cmbns <- genome.cmbns[!duplicated(genome.cmbns), ]
    genome.cmbns <- genome.cmbns[genome.cmbns$genome1 != genome.cmbns$genome2, ]

    pw.cta <- rbindlist(apply(genome.cmbns, 1, function(x){
      if (verbose)
        cat("\t",x[1], "<-->", x[2], "... ")
      gff.in <- prep_gff4mcscan(gff = gff,
                                genomeIDs = x)
      blast.in <- prep_blast4mcscan(map = map,
                                    genomeIDs = x)

      gff.file <- file.path(mcscan.dir, "xyz.gff")
      blast.file <- file.path(mcscan.dir, "xyz.blast")
      collin.file <- file.path(mcscan.dir, "xyz.collinearity")
      dcta.file <- file.path(mcscan.dir, "xyz.TCA.txt")

      write.table(gff.in,
                  file = gff.file,
                  row.names = F,
                  col.names = F,
                  quote = F, sep = "\t")
      write.table(blast.in,
                  file = blast.file,
                  row.names = F,
                  col.names = F,
                  quote = F, sep = "\t")

      com <- paste("MCScanX -a",
                   MCScanX.param,
                   file.path(mcscan.dir,"xyz"),
                   "&> /dev/null")
      system(com)

      com <- paste(DCTA.path,
                   "-b", blast.file,
                   "-g", gff.file,
                   "-c", collin.file,
                   "-o", dcta.file, "&> /dev/null")
      system(com)

      dtca <- load_DCTAs(dcta.file,
                         map = map)
      dtca$unique.genome <- paste(x[1], x[2])
      if (verbose)
        cat("found",
            length(unique(dtca$array.id)),
            "tandem arrays across",
            nrow(dtca), "blast hits\n")
      return(dtca)
    }))
    map$array.id <- with(map, paste0(block.id, "|", og1))

    map$unique.genome <- with(map, paste(genome1, genome2))
    u.id <- paste(map$og1, map$unique.genome)
    p.id <- paste(pw.cta$og, pw.cta$unique.genome)
    wh.in.u = which(!u.id %in% p.id)
    map$array.id[wh.in.u]<-NA
    return(map)
  }
  #######################################################

  #######################################################
  #######################################################
  find_collinearArrays <- cmpfun(find_collinearArrays)
  load_DCTAs <- cmpfun(load_DCTAs)
  propagate_connectedOGs <- cmpfun(propagate_connectedOGs)
  propagate_connectedDCTAs <- cmpfun(propagate_connectedDCTAs)
  prep_blast4mcscan <- cmpfun(prep_blast4mcscan)
  prep_gff4mcscan <- cmpfun(prep_gff4mcscan)
  #######################################################
  #######################################################

  #######################################################
  if(!method %in% c("cluster","MCScanX") |
     length(method) != 1){
    stop("method must be one of either cluster or MCScanX\n")
  }
  if(verbose)
    cat("Re-ranking blast hits ... ")
  map <- rerank_map(map = map,
                    gff = gff,
                    cull.gff2map = cull.gff2map)
  map$array.id <- with(map, paste0(block.id, "|", og1))
  if(verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if(!method %in% c("cluster","MCScanX") |
     length(method) != 1){
    stop("method must be one of either cluster or MCScanX\n")
  }
  #######################################################

  #######################################################
  if (method == "cluster") {
    if(verbose)
      cat("Dropping intra-genomic hits ... ")
    map <- map[with(map, genome1 != genome2), ]
    if(verbose)
      cat("Done!\n")
    if (verbose)
      cat("Clustering blast hits within orthogroups and blocks ... \n")
    out <- cluster_map2arrays(map = map,
                                      min.hits1 = 1,
                                      min.hits2 = 3)
    if(verbose)
      cat("\tFound",length(unique(out$array.id )),
          "tandem arrays across",
          sum(!is.na(out$array.id)),"blast hits\n\tDone!\n")

  }else{
    #######################################################

    #######################################################
    if(verbose)
      cat(" ... ")
    out = find_collinearArrays(gff = gff,
                                      map = map,
                                      rerank = T,
                                      mcscan.dir = mcscan.dir,
                                      MCScanX.param = MCScanX.param,
                                      DCTA.path = DCTA.path)
    if (verbose)
      cat("found",
          length(unique(out$og1[!is.na(out$array.id)])),
          "tandem arrays across",
          sum(!is.na(out$array.id)), "blast hits\n")
  }
  #######################################################
  return(out)
}
