#' @title Run the exonerate program
#'
#' @description
#' \code{pipe_exonerateOrphans} A simple wrapper to run orthofinder from R.
#'
#' @param orphan.blast Results from pipe_orphanBlast
#' @param dir.list List of directories, produced by check_environment
#' @param buffer numeric, the amount of space around the focal region to feed
#' to exonerate.
#' @param verbose Logical, should updates be printed?
#' @param ... Not currently in use
#' @details  ...

#' @return Nothing, writes results to the blast.dir directory
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
pipe_exonerateOrphans <- function(orphan.blast,
                                 dir.list,
                                 buffer = 1e3,
                                 verbose = T){

  ########################################################
  pipe_getfasta <- function(bed.dt, tmp.dir, ass.dir,
                           verbose = T){
    if(dir.exists(tmp.dir))
      unlink(tmp.dir, recursive = T)
    dir.create(tmp.dir)

    n <- nrow(bed.dt)
    nb <- ifelse(n < 1000, 100,
                 ifelse(n < 5000, 500,
                        ifelse(n < 10000, 1000, 5000)))

    fa.locs <- sapply(1:nrow(bed.dt), function(i){
      if (verbose)
        if (i %% nb == 0)
          cat(i,"/",nrow(bed.dt),"\n\t")

      x <- bed.dt[i,]
      bed <- x[,c("chr","start","end","reg.name")]
      ass.file <- file.path(assembly.dir,paste0(x$genome,".fa"))
      fa.file <- file.path(tmp.dir,paste0(x$reg.name,".fa"))
      bed.file <- file.path(tmp.dir,paste0(x$reg.name,".bed"))
      run_getfasta(bed = bed,
                   ass.file = ass.file,
                   fa.file = fa.file,
                   bed.file = bed.file)
      return(fa.file)
    })
    bed.dt$fa.file <- fa.locs
    return(bed.dt)
  }
  ########################################################
  write_peptideByGene <- function(pep.fastas,
                                 tmp.dir,
                                 bed.dt,
                                 verbose  = T){
    if (verbose) {
      n <- nrow(bed.dt)
      nb <- ifelse(n < 1000, 100,
                   ifelse(n < 5000, 500,
                          ifelse(n < 10000, 1000, 5000)))
    }

    ids <- unique(bed.dt$id)

    fa.loc <- sapply(1:length(ids), function(i){
      if (verbose)
        if (i %% nb == 0)
          cat(i,"/",length(ids),"\n\t")
      x <- ids[i]
      out.loc <- file.path(tmp.dir, paste0(x,".fa"))
      writeXStringSet(pep.fastas[x], filepath = out.loc)
      return(out.loc)
    })

    out <- data.table(id = ids, pep.file = fa.loc)
    setkey(out, id)
    setkey(bed.dt, id)
    ret <- merge(bed.dt, out)
    return(ret)
  }
  ########################################################
  load.annotations <- function(genomeIDs,
                               cds.dir,
                               peptide.dir,
                               assembly.dir){

    cds.fastas <- do.call(c, lapply(genomeIDs, function(i)
      readDNAStringSet(file.path(cds.dir,
                                 paste0(i, ".fa")))))

    pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
      readAAStringSet(file.path(peptide.dir,
                                paste0(i, ".fa")))))

    fais <- rbindlist(lapply(genomeIDs, function(i){
      tmp <- read.delim(file.path(assembly.dir,
                                  paste0(i, ".fa.fai")),
                        header = F,
                        stringsAsFactors = F,
                        col.names = c("chr",  "chr.length",
                                      "v1", "v2", "v3"))[, 1:2]
      tmp$genome <- i
      return(tmp)
    }))
    setkey(fais, genome, chr)
    return(list(peptide = pep.fastas,
                cds = cds.fastas,
                chrlen = fais))
  }
  ########################################################
  parse_blastLoc <- function(bl, fais, bed.list, genome.list){
    ## Process bed and genome lists to block metadata
    bd <- rbindlist(bed.list)
    bd$genome <- unlist(genome.list)
    genome.info <- lapply(bl$V1, function(x) strsplit(x,".", fixed = T)[[1]])
    setnames(bd,"id","block.id")
    setkey(bd, genome,chr)
    setkey(fais, genome,chr)
    bs <- merge(fais, bd)
    setnames(bs,4:5,c("genome.start","genome.end"))

    ## Convert blast to standard format and merge with metadata
    setnames(bl, c("block.id", "id", "perc.iden",
                   "align.length", "n.mismatch",
                   "n.gapOpen", "q.start", "q.end",
                   "s.start", "s.end",
                   "eval", "score"))
    setkey(bl,block.id)
    setkey(bs, block.id)
    blo <- merge(bs, bl)

    ## Generate start and end genomic locations for the hits
    blo$tend <- with(blo, genome.start+q.start)
    blo$tstart <- with(blo, genome.start+q.end)
    blo$length <- with(blo, tend - tstart)
    blo$start <- with(blo, ifelse(length < 0, tend, tstart))
    blo$end <- with(blo, ifelse(length > 0, tend, tstart))
    blo$tend <- NULL
    blo$tstart <- NULL
    blo$length <- with(blo, end - start)

    ## Reorder and return
    blo$neg.score <- blo$score * (-1)
    blo$neg.len <- blo$length * (-1)
    setkey(blo, id, neg.score, neg.len)
    return(blo)
  }
  ########################################################
  make_blastRegion <- function(bl,
                              max.dist2besthit = 5e3,
                              min.blk.dist2drop = 1e6){

    # -- pull the best hits
    bl[,rank := frank(neg.score, ties.method = "random"),
       by = list(id, block.id)]
    bs <- bl[bl$rank == 1,c("block.id","id","start","end")]
    setnames(bs, 3:4, c("best.start","best.end"))

    # -- Merge best with all hits to calc dist to best
    setkey(bl, id, block.id)
    setkey(bs, id, block.id)
    bo <- merge(bl, bs)
    bo$dist2best <- with(bo, abs(best.start - start))
    bo <- bo[bo$dist2best <= max.dist2besthit,]

    # --  Generate new coordinates for hits, spanning small gaps
    bi <- bo[,list(start = min(start),
                  end = max(end),
                  mean.score = mean(score),
                  mean.length = mean(abs(length))),
            by = list(genome, chr, block.id, id)]
    bi$length <- with(bi, end-start)

    # toss identical hits
    bi <- data.table(bi[!duplicated(bi[,c("genome","chr","id","start","end")]),])
    return(bi)
  }
  ########################################################
  add_buffer <- function(bl, fais, buffer){
    bo <- merge(fais, bl)
    bo$start <- bo$start - buffer
    bo$end <- bo$end + buffer
    bo$end[bo$end>bo$chr.length] <- bo$chr.length[bo$end>bo$chr.length]
    bo$start[bo$start<=0] <- 1
    bo$reg.name <- with(bo, paste0(id,"_",block.id))
    return(bo)
  }
  ########################################################
  pipe_exonerate <- function(locs,
                             verbose = T,
                             min.score = 20){
    if (verbose) {
      n <- nrow(locs)
      nb <- ifelse(n < 1000, 100,
                   ifelse(n < 5000, 500,
                          ifelse(n < 10000, 1000, 5000)))
    }

    exon.out <- lapply(1:nrow(locs), function(i){
      if (verbose)
        if (i %% nb == 0)
          cat(i,"/",nrow(locs),"\n\t")
      x <- locs[i]
      bedt <- data.frame(x[,c("chr","start","end","reg.name")])
      colnames(bedt)[4] <- "id"
      exon <- run_exonerate(ass.fasta.file = x$fa.file,
                           pep.fasta.file = x$pep.file,
                           bed = bedt,
                           min.score = min.score,
                           tmp.dir = tmp.dir)
      return(exon)
    })
    exon.cds <- do.call(c, lapply(exon.out, function(x) x$cds.seq))
    gff.gene <- rbindlist(lapply(exon.out, function(x) x$gff))
    gff.gene$ex.start <- with(gff.gene, as.numeric(start)+as.numeric(exon.start))
    gff.gene$ex.end <- with(gff.gene, as.numeric(start)+as.numeric(exon.end))
    gff.out <- gff.gene[,list(start = min(ex.start),
                              end = max(ex.end)),
                        by = list(id, chr, strand)]
    return(list(gff.region = gff.out,
                gff.cds = gff.gene,
                cds = exon.cds))
  }
  ########################################################
  calc_weightedId <- function(gff.cds, gff.region){
    gff.cds$length <- with(gff.cds, as.numeric(exon.end) - as.numeric(exon.start))
    gff.cds$identity <- sapply(gff.cds$align.info, function(x)
      as.numeric(gsub(" ","",gsub("identity","",strsplit(x,";", fixed = T)[[1]][3]))))

    gff.cds[,tot.len := sum(length),
            by = id]
    gff.cds$identity <- as.numeric(gff.cds$identity)
    gff.cds$tot.len <- as.numeric(gff.cds$tot.len)
    gff.cds$length <- as.numeric(gff.cds$length)
    gff.cds$wt.id <- with(gff.cds, identity * (length/tot.len))

    out <- gff.cds[,list(weighted.id = sum(wt.id),
                        cds.length = max(tot.len)),
                  by = id]
    setkey(out, id)
    setkey(gff.region, id)
    ret <- merge(gff.region, out)
    return(ret)
  }
  ########################################################
  ########################################################

  ########################################################
  has.blast <- sapply(orphan.blast$blast.results, function(x) length(x) > 1)
  blast.in <- rbindlist(orphan.blast$blast.results[has.blast])


  cds.dir = dir.list$cds
  peptide.dir = dir.list$peptide
  assembly.dir = dir.list$assembly
  block.dir = dir.list$block
  tmp.dir = dir.list$tmp

  bed.list = orphan.blast$bed.list
  genome.list = orphan.blast$genome.list
  ########################################################
  if (verbose)
    cat("Preparing blast data for exonerate ...\n")
  if (verbose)
    cat("\tImporting annotation data ... ")
  annot <- load.annotations(genomeIDs = genomeIDs,
                            cds.dir = cds.dir,
                            peptide.dir = peptide.dir,
                            assembly.dir = assembly.dir)
  fais <- annot$chrlen
  pep.fastas <- annot$peptide
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("\tParsing blast results ... ")
  blast <- parse_blastLoc(bed.list = bed.list,
                          genome.list = genome.list,
                          bl = blast.in,
                          fais = fais)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("\tFinding best blast hit regions ... ")
  hit.reg <- make_blastRegion(bl = blast)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("\tAdding",buffer,"bp buffer to ends of blast regions ... ")
  hits2exonerate <- add_buffer(bl = hit.reg,
                               fais = fais,
                               buffer = buffer)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Writing region sequences to file ... Completed:\n\t")
  locs <- pipe_getfasta(bed.dt= hits2exonerate,
                       tmp.dir = tmp.dir,
                       ass.dir = assembly.dir,
                       verbose = T)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Writing peptide sequences to file ... Completed:\n\t")
  locs <- write_peptideByGene(bed.dt= locs,
                             tmp.dir = tmp.dir,
                             pep.fastas = pep.fastas,
                             verbose  = T)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Running exonerate ... Completed:\n\t")
  exonerate.out <- pipe_exonerate(locs)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Compiling exonerate results ... ")
  exon.sum <- with(exonerate.out,
                  calc_weightedId(gff.cds =  gff.cds,
                                  gff.region = gff.region))
  if (verbose)
    cat("Done!\n")

  return(list(file.locations = locs,
         proc.blast.hits = hits2exonerate,
         exonerate.out = exonerate.out,
         exonerate.summary = exon.sum))

}
