#' @title Query genespace datasets
#'
#' @description
#' \code{query_genespace} Parse and output underlying data.
#'
#' @name query_genespace
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param genomeIDs character vector, specifying which genomes to use. Defaults
#' to all genomeIDs specification in gsParam.
#' @param pg the pangenome data.table
#' @param refChrom character string of length 1, the reference chromosome to
#' grab from the pangenome
#' @param startOrder numeric or integer specifying the start position for the
#' pangenome query, in linear order.
#' @param endOrder numeric or integer specifying the end position of the query
#' @param writeTo file/path or character string giving file name to write to
#' @param pavSynOnly logical, should PAV be calculated only for syntenic
#' entries?
#' @param cnv2keep integer vector of length matching genomeIDs, specifying the
#' copy number for pangenome entries to be printed. For example, with three
#' genomes, c(1,1,1) would return a single copy orthogroup. The behavior of
#' this argument depends on whether pavSynOnly. If TRUE, then the copy number
#' is only calculated on syntenic orthogroup CNVs. Otherwise, non-syntenic
#' orthologs flagged with * are also used in the copy number calculation.
#' @param pav2keep logical vector of lenght matching genomeIDs. Behaves the
#' same as cnv2keep, but with a logical vector specifying whether any gene is
#' present in that group. For example, c(TRUE, TRUE, TRUE) would return all
#' pangenome entries complete across the three genomes.
#' @param useGlobalPgOrder logical, should the start/endOrder be calculated
#' within the specified refChrom or should the global pangenome order be used?
#'
#' @details ...
#'
#' @return ...
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @note \code{query_genespace} is a generic name for the functions documented.
#' \cr
#' If called, \code{query_genespace} returns its own arguments.
#'
#' @title query_pangenome
#' @description
#' \code{query_pangenome} query_pangenome
#' @rdname query_genespace
#' @import data.table
#' @export
query_pangenome <- function(pg,
                            refChrom = NULL,
                            genomeIDs = NULL,
                            startOrder = NULL,
                            endOrder = NULL,
                            writeTo = NULL,
                            useGlobalPgOrder = FALSE,
                            cnv2keep = NULL,
                            pav2keep = NULL,
                            pavSynOnly = FALSE){

  # -- check the formatting of the pangenome
  if(nrow(pg) < 1)
    stop("can't parse pg; is this a data.table produced by pangenome?\n")
  if(!all(c("chr", "ord", "pgID", "og") %in% colnames(pg)))
    stop("can't find the needed column names in pg; was this produced by pangenome?")

  if(is.null(refChrom)){
    startOrder = NULL
    endOrder = NULL
  }
  if(!is.null(refChrom))
    if(!refChrom %in% pg$chr)
      stop("can't find", refChrom, "in the pangenome\n")

  # -- get genomeIDs in order, using the pangenome
  pg <- as.data.table(pg)
  if(is.null(genomeIDs))
    genomeIDs <- colnames(pg)[-c(1:4)]
  if(!all(genomeIDs %in% colnames(pg)))
    genomeIDs <- colnames(pg)[-c(1:4)]
  genomeIDs <- genomeIDs[!duplicated(genomeIDs)]

  # -- check and prepare the cnv2keep and/or pav2keep
  if(!is.null(cnv2keep)){
    cnv2keep <- as.integer(cnv2keep)
    if(length(cnv2keep) != length(genomeIDs) | any(is.na(cnv2keep)))
      cnv2keep <- NULL
    names(cnv2keep) <- genomeIDs
  }

  if(!is.null(pav2keep)){
    pav2keep <- as.logical(pav2keep)
    if(length(pav2keep) != length(pav2keep) | any(is.na(pav2keep)))
      pav2keep <- NULL
    names(pav2keep) <- genomeIDs
  }

  if(!is.null(pav2keep) & !is.null(cnv2keep))
    pav2keep <- NULL

  # -- get the starting/ending positions if not specified
  chr <- ord <- tmp <- NULL
  if(is.null(startOrder) & !is.null(refChrom))
    startOrder <- 0
  if(is.null(endOrder) & !is.null(refChrom))
    endOrder <- max(pg$ord, na.rm = T)

  # -- subset to the region of interest
  p <- data.table(pg)
  if(!is.null(refChrom)){
    p <- subset(p, chr == refChrom)
    if(!useGlobalPgOrder)
      p[,ord := 1:.N]
    p <- subset(p, ord >= startOrder & ord <= endOrder)
  }

  # -- if only syntenic ogs, parse these
  pc <- data.table(p)
  if(pavSynOnly){
    pc[,(genomeIDs) := lapply(.SD, function(x) lapply(x, function(y)
      y[!grepl("*", y, fixed = T)])),
      .SDcols = genomeIDs, by = "pgID"]
  }

  # -- calculate CNV
  pc[,(genomeIDs) := lapply(.SD, function(x) sapply(x, length)),
     .SDcols = genomeIDs, by = "pgID"]

  # -- calculate PAV
  pav <- data.table(pc)
  pav[,(genomeIDs) := lapply(.SD, function(x) x > 0),
      .SDcols = genomeIDs, by = "pgID"]

  # -- compress lists to character vectors
  pr <- data.table(p)
  pr[,(genomeIDs) := lapply(.SD, function(x) sapply(x, paste, collapse = ";")),
     .SDcols = genomeIDs, by = "pgID"]

  # -- subset to user defined PAV/CNV
  if(!is.null(pav2keep)){
    pav[,tmp := apply(.SD, 1, paste, collapse = ""), .SDcols = genomeIDs]
    wh <- which(pav$tmp == paste(pav2keep, collapse = ""))
    pav[,tmp := NULL]
  }else{
    if(!is.null(cnv2keep)){
      pc[,tmp := apply(.SD, 1, paste, collapse = ""), .SDcols = genomeIDs]
      wh <- which(pc$tmp == paste(cnv2keep, collapse = ""))
      pc[,tmp := NULL]
    }else{
      wh <- 1:nrow(pr)
    }
  }

  # -- write output and return
  if(is.null(writeTo)){
    if(is.null(refChrom)){
      writeTo <- "Pangenome.txt.gz"
    }else{
      writeTo <- sprintf(
        "Pangenome_Chr%s_Start%s_End%s.txt.gz",
        refChrom, min(pr$ord, na.rm = T), max(pr$ord, na.rm = T))
    }
  }
  fwrite(pr[wh,], file = writeTo, sep = "\t", quote = F)
  return(list(raw = p[wh,], pav = pav[wh,], cnts = pc[wh,]))
}


#' @title write_phytozome
#' @description
#' \code{write_phytozome} write_phytozome
#' @rdname query_genespace
#' @import data.table
#' @export
write_phytozome <- function(gsParam,
                            genomeIDs = NULL){

  g1 <- g2 <- refGenome <- isSelf <- ofID1 <- ofID2 <- blkAnchor <- orient <- NULL
  # -- get the output directory
  outdir <- file.path(gsParam$results, "phytozome")
  dir.create(outdir)

  # -- read in the block coordinates
  blksFile <- file.path(gsParam$results, "syntenicBlocks.txt.gz")
  blk <- fread(blksFile)

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs

  # -- load the annotated gff
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, showProgress = F, na.strings = c("", "NA"))

  # -- find the hit files
  pfs <- CJ(g1 = genomeIDs, g2 = genomeIDs)
  pfs[,fs := file.path(gsParam$paths$results,
                       sprintf("%s_%s_synHits.txt.gz", g1, g2))]
  pfs <- subset(pfs, file.exists(fs))
  fs <- pfs$fs

  # -- load the hits (either regions or blks)
  so <- rbindlist(mclapply(fs, mc.cores = gsParam$params$nCores, function(i){
    x <- fread(i,
               select = c("ofID1", "ofID2","blkID","blkAnchor","gen1"),
               showProgress = F,
               na.strings = c("NA", "-", ""))
    if(x$gen1 != refGenome){
      setnames(x, c("ofID2", "ofID1","blkID","blkAnchor","gen1"))
      x <- x[,c("ofID1", "ofID2","blkID","blkAnchor","gen1")]
    }
    x[,isSelf := any(ofID1 %in% ofID2), by = "blkID"]
    return(subset(x, blkAnchor & !isSelf)[,1:3])
  }))


  sogFile <- list.files(gsParam$results, pattern = "_synog.txt.gz", full.names = T)
  so <- rbindlist(lapply(sogFile, fread))
  so <- subset(so, ofID1 %in% gff$ofID & ofID2 %in% gff$ofID)
  ov <- gff$ord; names(ov) <- gff$ofID
  cv <- gff$chr; names(cv) <- gff$ofID
  sv <- gff$start; names(sv) <- gff$ofID
  ev <- gff$end; names(ev) <- gff$ofID
  gv <- gff$genome; names(gv) <- gff$ofID
  iv <- gff$id; names(iv) <- gff$ofID
  so[,`:=`(genome1 = gv[ofID1], genome2 = gv[ofID2],
           chr1 = cv[ofID1], chr2 = cv[ofID2],
           start1 = sv[ofID1], start2 = sv[ofID2],
           end1 = ev[ofID1], end2 = ev[ofID2])]
  so <- rbind(so, with(so, data.table(
    ofID1 = ofID2, ofID2 = ofID1, blkID = blkID,
    genome1 = genome2, genome2 = genome1, chr1 = chr2, chr2 = chr1,
    start1 = start2, start2 = start1, end1 = end2, end2 = end1)))
  so <- subset(so, !duplicated(so))
  bo <- rbind(
    with(subset(blk, orient == "+"), data.table(
      proteomeA = c(genome1, genome2),
      chromA = c(chr1, chr2),
      startA = c(bpStart1, bpStart2),
      endA = c(bpEnd1, bpEnd2),
      proteomeB = c(genome2, genome1),
      chromB = c(chr2, chr1),
      startB = c(bpStart2, bpStart1),
      endB = c(bpEnd2, bpEnd1),
      strand = c(orient, orient),
      blkID = c(blkID, blkID))),
    with(subset(blk, orient == "-"), data.table(
      proteomeA = c(genome1, genome2),
      chromA = c(chr1, chr2),
      startA = c(bpStart1, bpStart2),
      endA = c(bpEnd1, bpEnd2),
      proteomeB = c(genome2, genome1),
      chromB = c(chr2, chr1),
      startB = c(bpEnd2, bpEnd1),
      endB = c(bpStart2, bpStart1),
      strand = c(orient, orient),
      blkID = c(blkID, blkID))))
  bo <- subset(bo, !duplicated(bo))
  splb <- split(bo, by = "proteomeA")
  splo <- split(so, by = "genome1")
  nu <- lapply(names(splb), function(i){
    x <- splb[[i]]
    y <- splo[[i]]
    splx <- split(x, by = "proteomeB")
    sply <- split(y, by = "genome2")
    dirx <- file.path(outdir, i)
    dir.create(dirx)
    dir.create(file.path(dirx, "alignments"))
    for(j in names(splx)){
      bcoords <- splx[[j]]
      oc <- sply[[j]]
      splk <- split(oc, by = "blkID")
      fwrite(bcoords[,1:9, with = F],
             file = file.path(dirx, sprintf("%s_%s.csv", i, j)))
      for(k in 1:nrow(bcoords)){
        w <- bcoords[k,]
        print(w)
        z <- with(splk[[w$blkID]], data.table(
          geneA = iv[ofID1], startA = start1, endA = end1,
          geneB = iv[ofID2], startB = start2, endB = end2))
        fn <- with(w, sprintf("%s__%s:%s-%s_%s__%s:%s-%s.csv",
                              proteomeA, chromA, startA, endA,
                              proteomeB, chromB, startB, endB))
        fwrite(z, file = file.path(dirx, "alignments", fn))
      }
    }
  })
}

