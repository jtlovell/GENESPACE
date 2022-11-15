#' @title construct syntenic orthogroups from pairwise blast
#' @description
#' \code{build_synOGs} integrates many pairwise results from synteny into
#' vectors of syntenic orthogroups. Also can re-run orthofinder within blocks
#' and re-calculate orthogroups from those results.
#' @name build_synOGs
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param genome1 character string specifying the first genome to analyze
#' @param genome2 character string specifying the second genome to analyze
#' @param onlyInBuffer logical, should only inbuffer hits be retained?
#' @param overwrite logical, should results be overwritten?
#' \cr
#' If called, \code{build_synOGs} returns its own arguments.
#'
#' @details info here


#' @title build syntenic orthogroups
#' @description
#' \code{build_synOGs} main function to aggregate syntenic orthogroups into a
#' single vector
#' @rdname build_synOGs
#' @import data.table
#' @export
build_synOGs <- function(gsParam){

  # -- read in the combined bed file
  bed <- read_combBed(gsParam$synteny$combBed)

  cat("\t##############\n\tAggregating syntenic orthogroups ... ")
  ofID <- ofID1 <- ofID2 <- synOG <- bedFile <- NULL
  soh <- pull_synOgs(gsParam = gsParam)
  ic <- with(soh, clus_igraph(id1 = ofID1, id2 = ofID2))
  bed[,synOG := ic[ofID]]
  whna <- which(is.na(bed$synOG))
  mx <- max(bed$synOG, na.rm = T)
  bed$synOG[whna] <- paste((mx + 1): (mx + length(whna)))
  cat("Done!\n")

  if(gsParam$params$orthofinderInBlk){
    cat("\t##############\n\tRunning Orthofinder within syntenic regions\n")
    ofID <- ofID1 <- ofID2 <- inblkOG <- og <- NULL
    # -- get the graph of inBlock orthogroups
    synOgHits <- pull_inblkOgs(gsParam)

    ic <- with(soh, clus_igraph(id1 = ofID1, id2 = ofID2))
    bed[,inblkOG := ic[ofID]]
    whna <- which(is.na(bed$inblkOG))
    mx <- max(bed$inblkOG, na.rm = T)
    bed$inblkOG[whna] <- paste((mx + 1): (mx + length(whna)))
    bed[,og := inblkOG]
    write_combBed(x = bed, filepath = bedFile)

    cat("\t##############\n\tRe-annotating the blast files\n")
    gsParam <- annotate_blast(gsParam = gsParam)
    cat("\t##############\n\tRe-running synteny\n")
    gsParam <- synteny(gsParam)

  }else{
    inblkOG <- og <- NULL
    bed[,`:=`(inblkOG = NA, og = synOG)]
    write_combBed(x = bed, filepath = gsParam$synteny$combBed)
  }
  return(gsParam)
}

#' @title extract syntenic orthogroups from hits
#' @description
#' \code{pull_synOgs} engine to pull pairwise hits in syntenic orthogroups
#' @rdname build_synOGs
#' @import data.table
#' @export
pull_synOgs <- function(gsParam, onlyInBuffer = TRUE){
  blkID <- sameOg <- sameInblkOg <- inBuffer <- blkID <- NULL
  ##############################################################################
  # -- 1. Get metadata together
  md <- data.table(gsParam$synteny$blast)
  # -- for each line in metadata
  hitsInOgs <- rbindlist(lapply(1:nrow(md), function(i){
    allBlast <- fread(
      md$synHits[i], showProgress = F, na.strings = c("NA", ""))

    # -- subset to only hits in regions
    out <- subset(allBlast, !is.na(blkID) & sameOg)
    out[,sameInblkOg := TRUE]
    if(onlyInBuffer)
      out <- subset(out, inBuffer)
    out <- subset(out, (sameOg | sameInblkOg))
    out <- out[,c("ofID1", "ofID2")]
    return(out)
  }))
  return(hitsInOgs)
}

#' @title Run orthofinder in block
#' @description
#' \code{pull_inblkOgs} Loop through the pairwise combination of genomes and
#' run orthofinderInBlk
#' @rdname build_synOGs
#' @import data.table
#' @export
pull_inblkOgs <- function(gsParam,
                          overwrite = FALSE,
                          onlyInBuffer = TRUE){

  lab <- query <- target <- blkID <- sameOg <- sameInblkOg <- inBuffer <- NULL
  ##############################################################################
  # -- 1. Get metadata together
  md <- data.table(gsParam$annotBlastMd)
  md[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]
  # -- for each line in metadata
  hitsInOgs <- rbindlist(lapply(1:nrow(md), function(i){
    cat("\t...", md$lab[i])
    inblkOgs <- with(md[i,], run_orthofinderInBlk(
      gsParam = gsParam, genome1 = query, genome2 = target, overwrite = overwrite))
    out <- subset(inblkOgs, !is.na(blkID) & (sameOg | sameInblkOg))
    if(onlyInBuffer)
      out <- subset(out, inBuffer)
    out <- out[,c("ofID1", "ofID2")]
    return(out)
  }))

  return(hitsInOgs)
}

#' @title Run orthofinder within blocks
#' @description
#' \code{run_orthofinderInBlk} Main engine to run orthofinder in blocks
#' @rdname build_synOGs
#' @import data.table
#' @import parallel
#' @export
run_orthofinderInBlk <- function(gsParam,
                                 genome1,
                                 genome2,
                                 overwrite,
                                 allow){

  query <- target <- ofID1 <- ofID2 <- blkID <- rid <- uid1 <- uid2 <-
    ofID1 <- ofID2 <- sameHog <- hog1 <- hog2 <- sameInblkOg <- NULL


  blNames <- c(
    "ofID1", "ofID2", "pid", "length", "mismatches","gapopenings", "queryStart",
    "queryEnd", "subjectStart", "subjectEnd", "Evalue", "bitScore", "rid")
  blNamesR <- c(
    "ofID2", "ofID1", "pid", "length", "mismatches","gapopenings",
    "subjectStart", "subjectEnd", "queryStart", "queryEnd", "Evalue",
    "bitScore", "rid")
  ##############################################################################
  # -- 1. Get the data read in
  ########
  # -- 1.1 parse the annotated blast metadata

  md <- data.table(gsParam$annotBlastMd)
  x <- subset(md, query == genome1 & target == genome2)

  if(nrow(x) == 0){
    tmp <- genome2
    genome2 <- genome1
    genome1 <- tmp
  }

  if(nrow(x) > 1)
    stop(sprintf("genome1: %s and genome2: %s are not unique in gsParam\n",
                 genome1, genome2))
  ########
  # 1.2 read in and parse the annotated blast file
  allBlast <- read_synHits(x$synHits)

  # -- add unique identifier to geneIDs (in case of intragenomic)

  # -- subset to only hits in regions
  sb01 <- subset(allBlast, !is.na(blkID))

  cat(sprintf("%s synhits || %s in global OGs ", nrow(sb01), sum(sb01$sameOg)))
  dontRerun <- "sameInblkOg" %in% colnames(allBlast)
  if(dontRerun)
    dontRerun <- all(!is.na(allBlast$sameInblkOg)) & !overwrite
  if(dontRerun)
    cat(sprintf("|| %s in inBlk OGs (pre-existing file)\n",
                sum(allBlast$sameInblkOg)))

  sb01[,hasSelf := any(ofID1 %in% ofID2), by = "blkID"]
  onlySelf <- sum(!sb01$hasSelf) < gsParam$params$blkSize
  if(onlySelf)
    cat("|| no non-self blocks for inBlk Orthofinder\n")
  sb01 <- subset(sb01, !hasSelf)

  sb01[,`:=`(ofID1 = sprintf("%s_g1", ofID1),
             ofID2 = sprintf("%s_g2", ofID2))]

  if(!onlySelf && !dontRerun){
    # -- make a new variable `rid` that is a unique regionID vector
    sb01[,rid := sprintf("reg%s", as.numeric(as.factor(blkID)))]
    targetGenome <- sb01$genome2[1]
    queryGenome <- sb01$genome1[1]
    # -- subset to only regions with enough genes in each genome
    sb01[,`:=`(uid1 = uniqueN(ofID1), uid2 = uniqueN(ofID2)), by = "rid"]
    sb01 <- subset(sb01, uid1 >= 40 & uid2 >= 40)

    # -- make an ofID --> id dictionary
    di1 <- sb01$ofID1;  names(di1) <- sb01$id1
    di2 <- sb01$ofID2;  names(di2) <- sb01$id2

    ########
    # -- 1.3 read in and parse the peptides
    peps0 <- read_aaFasta(file.path(
      gsParam$paths$peptide, sprintf("%s.fa", queryGenome)))
    peps1 <- read_aaFasta(file.path(
      gsParam$paths$peptide, sprintf("%s.fa", targetGenome)))

    # -- subset peptides to only genes in the hits
    peps0 <- peps0[unique(sb01$id1)]
    peps1 <- peps1[unique(sb01$id2)]

    # -- replace names with unique ofIDs
    names(peps0) <- di1[names(peps0)]
    names(peps1) <- di2[names(peps1)]

    ########
    # -- 1.4 build blast-only column data.tables
    bl01 <- sb01[,blNames, with = F]

    # -- invert for reverse alignment
    bl10 <- sb01[,blNamesR, with = F]
    setnames(bl10, blNames)

    # -- if self hits, just re-name
    if(targetGenome == queryGenome){
      bl00 <- data.table(bl01)
      bl11 <- data.table(bl10)
    }else{
      # -- if not, read in and parse the intragenomic files
      p0md <- subset(md, query == target & query == queryGenome)
      p1md <- subset(md, query == target & query == targetGenome)
      bl00 <- fread(
        p0md$annotBlastFile, showProgress = F, na.strings = c("NA", ""),
        select = blNames[1:12])
      bl00[,`:=`(ofID1 = sprintf("%s_g1", ofID1),
                 ofID2 = sprintf("%s_g1", ofID2))]
      bl11 <- fread(
        p1md$annotBlastFile, showProgress = F, na.strings = c("NA", ""),
        select = blNames[1:12])
      bl11[,`:=`(ofID1 = sprintf("%s_g2", ofID1),
                 ofID2 = sprintf("%s_g2", ofID2))]
    }

    ##############################################################################
    # -- 2. set up the environment

    ########
    # -- 2.1 make the tmp directory
    tmpDir <- gsParam$paths$tmp
    if(dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    dir.create(tmpDir)

    ########
    # -- 2.2 create the directories for each region
    blkIDs <- unique(bl01$rid)
    tmpDirs <- sapply(blkIDs, function(j){
      pth <- file.path(tmpDir, j)
      if(dir.exists(pth))
        unlink(pth, recursive = T)
      dir.create(pth)
      dir.create(file.path(pth, "pep"))
      return(pth)
    })

    ########
    # -- 2.3 subset and write peptides
    spl01 <- split(sb01, by = "rid")
    for(j in blkIDs){
      y <- spl01[[j]]
      p0 <- peps0[unique(y$ofID1)]
      p1 <- peps1[unique(y$ofID2)]
      writeXStringSet(p0, filepath = file.path(tmpDirs[j], "pep", "g0.fa"))
      writeXStringSet(p1, filepath = file.path(tmpDirs[j], "pep", "g1.fa"))
    }

    ########
    # -- 2.4 make the orthofinder input files
    ofcall <- gsParam$shellCalls$orthofinder
    ofDirs <- sapply(names(spl01), function(j){
      pdir <- file.path(tmpDirs[j], "pep")
      odir <- file.path(tmpDirs[j], "orthofinder")
      if(dir.exists(odir))
        unlink(odir, recursive = T)
      ofComm <- sprintf("-f %s -op -o %s", pdir, odir)
      outp <- system2(
        ofcall,  ofComm, stdout = TRUE, stderr = TRUE)
      fp <- list.files(
        path = odir, pattern = "SequenceIDs.txt", full.names = T, recursive = T)
      return(fp)
    })

    ########
    # -- 2.5 read in sequenceIDs, re-name blasts, write blasts
    write_blast <- function(blastHits, filepath){
      fwrite(
        blastHits, file = filepath,
        sep = "\t", quote = F, row.names = F,
        col.names = F, showProgress = F)
    }

    bl <- sapply(names(spl01), function(j){
      if(file.exists(ofDirs[j])){
        sids <- read_orthofinderSequenceIDs(ofDirs[j])
        idv <- sids$ofID; names(idv) <- sids$id
        y <- data.table(spl01[[j]])
        u1 <- unique(y$ofID1)
        u2 <- unique(y$ofID2)
        y01 <- subset(bl01, ofID1 %in% u1 & ofID2 %in% u2)
        y10 <- subset(bl10, ofID1 %in% u2 & ofID2 %in% u1)
        y00 <- subset(bl00, ofID1 %in% u1 & ofID2 %in% u1)
        y11 <- subset(bl11, ofID1 %in% u2 & ofID2 %in% u2)
        y01[,`:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
        y10[,`:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
        y00[,`:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
        y11[,`:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
        write_blast(y01, file.path(dirname(ofDirs[j]), "Blast0_1.txt.gz"))
        write_blast(y10, file.path(dirname(ofDirs[j]), "Blast1_0.txt.gz"))
        write_blast(y00, file.path(dirname(ofDirs[j]), "Blast0_0.txt.gz"))
        write_blast(y11, file.path(dirname(ofDirs[j]), "Blast1_1.txt.gz"))
      }
    })

    ##############################################################################
    # -- 3. Run orthofinder within block
    ########
    # -- 3.1 make the orthofinder commands
    comms <- sapply(names(spl01), function(j){
      if(file.exists(ofDirs[j]))
        return(sprintf("-b %s -t 1 -a 1 -X", dirname(ofDirs[j])))
    })


    ########
    # -- 3.2 call orthofinder and get the paths to the HOG files
    allowInBlkOGs <- TRUE
    hogf <- mclapply(names(spl01), mc.cores = gsParam$params$nCores, function(j){
      if(file.exists(ofDirs[j])){
        outp <- system2(
          ofcall,
          comms[j],
          stdout = TRUE, stderr = TRUE)
        HOG <- list.files(
          dirname(ofDirs[j]),
          pattern = "N0.tsv",
          recursive = T, full.names = T)
        if(length(HOG) ==  1){
          tmp <- parse_hogs(filepath = HOG, genomeIDs = c("g0", "g1"))
          setnames(tmp, 1, "ogID")
          return(tmp)
        }else{
          if(allowInBlkOGs){
            OG <- list.files(
              dirname(ofDirs[j]),
              pattern = "Orthogroups.tsv",
              recursive = T, full.names = T)
            return(parse_ogs(filepath = OG, genomeIDs = c("g0", "g1")))
          }
        }
      }
    })
    names(hogf) <- names(spl01)

    ########
    # -- 3.3 merge hogs with blast files
    hogout <- rbindlist(lapply(names(spl01), function(j){
      y <- spl01[[j]]
      z <- hogf[[j]]
      if(length(z) == 3){
        hogv <- z$ogID; names(hogv) <- z$id
        y[,`:=`(hog1 = hogv[ofID1], hog2 = hogv[ofID2])]
        y[,sameHog := hog1 == hog2]
        return(subset(y, sameHog)[,c("ofID1", "ofID2")])
      }
    }))

    ########
    # -- 3.4 get unique pairs of genes in the same inblk HOG
    u <- with(hogout, paste(gsub("_g1$", "", ofID1), gsub("_g2$", "", ofID2)))
    allBlast[,sameInblkOg := paste(ofID1, ofID2) %in% u]

    cat(sprintf("|| %s in inBlk OGs\n", sum(allBlast$sameInblkOg)))
  }

  ########
  # -- 3.5 write the blast file, with the new column `sameInblkOg`
  write_synBlast(allBlast, filepath = x$annotBlastFile)
  return(allBlast)
}
