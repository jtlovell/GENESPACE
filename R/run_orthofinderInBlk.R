#' @title run_orthofinderInBlk
#' @description
#' \code{run_orthofinderInBlk} run_orthofinderInBlk
#' @name run_orthofinderInBlk
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param genome1 character string specifying the first genome to analyze
#' @param genome2 character string specifying the second genome to analyze
#' @param onlyInBuffer logical, should only inbuffer hits be retained?
#' @param overwrite logical, should results be overwritten?
#' \cr
#' If called, \code{run_orthofinderInBlk} returns its own arguments.
#'
#' @details info here

#' @title Run orthofinder in block
#' @description
#' \code{pull_inblkOgs} Loop through the pairwise combination of genomes and
#' run orthofinderInBlk
#' @rdname run_orthofinderInBlk
#' @import data.table
#' @export
run_orthofinderInBlk <- function(gsParam,
                                 overwrite = FALSE,
                                 onlyInBuffer = TRUE){

  lab <- query <- target <- blkID <- sameOg <- sameInblkOg <- inBuffer <-
    inblkOG <- ofID <- og <- NULL
  ##############################################################################
  # -- 1. Get metadata together
  bed <- read_combBed(gsParam$synteny$combBed)
  md <- data.table(gsParam$synteny$blast)
  if(!"lab" %in% colnames(md))
    md[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]

  ##############################################################################
  # -- 2 for each line in metadata
  hitsInOgs <- rbindlist(lapply(1:nrow(md), function(i){
    cat(sprintf("\t...%s ", md$lab[i]))
    inblkOgs <- with(md[i,], ofInBlk_engine(
      gsParam = gsParam, genome1 = query, genome2 = target,
      overwrite = overwrite))
    out <- subset(inblkOgs, !is.na(blkID) & (sameOg | sameInblkOg))
    if(onlyInBuffer)
      out <- subset(out, inBuffer)
    out <- out[,c("ofID1", "ofID2")]
    return(out)
  }))

  ##############################################################################
  # -- 3. add vector to the bed file
  ic <- with(hitsInOgs, clus_igraph(id1 = ofID1, id2 = ofID2))
  bed[,inblkOG := ic[ofID]]

  # -- convert NAs to unique orthogroups
  whna <- which(is.na(bed$inblkOG))
  mx <- max(bed$inblkOG, na.rm = T)
  bed$inblkOG[whna] <- paste((mx + 1): (mx + length(whna)))

  # -- replace og with inblkOG and write
  bed[,og := inblkOG]
  write_combBed(x = bed, filepath = gsParam$synteny$combBed)

  return(gsParam)
}

#' @title Run orthofinder within blocks
#' @description
#' \code{run_orthofinderInBlk} Main engine to run orthofinder in blocks
#' @rdname run_orthofinderInBlk
#' @import data.table
#' @import parallel
#' @export
ofInBlk_engine <- function(gsParam,
                           genome1,
                           genome2,
                           overwrite){

  query <- target <- ofID1 <- ofID2 <- blkID <- rid <- uid1 <- uid2 <- regID <-
    ofID1 <- ofID2 <- sameHog <- hog1 <- hog2 <- sameInblkOg <- hasSelf <- NULL


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
  nCores <- gsParam$params$nCores
  md <- data.table(gsParam$synteny$blast)
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
    cat("|| no non-self blocks\n")
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
        p0md$synHits, showProgress = F, na.strings = c("NA", ""),
        select = blNames[1:12])
      bl00[,`:=`(ofID1 = sprintf("%s_g1", ofID1),
                 ofID2 = sprintf("%s_g1", ofID2))]
      bl11 <- fread(
        p1md$synHits, showProgress = F, na.strings = c("NA", ""),
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

    ofDirs <- rbindlist(mclapply(names(spl01), mc.cores = nCores, function(j){
      pdir <- file.path(tmpDirs[j], "pep")
      odir <- file.path(tmpDirs[j], "orthofinder")
      if(dir.exists(odir))
        unlink(odir, recursive = T)
      ofComm <- sprintf("-f %s -op -o %s", pdir, odir)
      outp <- system2(
        ofcall, ofComm, stdout = TRUE, stderr = TRUE)
      fp <- list.files(
        path = odir, pattern = "SequenceIDs.txt", full.names = T, recursive = T)
      out <- data.table(regID = j, path = dirname(fp))
      return(out)
    }))
    # names(ofDirs) <- names(spl01)

    ########
    # -- 2.5 read in sequenceIDs, re-name blasts, write blasts, make command
    write_blast <- function(blastHits, filepath){
      fwrite(
        blastHits, file = filepath,
        sep = "\t", quote = F, row.names = F,
        col.names = F, showProgress = F)
    }
    read_hogog <- function(path, genomeIDs, allowInBlkOGs){
      HOG <- list.files(path, pattern = "N0.tsv", recursive = T, full.names = T)
      if(length(HOG) ==  1){
        ogOut <- parse_hogs(filepath = HOG)
        setnames(ogOut, 1, "ogID")
        return(ogOut)
      }else{
        if(allowInBlkOGs){
          OG <- list.files(
            path, pattern = "Orthogroups.tsv",
            recursive = T, full.names = T)
          if(length(OG) == 1){
            ogOut <- parse_ogs(filepath = OG, genomeIDs = genomeIDs)
          }else{
            ogOut <- NULL
          }
        }
      }
      return(ogOut)
    }

    inblkHOGs <- rbindlist(mclapply(1:nrow(ofDirs), mc.cores = nCores, function(k){
      print(k)
      j <- ofDirs$regID[k]
      regj <- ofDirs$path[k]
      sidf <- file.path(regj, "SequenceIDs.txt")
      if(file.exists(sidf)){
        sids <- read_orthofinderSequenceIDs(sidf)
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
        write_blast(y01, file.path(regj, "Blast0_1.txt.gz"))
        write_blast(y10, file.path(regj, "Blast1_0.txt.gz"))
        write_blast(y00, file.path(regj, "Blast0_0.txt.gz"))
        write_blast(y11, file.path(regj, "Blast1_1.txt.gz"))

        comm <- sprintf("-b %s -t 1 -a 1 -X", regj)

        outp <- system2(ofcall, comm, stdout = TRUE, stderr = TRUE)
        ogout <- read_hogog(
          path = regj, genomeIDs = c("g0", "g1"), allowInBlkOGs = TRUE)
        if(!is.null(ogout))
          ogout[,regID := j]
      }else{
        ogout <- NULL
      }
      return(ogout)
    }))

    ########
    # -- 3.3 merge hogs with blast files
    splhog <- split(inblkHOGs, by = "regID")
    hogout <- rbindlist(lapply(names(splhog), function(j){
      y <- spl01[[j]]
      z <- splhog[[j]]
      if(length(z) == 4){
        hogv <- z$ogID; names(hogv) <- z$id
        y[,`:=`(hog1 = hogv[ofID1], hog2 = hogv[ofID2])]
        y[,sameHog := hog1 == hog2]
        return(subset(y, sameHog)[,c("ofID1", "ofID2")])
      }else{
        return(NULL)
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
  write_synBlast(allBlast, filepath = x$synHits)
  return(allBlast)
}
