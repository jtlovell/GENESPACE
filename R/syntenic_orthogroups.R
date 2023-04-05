#' @title construct syntenic orthogroups from pairwise blast
#' @name syntenic_orthogroups
#' @description
#' \code{syntenic_orthogroups} integrates many pairwise results from synteny into
#' vectors of syntenic orthogroups. Also can re-run orthofinder within blocks
#' and re-calculate orthogroups from those results.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param updateArrays logical, should arrays be updated?
#' @param overwrite logical, should results be overwritten?
#' @param genome1 character vector with the query genome for orthofinderInBlk
#' @param genome2 character vector with the target genome for orthofinderInBlk

#' @details info here

#' @title extract syntenic orthogroups
#' @description
#' \code{syntenic_orthogroups} combine blast hits to cluster hits into
#' syntenic orthogroups.
#' @rdname syntenic_orthogroups
#' @import data.table
#' @export
syntenic_orthogroups <- function(gsParam, updateArrays){

  pull_synOgs <- function(gsParam){
    blkID <- sameOG <- sameInblkOG <- inBuffer <- blkID <- NULL
    md <- data.table(gsParam$synteny$blast)
    nCores <- gsParam$params$nCores
    hitsInOgs <- rbindlist(mclapply(1:nrow(md), mc.cores = nCores, function(i)
      subset(read_synHits(
        md$synHits[i]), sameOG & !noAnchor)[,c("ofID1", "ofID2")]))
    return(hitsInOgs)
  }

  get_ogsFromHits <- function(ofID1, ofID2, bed, colName){

    # -- convert to
    soh <- data.table(ofID1 = ofID1, ofID2 = ofID2)
    ofID1 <- ofID2 <- NULL

    # -- get a vector of graph clusters
    ic <- with(soh, clus_igraph(id1 = ofID1, id2 = ofID2))

    bed[[colName]] <- ic[bed$ofID]

    whna <- which(is.na(bed[[colName]]))
    mx <- max(as.numeric(bed[[colName]]), na.rm = T)
    bed[[colName]][whna] <- paste((mx + 1):(mx + length(whna)))

    return(bed)
  }

  noAnchor <- isArrayRep <- ofID <- newArrayID <- arrayID <- og <- NULL

  # -- read in the combined bed file
  bedf <- file.path(gsParam$paths$results, "combBed.txt")
  bed <- read_combBed(bedf)
  sogs <- pull_synOgs(gsParam)

  bedr <- subset(bed, isArrayRep)
  bedr <- get_ogsFromHits(
    ofID1 = sogs$ofID1, ofID2 = sogs$ofID2, colName = "og", bed = bedr)

  ar <- bedr[,c("arrayID", "og")]
  ar <- subset(ar, !duplicated(arrayID))

  # -- update the array reps

  beda <- data.table(bed)
  beda[,og := NULL]
  beda <- subset(beda, !duplicated(beda))
  ar <- subset(ar, !duplicated(ar))
  bed <- merge(beda, ar, by = "arrayID", all = T, allow.cartesian = TRUE)

  if(updateArrays){
    reps <- add_array2bed(
      bed = subset(bed, !noAnchor & isArrayRep),
      synBuff = gsParam$params$synBuff,
      maxIter = 1,
      reorder = T)

    reps <- add_arrayReps2bed(reps)
    arrayMap <- subset(bed, !noAnchor & isArrayRep)[,c("ofID", "arrayID")]
    newMap <- subset(reps)[,c("ofID", "arrayID")]
    setnames(newMap, "arrayID", "newArrayID")
    arrayMap <- subset(arrayMap, !duplicated(arrayMap))
    newMap <- subset(newMap, !duplicated(newMap))
    arrayMap <- merge(arrayMap, newMap, by = "ofID", all = T, allow.cartesian = T)
    arrayMap <- arrayMap[,c("arrayID", "newArrayID")]
    arrayMap <- subset(arrayMap, !duplicated(arrayMap))
    bed <- subset(bed, !duplicated(bed))
    bed <- merge(bed, arrayMap, by = "arrayID", all.x = T, allow.cartesian = T)
    bed[,isArrayRep := ofID %in% reps$ofID[reps$isArrayRep]]
    bed[,`:=`(arrayID = newArrayID, newArrayID = NULL)]
  }

  write_combBed(x = bed, filepath = bedf)
  return(gsParam)
}

#' @title Run orthofinder in block
#' @description
#' \code{run_orthofinderInBlk} Loop through the pairwise combination of genomes
#' and run orthofinderInBlk
#' @rdname syntenic_orthogroups
#' @import data.table
#' @export
run_orthofinderInBlk <- function(gsParam,
                                 overwrite = FALSE){

  lab <- query <- target <- blkID <- sameOG <- sameInblkOG <-
    inblkOG <- ofID <- og <- ofID1 <- ofID2 <- NULL
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
      gsParam = gsParam,
      genome1 = query,
      genome2 = target,
      overwrite = overwrite))
    if(is.data.frame(inblkOgs)){
      if("sameInblkOG" %in% colnames(inblkOgs)){
        out <- subset(inblkOgs, (sameInblkOG))
        out <- out[,c("ofID1", "ofID2")]
        u <- with(out, paste(c(ofID1, ofID2), c(ofID2, ofID1)))
        hits <- read_synHits(md$synHits[i])
        cat(sprintf("n syn / syn OG / inblk OG hits = %s / %s",
                    nrow(hits), sum(hits$sameOG)))
        hits[,sameOG := sameOG | paste(ofID1, ofID2) %in% u]
        write_synHits(hits, md$synHits[i])
        cat(sprintf(" / %s\n",
                    sum(hits$sameOG)))
        return(out)
      }
    }
  }))
  fwrite(
    hitsInOgs,
    file = file.path(gsParam$paths$results, "allInBlkOgHits.txt.gz"),
    sep = "\t")

  return(gsParam)
}

#' @title Run orthofinder within blocks
#' @description
#' \code{ofInBlk_engine} Main engine to run orthofinder in blocks
#' @rdname syntenic_orthogroups
#' @import data.table
#' @import parallel
#' @export
ofInBlk_engine <- function(gsParam,
                           genome1,
                           genome2,
                           overwrite){

  query <- target <- ofID1 <- ofID2 <- blkID <- rid <- uid1 <- uid2 <- regID <-
    ofID1 <- ofID2 <- sameHog <- hog1 <- hog2 <- sameInblkOG <- hasSelf <-
    isSyntenic <- noAnchor <- isArrayRep1 <- isArrayRep2 <- NULL


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
  allBlast <- subset(
    read_allBlast(x$allBlast),
    !is.na(regID) & !noAnchor & isArrayRep1 & isArrayRep2)


  ########
  # 1.3 check it the run has already happened, exit if so
  sb01 <- data.table(allBlast)
  dontRerun <- any(!is.na(allBlast$sameInblkOG))
  if(dontRerun)
    dontRerun <- all(!is.na(allBlast$sameInblkOG)) & !overwrite
  if(dontRerun){
    cat("found existing run, not rerunning\n")
    return(NULL)
  }

  ########
  # 1.4 Check if it is only self hits, exit if so
  sb01[,hasSelf := any(ofID1 %in% ofID2), by = "regID"]
  onlySelf <- sum(!sb01$hasSelf) < gsParam$params$blkSize
  sb01 <- subset(sb01, !hasSelf)
  sb01[,`:=`(ofID1 = sprintf("%s_g1", ofID1),
             ofID2 = sprintf("%s_g2", ofID2))]
  if(onlySelf){
    cat("no non-self blocks\n")
    return(NULL)
  }

  ########
  # 1.5 Check if the regions are too small, exit if so
  sb01[,rid := sprintf("reg%s", as.numeric(as.factor(regID)))]
  targetGenome <- sb01$genome2[1]
  queryGenome <- sb01$genome1[1]
  sb01[,`:=`(uid1 = uniqueN(ofID1), uid2 = uniqueN(ofID2)), by = "rid"]
  sb01md <- sb01[,list(n = (uid1[1] + uid2[1])/2), by = "rid"]
  propPass <- sum(sb01md$n[sb01md$n >= 40])/sum(sb01md$n)
  if(propPass < 0.5){
    cat(sprintf(
      "<50%% (%s%%) of syn. hits in regions with >=40 genes, not running\n",
      round(propPass,3)*100))
    return(NULL)
  }

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
      p0md$allBlast, showProgress = F, na.strings = c("NA", ""),
      select = blNames[1:12])
    bl00[,`:=`(ofID1 = sprintf("%s_g1", ofID1),
               ofID2 = sprintf("%s_g1", ofID2))]
    bl11 <- fread(
      p1md$allBlast, showProgress = F, na.strings = c("NA", ""),
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
  on.exit(expr = unlink(list.files(tmpDir, full.names = T), recursive = T))

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
  allBlast[,sameInblkOG := paste(ofID1, ofID2) %in% u]

  out <- subset(allBlast, sameInblkOG)
  return(out)
}

