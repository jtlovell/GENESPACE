#' @title build_synOGs
#' @description
#' \code{build_synOGs} build_synOGs
#' @name build_synOGs
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param genome1 xx
#' @param genome1 xx
#' @param hits data.table containing the blast hits, also stored in /synHits
#' \cr
#' If called, \code{build_synOGs} returns its own arguments.
#'
#' @details info here

#' @title build_synOGs
#' @description
#' \code{build_synOGs} build_synOGs
#' @rdname build_synOGs
#' @import data.table
#' @export
build_synOGs <- function(gsParam){
  "xx"
}


#' @title run_orthofinderInBlk
#' @description
#' \code{run_orthofinderInBlk} run_orthofinderInBlk
#' @rdname build_synOGs
#' @import data.table
#' @export
run_orthofinderInBlk <- function(gsParam){
  ##############################################################################
  # -- 1. Get metadata together
  md <- data.table(gsParam$annotBlastMd)
  md[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]
  # -- for each line in metadata
  hitsInOgs <- rbindlist(lapply(1:nrow(md), function(i){
    cat(md$lab[i])
    inblkOgs <- with(md[i,], add_inblkHogs2hits(
      gsParam = gsParam, genome1 = query, genome2 = target))
    out <- subset(inblkOgs, !is.na(regID) & (sameOg | sameInblkOg))
    out <- out[,c("ofID1", "ofID2", "sameOg", "sameInblkOg", "inBuffer")]
    return(out)
  }))

  return(hitsInOgs)
}

#' @title add_inblkHogs2hits
#' @description
#' \code{add_inblkHogs2hits} add_inblkHogs2hits
#' @rdname build_synOGs
#' @import data.table
#' @export
add_inblkHogs2hits <- function(gsParam, genome1, genome2){

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
  print(x)
  if(nrow(x) > 1)
    stop(sprintf("genome1: %s and genome2: %s are not unique in gsParam\n",
                 genome1, genome2))

  ########
  # 1.2 read in and parse the annotated blast file
  allBlast <- fread(
    x$annotBlastFile, showProgress = F, na.strings = c("NA", ""))

  # -- add unique identifier to geneIDs (in case of intragenomic)
  allBlast[,`:=`(ofID1 = sprintf("%s_g1", ofID1),
                 ofID2 = sprintf("%s_g2", ofID2))]

  # -- subset to only hits in regions
  sb01 <- subset(allBlast, !is.na(regID))
  cat(sprintf("%s synhits || %s in global OGs ", nrow(sb01), sum(sb01$sameOg)))

  # -- make a new variable `rid` that is a unique regionID vector
  sb01[,rid := sprintf("reg%s", as.numeric(as.factor(regID)))]

  # -- subset to only regions with enough genes in each genome
  sb01[,`:=`(uid1 = uniqueN(ofID1), uid2 = uniqueN(ofID2)), by = "rid"]
  sb01 <- subset(sb01, uid1 >= 40 & uid2 >= 40)

  # -- make an ofID --> id dictionary
  di <- with(sb01, c(ofID1, ofID2))
  names(di) <- with(sb01, c(id1, id2))

  ########
  # -- 1.3 read in and parse the peptides
  peps0 <- read_aaFasta(file.path(
    gsParam$paths$peptide, sprintf("%s.fa", x$query)))
  peps1 <- read_aaFasta(file.path(
    gsParam$paths$peptide, sprintf("%s.fa", x$target)))

  # -- subset peptides to only genes in the hits
  peps0 <- peps0[unique(sb01$id1)]
  peps1 <- peps1[unique(sb01$id2)]

  # -- replace names with unique ofIDs
  names(peps0) <- di[names(peps0)]
  names(peps1) <- di[names(peps1)]

  ########
  # -- 1.4 build blast-only column data.tables
  bl01 <- sb01[,blNames, with = F]

  # -- invert for reverse alignment
  bl10 <- sb01[,blNamesR, with = F]
  setnames(bl10, blNames)

  # -- if self hits, just re-name
  if(x$query == x$target){
    bl00 <- data.table(bl01)
    bl11 <- data.table(bl10)
  }else{
    # -- if not, read in and parse the intragenomic files
    p0md <- subset(md, query == target & query == x$query)
    p1md <- subset(md, query == target & query == x$target)
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
  if(!dir.exists(tmpDir))
    dir.create(tmpDir)

  ########
  # -- 2.2 create the directories for each region
  regIDs <- unique(bl01$rid)
  tmpDirs <- sapply(regIDs, function(j){
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
  for(j in regIDs){
    y <- spl01[[j]]
    p0 <- peps0[unique(y$ofID1)]
    p1 <- peps1[unique(y$ofID2)]
    writeXStringSet(p0, file = file.path(tmpDirs[j], "pep", "g0.fa"))
    writeXStringSet(p1, file = file.path(tmpDirs[j], "pep", "g1.fa"))
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
  })

  ##############################################################################
  # -- 3. Run orthofinder within block
  ########
  # -- 3.1 make the orthofinder commands
  comms <- sapply(names(spl01), function(j)
    sprintf("-b %s -t 1 -a 1 -X", dirname(ofDirs[j])))

  ########
  # -- 3.2 call orthofinder and get the paths to the HOG files
  hogf <- mclapply(names(spl01), mc.cores = 10, function(j){
    outp <- system2(
      ofcall,
      comms[j],
      stdout = TRUE, stderr = TRUE)
    HOG <- find_ofFiles(
      orthofinderDir = file.path(dirname(ofDirs[j]), "OrthoFinder"))$hogs
    return(HOG)
  })
  hogf <- unlist(hogf); names(hogf) <- names(spl01)

  ########
  # -- 3.3 merge hogs with blast files
  hogout <- rbindlist(lapply(names(spl01), function(j){
    y <- spl01[[j]]
    z <- parse_hogs(path = hogf[j], genomeIDs = c("g0", "g1"))
    hogv <- z$hogID; names(hogv) <- z$id
    y[,`:=`(hog1 = hogv[ofID1], hog2 = hogv[ofID2])]
    y[,sameHog := hog1 == hog2]
    return(subset(y, sameHog)[,c("ofID1", "ofID2")])
  }))

  ########
  # -- 3.4 get unique pairs of genes in the same inblk HOG
  u <- with(hogout, paste(ofID1, ofID2))
  allBlast[,sameInblkOg := paste(ofID1, ofID2) %in% u]
  allBlast[,`:=`(ofID1 = gsub("_g1$", "", ofID1),
                 ofID2 = gsub("_g2$", "", ofID2))]

  cat(sprintf("|| %s in inBlk OGs\n", sum(allBlast$sameInblkOg)))
  ########
  # -- 3.5 write the blast file, with the new column `sameInblkOg`
  fwrite(
    allBlast, file = x$annotBlastFile,
    quote = F, sep = "\t", showProgress = FALSE)
  return(allBlast)
}
