#' @title Build orthofinder database for GENESPACE
#'
#' @name run_orthofinder
#'
#' @description
#' \code{run_orthofinder} GENESPACE routines for running orthofinder within R,
#' generating code to run outside of R, or re-running orthofinder on pre-
#' computed blast results within syntenic regions.
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param overwrite logical, should results be overwritten?
#' @param onlyCheckRun logical, should nothing be done but see if there is a run
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param genomeIDs character vector with the genomes to include in the run
#' @param minGenes4of integer specifying the minimum number of genes needed to
#' run orthofinder.
#' @details ...
#' @return ...
#'
#' @examples
#' \dontrun{
#'
#' runwd <- file.path(getwd(), "testGenespace")
#' make_exampleDataDir(writeDir = runwd)
#'
#' gpar <- init_genespace(
#'   genomeIDs = c("human","chimp","rhesus"),
#'   speciesIDs = c("human","chimp","rhesus"),
#'   versionIDs = c("human","chimp","rhesus"),
#'   ploidy = rep(1,3),
#'   diamondMode = "fast",
#'   orthofinderMethod = "fast",
#'   wd = runwd,
#'   nCores = 4,
#'   minPepLen = 50,
#'   gffString = "gff",
#'   pepString = "pep",
#'   path2orthofinder = "orthofinder",
#'   path2mcscanx = "~/MCScanX",
#'   rawGenomeDir = file.path(runwd, "rawGenomes"))
#'
#' parse_annotations(
#'   gsParam = gpar,
#'   gffEntryType = "gene",
#'   gffIdColumn = "locus",
#'   gffStripText = "locus=",
#'   headerEntryIndex = 1,
#'   headerSep = " ",
#'   headerStripText = "locus=")
#'
#' # -- run_orthofinder is separated from the rest of the pipeline to allow
#' # the user to run it externally if desired. Depending on your system config,
#' # you may be able to run it in fast mode through R (if orthofinder is in the
#' # path), or you may need to run on a separate environment. If the latter case,
#' # run_orthofinder will output a warning and the necessary commands.
#'
#' gpar <- run_orthofinder(gsParam = gpar, overwrite = F)
#'
#' }
#'
#' @rdname run_orthofinder
#' @import data.table
#' @export
run_orthofinder <- function(gsParam,
                            overwrite = FALSE){

  ##############################################################################
  ##############################################################################
  ##############################################################################
  default_ofDb <- function(gsParam){
    if(all(is.na(gsParam$params$synteny)))
      stop("must run set_syntenyParams first\n")

    if(gsParam$params$verbose)
      cat("\tCleaning out orthofinder directory and prepping run\n")
    ############################################################################
    # 1. clean out peptide directory of unused fastas, if necessary
    drop_unusedPeptides(gsParam)

    ############################################################################
    # 2. Remove existing orthofinder directory if it exists
    if(dir.exists(gsParam$paths$orthofinder)){
      unlink(gsParam$paths$orthofinder, recursive = T)
    }

    ############################################################################
    # 3. get command
    if(is.na(gsParam$paths$orthofinderCall)){
      dontRun <- TRUE
      p2of <- "orthofinder"
    }else{
      dontRun <- FALSE
      p2of <- gsParam$paths$orthofinderCall
    }
    if(gsParam$params$verbose & !dontRun)
      cat("\tRunning full orthofinder on pre-computed blast",
          "\n\t##################################################",
          "\n\t##################################################\n")
    com <- sprintf(
      "%s -f %s -t %s -a 1 -X -o %s",
      p2of,
      dirname(gsParam$paths$peptide[1]),
      gsParam$params$nCores,
      gsParam$paths$orthofinder)

    ############################################################################
    # 4. run it
    if(dontRun){
      cat(
        "\tCould not find valid orthofinder executable in the path\n",
        "\tRun the following command outside of R (assuming orthofinder is in the path):",
        "\n################\n",
        sprintf("cd %s", gsParam$paths$orthofinder),
        com,
        "\n################\n", sep = "")
    }else{
      system(com)
    }

    return(com)
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  fast_ofDb <- function(gsParam){
    runBlast <- genome1 <- genome2 <- db2 <- fa1 <- blFile <- NULL
    ############################################################################
    # Ad hoc internal functions
    ############################################################################
    # a. invert blast file
    invert_blast <- function(fileIn, fileOut){
      tmp <- fread(fileIn, verbose = F, showProgress = F)
      tmp <- tmp[,c(2,1,3:6,8,7,10,9,11,12)]
      fwrite(tmp, sep = "\t",
             quote = FALSE,
             col.names = FALSE,
             row.names = FALSE,
             file = fileOut,
             showProgress = FALSE,
             verbose = FALSE)
    }
    ############################################################################
    # b. move and reorganize orthofinder input files
    reorg_ofInput <- function(ofDir){
      origF <- list.files(ofDir, full.names = TRUE)

      ofTmp <- dirname(list.files(
        path = ofDir,
        pattern = "diamondDBSpecies0.dmnd",
        recursive = T,
        full.names = T))
      ofFiles <- list.files(path = ofTmp, full.names = TRUE)
      for(i in ofFiles)
        file.copy(from = i, to = ofDir, overwrite = T)
      unlink(origF, recursive = T, force = T)
    }
    ############################################################################
    # c. add blast metadata / calls
    add_blastInfo2syn <- function(gsParam, ofSpeciesIDs){
      p <- data.table(gsParam$params$synteny)
      diamondMode <- gsParam$params$diamondMode
      ofd <- gsParam$paths$orthofinder
      p[,`:=`(db1 = file.path(ofd, sprintf("diamondDBSpecies%s.dmnd",si[genome1])),
              db2 = file.path(ofd, sprintf("diamondDBSpecies%s.dmnd",si[genome2])),
              fa1 = file.path(ofd, sprintf("Species%s.fa",si[genome1])),
              fa2 = file.path(ofd, sprintf("Species%s.fa",si[genome2])),
              blFile = file.path(ofd, sprintf("Blast%s_%s.txt.gz",si[genome1], si[genome2])),
              invertFile = file.path(ofd, sprintf("Blast%s_%s.txt.gz",si[genome2], si[genome1])))]
      dm <- ifelse(gsParam$params$diamondMode == "--fast", "", gsParam$params$diamondMode)
      p[,com := sprintf(
        "%s blastp %s --quiet -e %s -p %s --compress 1 -d %s -q %s -o %s",
        "diamond", dm, .1, gsParam$params$nCores, db2, fa1, blFile)]
      return(p)
    }

    ############################################################################
    # 1. Remove existing orthofinder directory if it exists
    if(dir.exists(gsParam$paths$orthofinder)){
      unlink(gsParam$paths$orthofinder, recursive = T)
    }

    ############################################################################
    # 2. convert to orthofinder
    com <- sprintf(
      "%s -f %s -t %s -a 1 -op -o %s 1>/dev/null 2>&1",
      gsParam$paths$orthofinderCall,
      dirname(gsParam$paths$peptide[1]),
      gsParam$params$nCores,
      gsParam$paths$orthofinder)
    system(com)

    ############################################################################
    # 3. place orthofinder input files in paths$orthofinder
    reorg_ofInput(gsParam$paths$orthofinder)
    si <- read_orthofinderSpeciesIDs(gsParam$paths$orthofinder)

    ############################################################################
    # 4. get blast parameters
    p <- add_blastInfo2syn(gsParam = gsParam, ofSpeciesIDs = si)

    ############################################################################
    # 5. Run blasts
    pwp <- subset(p, runBlast)
    for(i in 1:nrow(pwp)){
      if(gsParam$params$verbose)
        with(pwp[i,], cat(sprintf("\t\tRunning %s/%s (%s vs. %s)\n",
                                  i, nrow(pwp), genome1, genome2)))
      system(pwp$com[i])
    }

    ############################################################################
    # 6. Invert blasts if necessary
    if(gsParam$params$verbose)
      cat("\t\tDone!\n\tInverting intergenomic files ... ")
    ip <- subset(p, !runBlast & genome1 != genome2)
    for(i in 1:nrow(ip)){
      invert_blast(fileIn = ip$invertFile[i], fileOut = ip$blFile[i])
    }

    ############################################################################
    # 7. Run orthofinder
    if(gsParam$params$verbose)
      cat("Done!\n\tRunning orthofinder -og on pre-computed blast:\n")
    com <- with(gsParam, sprintf(
      "%s -b %s -t %s -a 1 -X -og",
      paths$orthofinderCall, paths$orthofinder, params$nCores, params$nCores))

    system(com)
    return(com)
  }

  ##############################################################################
  ##############################################################################
  ##############################################################################

  if(is.logical(gsParam$paths$orthofinderCall))
    if(!gsParam$paths$orthofinderCall && !is.na(gsParam$paths$orthofinderCall))
      gsParam$paths$orthofinderCall <- NA

  # set the synteny parameters
  if(is.data.table(gsParam$params$synteny))
    if(!all(genomeIDs %in% gsParam$params$synteny$genome1))
      gsParam$params$synteny <- NULL
  if(!is.data.table(gsParam$params$synteny)){
    cat("Synteny Parameters have not been set! Setting to defaults\n")
    gsParam <- set_syntenyParams(gsParam)
  }

  beenRun <- find_orthofinderResults(gsParam, onlyCheckRun = T)
  if(beenRun & !overwrite){
    warning("orthofinder run exists & !overwrite, so not running")
    gsParam <- find_orthofinderResults(gsParam, onlyCheckRun = F)
  }else{
    if(is.na(gsParam$paths$orthofinderCall)){
      com <- default_ofDb(gsParam)
    }else{
      if(gsParam$params$orthofinderMethod == "fast"){
        if(gsParam$params$verbose & gsParam$params$diamondMode == "fast")
          cat("\tRunning 'draft' a.k.a 'fast' genespace orthofinder method with 'fast' diamond mode",
              "\n\t############################################################",
              "\n\t***NOTE***\n\tThis method should only be used for:",
              "\n\t\t(1) closely related diploid species or",
              "\n\t\t(2) visualization/genome QC purposes",
              "\n\tIf you are building a multi-species or polyploid pangenome from global orthogroups ... \n\t\tcancel this and rerun with:\n\t\torthofinderMethod = 'default' or \n\t\tdiamondMode != 'fast'",
              "\n\t############################################################\n")
        if(gsParam$params$verbose & gsParam$params$diamondMode != "fast")
          cat("\tRunning 'draft' a.k.a 'fast' genespace orthofinder method",
              "\n\t############################################################",
              "\n\t***NOTE***\n\tThis method should only be used for:",
              "\n\t\t(1) closely related diploid species,",
              "\n\t\t(2) visualization/genome QC purposes, or",
              "\n\t\t(3) inferring orthogroups WITHIN syntenic regions",
              "\n\t############################################################\n")
        com <- fast_ofDb(gsParam)
      }else{
        if(gsParam$params$verbose)
          cat("\tRunning 'defualt' genespace orthofinder method",
              "\n\t############################################################\n")
        com <- default_ofDb(gsParam)
      }
    }
  }
  return(gsParam)
}

#' @title blkwise_orthofinder
#' @description
#' \code{blkwise_orthofinder} blkwise_orthofinder
#' @rdname run_orthofinder
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom Biostrings readAAStringSet
#' @export
blkwise_orthofinder <- function(gsParam,
                                gff,
                                genomeIDs = NULL,
                                minGenes4of = 40){

  ##############################################################################
  ##############################################################################
  read_hits4of <- function(gsParam, genome1, genome2){

    read_invertBlast <- function(gsParam, genome1, genome2, ofSpId, invert = T){
      V1 <- V2 <- V12 <- NULL
      h <- read_blast(
        path = gsParam$paths$blastDir, onlyIDScore = F,
        ofID1 = ofSpId[genome1],
        ofID2 = ofSpId[genome2])
      if(invert){
        h1 <- h[,c(2,1,3:6,8,7,10,9,11,12)]
        setnames(h1, colnames(h))
        h <- rbind(h, h1)
        setorder(h, -V12)
        h <- subset(h, !duplicated(paste(V1, V2)))
      }
      return(h)
    }

    ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
    if(genome1 == genome2){
      h <- read_invertBlast(
        gsParam = gsParam, genome1 = genome1, genome2 = genome2, ofSpId = ofSpId)
    }else{
      h <- rbind(
        read_invertBlast(
          gsParam = gsParam, genome1 = genome1, genome2 = genome1,
          ofSpId = ofSpId, invert = F),
        read_invertBlast(
          gsParam = gsParam, genome1 = genome2, genome2 = genome2,
          ofSpId = ofSpId, invert = F),
        read_invertBlast(
          gsParam = gsParam, genome1 = genome1, genome2 = genome2,
          ofSpId = ofSpId, invert = F),
        read_invertBlast(
          gsParam = gsParam, genome1 = genome2, genome2 = genome1,
          ofSpId = ofSpId, invert = F))
    }
    return(h)
  }
  ##############################################################################
  ##############################################################################
  run_ofInReg <- function(gsParam, genome1, genome2, blks, gff, pepspl){

    V1 <- V2 <- ord <- ofID1 <- ofID2 <- synOg1 <- synOg2 <- og <- regID <- NULL
    cls <- og <- gen1 <- gen2 <- nHits1 <- genome <- lgRegionID <- NULL
    nCores <- gsParam$params$nCores
    geno1 <- genome1[1]
    geno2 <- genome2[1]
    genome1 <- genome2 <- NULL

    # -- subset the blks object to the genomes of interst
    b <- subset(blks, gen1 == geno1 & gen2 == geno2)
    setorder(b, -nHits1)

    # -- read in the hits
    fs <- file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz",geno1, geno2))
    hits <- fread(
      fs,
      na.strings = c("NA", ""),
      showProgress = F)

    # -- read the raw blasts
    h <- read_hits4of(
      gsParam = gsParam,
      genome1 = geno1,
      genome2 = geno2)

    # -- parse the gff
    gs1 <- split(subset(gff, genome == geno1), by = "chr")
    gs2 <- split(subset(gff, genome == geno2), by = "chr")

    # -- for each block ...
    synOgHits <- rbindlist(mclapply(1:nrow(b), mc.cores = nCores, function(j){
      bj <- b[j,]
      hitj <- subset(hits, lgRegionID == bj$blkID)
      gj1 <- subset(
        gs1[[as.character(bj$chr1)]],
        ord >= bj$startOrd1 & ord <= bj$endOrd1)
      gj2 <- subset(
        gs2[[as.character(bj$chr2)]],
        ord >= bj$minOrd2 & ord <= bj$maxOrd2)
      uj <- unique(c(gj1$ofID, gj2$ofID))
      hj <- subset(h, V1 %in% uj & V2 %in% uj)
      hij <- subset(hits, ofID1 %in% uj & ofID2 %in% uj)
      tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4inBlkTMPdir", j))
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      ogdt <- run_ofFromObj(
        blast00 = subset(hj, V1 %in% gj1$ofID & V2 %in% gj1$ofID),
        blast01 = subset(hj, V1 %in% gj1$ofID & V2 %in% gj2$ofID),
        blast10 = subset(hj, V1 %in% gj2$ofID & V2 %in% gj1$ofID),
        blast11 = subset(hj, V1 %in% gj2$ofID & V2 %in% gj2$ofID),
        pep0 = pepspl[[geno1]],
        pep1 = pepspl[[geno2]],
        writeDir = tmpDir)
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      ogv <- ogdt$og; names(ogv) <- ogdt$ofID
      hij[,`:=`(synOg1 = ogv[ofID1], synOg2 = ogv[ofID2])]
      return(subset(hij, synOg1 == synOg2 | !(is.na(og) | is.na(regID))))
    }))
    synOgHits <- subset(synOgHits, !duplicated(paste(ofID1, ofID2)))
    synOgHits[,cls := ifelse(is.na(og), "onlyInBlk", ifelse(synOg1 != synOg2, "onlyGlob", "both"))]
    return(synOgHits)
  }
  ##############################################################################
  ##############################################################################
  run_ofFromObj <- function(blast00,
                            blast01,
                            blast10,
                            blast11,
                            pep0,
                            pep1,
                            writeDir){

    ofID <- ofID1 <- ofID2 <- id <- ofID <- og <- Orthogroup <- NULL
    if(dir.exists(writeDir))
      stop(sprintf("%s exists. Specify non-existing directory\n",
                   writeDir))
    dir.create(writeDir)

    if(colnames(blast00)[1] != "ofID1"){
      setnames(blast00, 1:2, c("ofID1", "ofID2"))
      setnames(blast01, 1:2, c("ofID1", "ofID2"))
      setnames(blast10, 1:2, c("ofID1", "ofID2"))
      setnames(blast11, 1:2, c("ofID1", "ofID2"))
    }

    # -- make gene ID dictionaries
    id0 <- unique(c(blast00$ofID1, blast00$ofID2, blast01$ofID1, blast10$ofID2))
    id1 <- unique(c(blast11$ofID1, blast11$ofID2, blast01$ofID2, blast10$ofID1))
    id0 <- id0[order(id0)]
    id1 <- id1[order(id1)]
    names(id0) <- sprintf("0_%s", (1:length(id0))-1)
    names(id1) <- sprintf("1_%s", (1:length(id1))-1)

    # -- ensure that all ids are in the peptide files
    id0 <- id0[id0 %in% names(pep0)]
    id1 <- id1[id1 %in% names(pep1)]

    # -- rename peptides and invert dictionary
    p0 <- pep0[id0]; names(p0) <- names(id0)
    p1 <- pep1[id1]; names(p1) <- names(id1)
    di0 <- names(id0); di1 <- names(id1); names(di0) <- id0; names(di1) <- id1

    # -- write the peptide files / fake diamond dbs
    writeXStringSet(p0, filepath = file.path(writeDir, "Species0.fa"))
    writeXStringSet(p1, filepath = file.path(writeDir, "Species1.fa"))
    cat("NA", file = file.path(writeDir, "diamondDBSpecies0.dmnd"))
    cat("NA", file = file.path(writeDir, "diamondDBSpecies1.dmnd"))

    # -- write the species and sequence IDs
    sid <- data.table(
      of = paste0(c(names(id0), names(id1)), ":"),
      id = c(id0, id1))
    fwrite(
      sid, file = file.path(writeDir, "SequenceIDs.txt"),
      sep = " ", quote = F, row.names = F, col.names = F)
    cat(
      c("0: species1.fa", "1: species2.fa"),
      sep = "\n", file = file.path(writeDir, "SpeciesIDs.txt"))

    # -- rename the blast files
    bl00 <- subset(blast00, ofID1 %in% names(di0) & ofID2 %in% names(di0))
    bl01 <- subset(blast01, ofID1 %in% names(di0) & ofID2 %in% names(di1))
    bl10 <- subset(blast10, ofID1 %in% names(di1) & ofID2 %in% names(di0))
    bl11 <- subset(blast11, ofID1 %in% names(di1) & ofID2 %in% names(di1))
    bl00[,`:=`(ofID1 = di0[ofID1], ofID2 = di0[ofID2])]
    bl01[,`:=`(ofID1 = di0[ofID1], ofID2 = di1[ofID2])]
    bl10[,`:=`(ofID1 = di1[ofID1], ofID2 = di0[ofID2])]
    bl11[,`:=`(ofID1 = di1[ofID1], ofID2 = di1[ofID2])]
    bl00 <- subset(bl00, complete.cases(bl00[,1:12]))[,1:12]
    bl01 <- subset(bl01, complete.cases(bl01[,1:12]))[,1:12]
    bl10 <- subset(bl10, complete.cases(bl10[,1:12]))[,1:12]
    bl11 <- subset(bl11, complete.cases(bl11[,1:12]))[,1:12]

    # -- write the blasts
    fwrite(
      bl00, file = file.path(writeDir, "Blast0_0.txt.gz"),
      sep = "\t", quote = F, row.names = F, col.names = F)
    fwrite(
      bl01, file = file.path(writeDir, "Blast0_1.txt.gz"),
      sep = "\t", quote = F, row.names = F, col.names = F)
    fwrite(
      bl10, file = file.path(writeDir, "Blast1_0.txt.gz"),
      sep = "\t", quote = F, row.names = F, col.names = F)
    fwrite(
      bl11, file = file.path(writeDir, "Blast1_1.txt.gz"),
      sep = "\t", quote = F, row.names = F, col.names = F)

    # -- run orthofinder
    com <- sprintf("orthofinder -b %s -og -a 1 -t 1  1>/dev/null 2>&1", writeDir)
    system(com)

    # -- find the files
    ogf <- order_filesByMtime(
      path = writeDir,
      recursive = T,
      pattern = "Orthogroups.tsv")[1]

    # -- read the orthogroups.tsv file and process
    ogdt <- fread(ogf, showProgress = F, verbose = F)
    ogdt <- melt(
      ogdt, id.vars = "Orthogroup", variable.name = "genome", value.name = "id")
    ogdt <- ogdt[,list(id = strsplit(id, ",")[[1]]), by = c("Orthogroup", "genome")]
    ogdt[,`:=`(genome = NULL, ofID = trimws(id), id = NULL,
               og = trimws(Orthogroup), Orthogroup = NULL)]
    ogdt <- subset(ogdt, !duplicated(ogdt))
    hasDup <- subset(ogdt, ofID %in% subset(ogdt, duplicated(ofID))$ofID)
    if(nrow(hasDup) > 1){
      m <- merge(hasDup, hasDup, by = "ofID", all = T, allow.cartesian = T)
      ci <- clus_igraph(m$og.x, m$og.y)
      ci <- ci[!duplicated(names(ci))]
      ogdt[,og := ifelse(og %in% names(ci), ci[og], og)]
      ogdt <- subset(ogdt, !duplicated(ogdt))
    }
    ogdt[,og := as.integer(factor(og, unique(og)))]

    nog <- c(id0, id0)
    nog <- nog[!nog %in% ogdt$ofID]

    # -- return data.table of ogs
    if(length(nog) > 0){
      ogdt <- rbind(ogdt, data.table(
        ofID = nog, og = sprintf("NOG_%s",1:length(nog))))
    }

    return(ogdt)
  }
  ##############################################################################
  ##############################################################################
  count_expectn <- function(hits, gff){
    ord1 <- ord2 <- genome <- chr <- start <- end <- ofID <- nReg <- regID <- NULL
    b <- hits[,list(start = min(ord1), end = max(ord1)),
              by = c("gen1","chr1","regID")]
    if(hits$gen1[1] != hits$gen2[1] ){
      b2 <- hits[,list(start = min(ord2), end = max(ord2)),
                 by = c("gen2","chr2","regID")]
      setnames(b2, colnames(b))
      b <- rbind(b, b2)
    }
    b <- with(b, data.table(
      genome = gen1, chr = chr1, regID = regID, start = start, end = end))
    g <- with(subset(gff, genome %in% c(hits$gen1[1], hits$gen2[1])), data.table(
      genome = genome, chr = chr, start = ord, end = ord, ofID = ofID))
    setkey(g, genome, chr, start, end)
    setkey(b, genome, chr, start, end)
    f <- subset(foverlaps(g, b), !is.na(ofID))[,c("genome","regID","ofID")]
    f[,nReg := uniqueN(regID[!is.na(regID)]), by = "ofID"]
    return(f)
  }
  ##############################################################################
  ##############################################################################
  genome1 <- genome2 <- gnum1 <- gnum2 <- lgRegionID <- regBuffer <- n2 <- NULL
  n1 <- V2 <- V1 <- blkID <- hasOgGlob <- og <- regID <- ofID <- regAnchor <- NULL
  ofID1 <- ofID2 <- synOg1 <- synOg2 <- hasOgPw <- NULL
  ##############################################################################
  # 1.Checking
  ##############################################################################
  # -- get various parameters
  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  verbose <- gsParam$params$verbose
  nCores <- gsParam$params$nCores
  # -- make sure orthofinder has been run and blast results exist
  if(is.na(gsParam$paths$blastDir))
    gsParam <- find_orthofinderResults(gsParam)

  # -- read in gff
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  # -- get the synParams in
  synp <- data.table(gsParam$params$synteny)
  synp[,`:=`(gnum1 = match(genome1, genomeIDs),
             gnum2 = match(genome2, genomeIDs))]
  synp <- subset(synp, gnum1 <= gnum2)
  setkey(synp, gnum1, gnum2)

  # -- get orthofinder species IDs
  ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)

  # -- read in and split up the peptide files
  pepspl <- sapply(genomeIDs, USE.NAMES = T, simplify = F, function(i)
    readAAStringSet(file.path(gsParam$paths$blastDir,
                              sprintf("Species%s.fa", ofSpId[i]))))

  # -- run orthofinder within each region
  p5 <- function(x)
    paste(c(x, rep(" ", max(0, 7 - nchar(x)))), collapse = "")
  if(verbose)
    cat(sprintf(
      "\tRunning orthofinder by region ...\n\tGenome1-Genome2| Copy Number: %s%s%s%s\n",
      p5("absent"),p5("1x"),p5("2x"),p5("2+x")))

  blnk <- paste(rep(" ", 16), collapse = "")
  synOgInBlkHits <- rbindlist(lapply(1:nrow(synp), function(i){
    geno1 <- synp$genome1[i]
    geno2 <- synp$genome2[i]
    if(verbose)
      cat(sprintf("\t%s-%s|", pull_strWidth(geno1, 7),pull_strWidth(geno2, 7)))

    # -- load the syntenic hits
    fs <- file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz",geno1, geno2))
    hits <- subset(fread(
      fs,
      na.strings = c("NA", ""),
      showProgress = F),
      !is.na(lgRegionID) & regBuffer)

    # -- calculate block coordinates from lgRegions (aggregated regions)
    hits[,`:=`(blkID = lgRegionID, blkAnchor = regAnchor, blkBuffer = regBuffer)]
    bv <- with(hits, c(blkID, blkID))
    names(bv) <- with(hits, c(paste(ofID1, ofID2), paste(ofID2, ofID1)))
    bv <- bv[!duplicated(names(bv))]
    expn <- count_expectn(hits = hits, gff = gff)
    tb <- subset(expn, !duplicated(paste(ofID, regID)))
    tb <- tb[,list(n = sum(!is.na(regID))), by = "ofID"]
    if(verbose)
      cat(sprintf(" expected CN: %s%s%s%s\n",
                  p5(sum(tb$n == 0)), p5(sum(tb$n == 1)), p5(sum(tb$n == 2)), p5(sum(tb$n > 2))))
    oggl <- with(subset(hits, !is.na(og)), unique(paste(
      c(regID, regID), c(ofID1, ofID2))))
    expn[,hasOgGlob := paste(regID, ofID) %in% oggl]
    tb <- expn[,list(n = sum(!is.na(regID) & hasOgGlob)), by = "ofID"]
    if(verbose)
      cat(sprintf("\t%sglob. hit CN: %s%s%s%s\n",
                  blnk,  p5(sum(tb$n == 0)), p5(sum(tb$n == 1)), p5(sum(tb$n == 2)), p5(sum(tb$n > 2))))

    # -- read raw blast hits
    h <- read_hits4of(
      gsParam = gsParam,
      genome1 = geno1,
      genome2 = geno2)
    h[,blkID := bv[paste(V1, V2)]]
    h <- subset(h, !is.na(blkID))
    hc <- h[,list(n1 = uniqueN(V1), n2 = uniqueN(V2)), by = "blkID"]
    hc <- subset(hc, n1 >= minGenes4of & n2 >= minGenes4of)
    setorder(hc, -n1, -n2)
    # -- split hits by blocks with enuf genes
    splh <- sapply(split(h, by = "blkID"), simplify = F, USE.NAMES = T, function(x) x[,1:12])

    # -- run of in each block
    synOgHits <- rbindlist(mclapply(hc$blkID, mc.cores = nCores, function(j){
      hj <- splh[[j]]
      rj <- subset(hits, lgRegionID == j)
      tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4inBlkTMPdir", j))
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      ogdt <- run_ofFromObj(
        blast00 = subset(h, V1 %in% rj$ofID1 & V2 %in% rj$ofID1),
        blast01 = subset(h, V1 %in% rj$ofID1 & V2 %in% rj$ofID2),
        blast10 = subset(h, V1 %in% rj$ofID2 & V2 %in% rj$ofID1),
        blast11 = subset(h, V1 %in% rj$ofID2 & V2 %in% rj$ofID2),
        pep0 = pepspl[[geno1]],
        pep1 = pepspl[[geno2]],
        writeDir = tmpDir)
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      ogv <- ogdt$og; names(ogv) <- ogdt$ofID
      rj[,`:=`(synOg1 = ogv[ofID1], synOg2 = ogv[ofID2])]
      return(subset(rj, !is.na(og) | synOg1 == synOg2))
    }))

    ogpw <- with(subset(synOgHits, synOg1 == synOg2), unique(paste(
      c(regID, regID), c(ofID1, ofID2))))

    expn[,hasOgPw := paste(regID, ofID) %in% ogpw]
    tb <- expn[,list(n = sum(!is.na(regID) & hasOgPw)), by = "ofID"]
    if(verbose)
      cat(sprintf("\t%sinblk hit CN: %s%s%s%s\n",
                  blnk,  p5(sum(tb$n == 0)), p5(sum(tb$n == 1)), p5(sum(tb$n == 2)), p5(sum(tb$n > 2))))

    return(with(subset(synOgHits, synOg1 == synOg2), data.table(
      ofID1 = ofID1, ofID2 = ofID2, synOg = synOg1)))
  }))

  return(synOgInBlkHits)
}
