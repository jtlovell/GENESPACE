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
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param genomeIDs character vector with the genomes to include in the run
#' @param minGenes4of integer specifying the minimum number of genes needed to
#' run orthofinder.
#' @param quietOrthofinder logical, should orthofinder-generated output be
#' printed ot the console?
#' @param blast00 data.table containing the blast hits of genome1 to genome1
#' @param blast01 data.table containing the blast hits of genome1 to genome2
#' @param blast10 data.table containing the blast hits of genome2 to genome1
#' @param blast11 data.table containing the blast hits of genome2 to genome2
#' @param pep0 aastringset containing peptides for genome1
#' @param pep1 aastringset containing peptides for genome2
#' @param writeDir file path pointing to a directory within which to run f
#' @param genome1 character string specifying genome1
#' @param genome2 character string specifying genome2
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
                            overwrite = FALSE,
                            quietOrthofinder = FALSE){

  ##############################################################################
  ##############################################################################
  ##############################################################################
  default_ofDb <- function(gsParam, quietOrthofinder){
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
    quiet <- ifelse(quietOrthofinder, "1>/dev/null 2>&1", "")
    com <- sprintf(
      "%s -f %s -t %s -a 1 -X -o %s %s",
      p2of,
      dirname(gsParam$paths$peptide[1]),
      gsParam$params$nCores,
      gsParam$paths$orthofinder,
      quiet)

    ############################################################################
    # 4. run it
    if(dontRun){
      cat(
        "\tCould not find valid orthofinder executable in the path\n",
        "\tRun the following command outside of R (assuming orthofinder is in the path):",
        "\n################\n",
        sprintf("cd %s\n", dirname(gsParam$paths$orthofinder)),
        com,
        "\n################\n", sep = "")
    }else{
      system(com)
    }

    return(com)
  }

  ##############################################################################
  # --- invert blast file
  invert_blast <- function(fileIn, fileOut){
    tmp <- fread(fileIn, verbose = F, showProgress = F)
    tmp <- tmp[,c(2,1,3:6,8,7,10,9,11,12)]
    fwrite(
      tmp, sep = "\t", quote = FALSE,  col.names = FALSE, row.names = FALSE,
      file = fileOut, showProgress = FALSE, verbose = FALSE)
  }

  ##############################################################################
  # -- move and reorganize orthofinder input files
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

  ##############################################################################
  # -- add blast metadata / calls
  add_blastInfo2syn <- function(gsParam){
    genome1 <- genome2 <- db2 <- fa1 <- blFile <- NULL
    si <- read_orthofinderSpeciesIDs(gsParam$paths$orthofinder)
    p <- data.table(gsParam$params$synteny)
    diamondMode <- gsParam$params$diamondMode
    ofd <- gsParam$paths$orthofinder
    p[,`:=`(
      db1 = file.path(ofd, sprintf("diamondDBSpecies%s.dmnd", si[genome1])),
      db2 = file.path(ofd, sprintf("diamondDBSpecies%s.dmnd", si[genome2])),
      fa1 = file.path(ofd, sprintf("Species%s.fa", si[genome1])),
      fa2 = file.path(ofd, sprintf("Species%s.fa", si[genome2])),
      blFile = file.path(ofd, sprintf("Blast%s_%s.txt.gz",
                                      si[genome1], si[genome2])),
      invertFile = file.path(ofd, sprintf("Blast%s_%s.txt.gz",
                                          si[genome2], si[genome1])))]
    dm <- ifelse(
      gsParam$params$diamondMode == "--fast", "", gsParam$params$diamondMode)
    p[,com := sprintf(
      "%s blastp %s --quiet -e %s -p %s --compress 1 -d %s -q %s -o %s",
      "diamond", dm, .1, gsParam$params$nCores, db2, fa1, blFile)]
    return(p)
  }

  ##############################################################################
  # -- fast pairwise (non-reciprocal) blast hits
  fast_ofDb <- function(gsParam, quietOrthofinder){

    # -- Remove existing orthofinder directory if it exists
    if(dir.exists(gsParam$paths$orthofinder))
      unlink(gsParam$paths$orthofinder, recursive = T)

    # -- convert to orthofinder
    com <- sprintf(
      "%s -f %s -t %s -a 1 -op -o %s 1>/dev/null 2>&1",
      gsParam$paths$orthofinderCall,
      dirname(gsParam$paths$peptide[1]),
      gsParam$params$nCores,
      gsParam$paths$orthofinder)
    system(com)

    # -- place orthofinder input files in paths$orthofinder
    reorg_ofInput(gsParam$paths$orthofinder)

    # -- get blast parameters
    p <- add_blastInfo2syn(gsParam = gsParam)

    # -- Run blasts
    runBlast <- NULL
    pwp <- subset(p, runBlast)
    for(i in 1:nrow(pwp)){
      if(gsParam$params$verbose)
        with(pwp[i,], cat(sprintf(
          "\t\tRunning %s/%s (%s vs. %s)\n",
          i, nrow(pwp), genome1, genome2)))
      system(pwp$com[i])
    }

    # -- Invert blasts if necessary
    if(gsParam$params$verbose)
      cat("\t\tDone!\n\tInverting intergenomic files ... ")
    runBlast <- genome1 <- genome2 <- NULL
    ip <- subset(p, !runBlast & genome1 != genome2)
    for(i in 1:nrow(ip))
      invert_blast(fileIn = ip$invertFile[i], fileOut = ip$blFile[i])

    # -- Run orthofinder
    if(gsParam$params$verbose)
      cat("Done!\n\tRunning orthofinder -og on pre-computed blast:\n")
    quiet <- ifelse(quietOrthofinder, "1>/dev/null 2>&1", "")
    com <- with(gsParam, sprintf(
      "%s -b %s -t %s -a 1 -X -og %s",
      paths$orthofinderCall,
      paths$orthofinder,
      params$nCores,
      quiet))
    system(com)
    return(com)
  }

  ##############################################################################
  ##############################################################################
  # -- check the orthofinder call
  if(is.logical(gsParam$paths$orthofinderCall))
    if(!gsParam$paths$orthofinderCall && !is.na(gsParam$paths$orthofinderCall))
      gsParam$paths$orthofinderCall <- NA

  # -- set the synteny parameters
  if(is.data.table(gsParam$params$synteny))
    if(!all(gsParam$genomes$genomeIDs %in% gsParam$params$synteny$genome1))
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
      com <- default_ofDb(
        gsParam,
        quietOrthofinder = quietOrthofinder)
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
        com <- fast_ofDb(
          gsParam,
          quietOrthofinder = quietOrthofinder)
      }else{
        if(gsParam$params$verbose)
          cat("\tRunning 'defualt' genespace orthofinder method",
              "\n\t############################################################\n")
        com <- default_ofDb(
          gsParam,
          quietOrthofinder = quietOrthofinder)
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
#' @importFrom Biostrings readAAStringSet
#' @export
blkwise_orthofinder <- function(gsParam,
                                gff,
                                overwrite = FALSE,
                                genomeIDs = NULL,
                                minGenes4of = 40){
  setDTthreads(1)
  inblkOG <- NULL
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
  genome1 <- genome2 <- runBlast <- gnum1 <- gnum2 <- NULL
  synp <- data.table(gsParam$params$synteny)
  synp[,`:=`(gnum1 = match(genome1, genomeIDs),
             gnum2 = match(genome2, genomeIDs))]
  synp <- subset(synp, runBlast)
  setkey(synp, gnum1, gnum2)

  # -- get orthofinder species IDs
  ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)

  # -- read in and split up the peptide files
  pepspl <- sapply(genomeIDs, USE.NAMES = T, simplify = F, function(i)
    readAAStringSet(file.path(gsParam$paths$blastDir,
                              sprintf("Species%s.fa", ofSpId[i]))))

  # -- run orthofinder within each region
  if(verbose)
    cat("Running orthofinder by region ... \n\tgenome combinat. : n. non-self genes, nOGs global/syntenic/inblk\n")
  genome <- chr <- start <- end <- isArrayRep <- ofID <- genome <- NULL
  setkey(gff, genome, chr, start, end)

  arrep <- with(subset(gff, isArrayRep), split(ofID, genome))

  synOgInBlkHits <- rbindlist(lapply(1:nrow(synp), function(i){

    geno1 <- synp$genome1[i]
    geno2 <- synp$genome2[i]

    outf <- file.path(
      gpar$paths$results,
      sprintf("%s_%s_inblkOGs.txt.gz",  geno1, geno2))

    if(verbose)
      cat(sprintf(
        "\t%s-%s: ",
        pull_strWidth(geno1, 8), pull_strWidth(geno2, 8)))

    # -- load the syntenic hits
    fs <- file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz", geno1, geno2))
    hits <- fread(fs, na.strings = c("NA", ""), showProgress = F,
                  select = c("ofID1","ofID2", "regID","inBuffer"))

    # -- subset to hits in the large regions
    regID <- inBuffer <- NULL
    hits <- subset(hits, !is.na(regID) & inBuffer)

    # -- drop self region hits
    isSelf <- ofID1 <- ofID2 <- NULL
    hits[,isSelf := any(ofID1 == ofID2), by = "regID"]
    hits <- subset(hits, !isSelf)

    # --  subset hits to array Reps and drop regions smaller than min genes
    ofID1 <- ofID2 <- n1 <- n2 <- genome <- NULL
    hits <- subset(hits, ofID1 %in% arrep[[geno1]] & ofID2 %in% arrep[[geno2]])

    hits[,`:=`(n1 = uniqueN(ofID1), n2 = uniqueN(ofID2)), by = "regID"]
    hits <- subset(hits, n1 >= minGenes4of & n2 >= minGenes4of)

    g <- subset(gff, genome %in% c(geno1, geno2))

    if(file.exists(outf) && !overwrite){
      inblkOgDt <- fread(outf, showProgress = F, na.strings = c("NA", ""))
      if(verbose)
        cat("using pre-calc. data ... ")
    }else{
      if(nrow(hits) < minGenes4of){
        if(verbose)
          cat("no non-self syn. regions\n")
        return(NULL)
      }else{
        # -- subset to the genomes of interest and report updates
        ofID1 <- ofID2 <- arrep <- V1 <- V2 <- genome <- ofID <- synOG <-
          globOG <- NULL

        # -- read in the hits
        h <- read_hits4of(
          gsParam = gsParam,
          genome1 = geno1,
          genome2 = geno2)
        u12 <- with(hits, paste(ofID1, ofID2))
        u21 <- with(hits, paste(ofID2, ofID1))
        h00 <- subset(h, V1 %in% arrep[[geno1]] & V1 == V2)
        h11 <- subset(h, V1 %in% arrep[[geno2]] & V1 == V2)
        h01 <- subset(h, paste(V1, V2) %in% u12)
        h10 <- subset(h, paste(V1, V2) %in% u21)
        hspl <- split(hits, by = "regID")

        # -- split hits by lgRegs
        inblkOgDt <- rbindlist(mclapply(names(hspl), mc.cores = nCores, mc.preschedule = F, function(j){
          tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4inBlkTMPdir", j))
          if(dir.exists(tmpDir))
            unlink(tmpDir, recursive = T)
          on.exit(expr = unlink(tmpDir, recursive = T))

          out <- data.table(hspl[[j]])
          u1 <- unique(out$ofID1)
          u2 <- unique(out$ofID2)

          # -- run orthofinder from these hits
          V1 <- V2 <- regID <- og <- isInblkOg <- ofID1 <- ofID2 <- NULL
          ogdt <- run_ofFromObj(
            blast00 = subset(h00, V1 %in% u1 & V2 %in% u1),
            blast01 = subset(h01, V1 %in% u1 & V2 %in% u2),
            blast10 = subset(h10, V1 %in% u2 & V2 %in% u1),
            blast11 = subset(h11, V1 %in% u2 & V2 %in% u2),
            pep0 = pepspl[[geno1]],
            pep1 = pepspl[[geno2]],
            writeDir = tmpDir)
          ogdt[,og := as.numeric(as.factor(og))]
          ogv <- ogdt$og; names(ogv) <- ogdt$ofID
          u <- unique(c(u1, u2))
          uo <- u[!u %in% names(ogv)]
          uv <- (max(ogv) + 1):(max(ogv) + length(uo))
          names(uv) <- uo
          ogv <- c(ogv, uv)
          out[,isInblkOg := ogv[ofID1] == ogv[ofID2]]
          unlink(tmpDir, recursive = T)
          return(out[,c(1:2,8)])
        }))

        ofID1 <- ofID2 <- ofID <- NULL
        ic <- with(subset(inblkOgDt, isInblkOg), clus_igraph(
          id1 = c(ofID1, ofID2), id2 = c(ofID2, ofID1)))
        ic <- ic[!duplicated(names(ic))]
        uc <- with(hits, unique(c(ofID1, ofID2)))
        uc <- uc[!uc %in% names(ic)]

        fwrite(inblkOgDt, file = outf, sep = "\t",
               quote = FALSE, showProgress = FALSE)
      }
      if(verbose)
        with(subset(g, ofID %in% c(inblkOgDt$ofID1, inblkOgDt$ofID2)), cat(
          sprintf("%s genes, %s / %s ",
                  uniqueN(ofID), uniqueN(globOG), uniqueN(synOG))))
      if(verbose)
        cat(sprintf("/ %s\n", uniqueN(ic) + length(uc)))
      return(inblkOgDt)
    }
  }))

  isInblkOg <- arrayID <- a1 <- a2 <- ofID1 <- ofID2 <- inBlkOG <- arrv <- NULL
  ib <- subset(synOgInBlkHits, isInblkOg)
  gff[,arrv := as.character(as.numeric(as.factor(arrayID)))]
  av <- gff$arrv; names(av) <- gff$ofID
  ib[,`:=`(a1 = as.character(av[ofID1]),
           a2 = as.character(av[ofID2]))]
  ic <- with(ib, clus_igraph(id1 = a1, id2 = a2))
  ic <- ic[!duplicated(names(ic))]
  icn <- gff$arrv[!gff$arrv %in% names(ic)]
  icv <- max(ic + 1):(max(ic) + length(icn))
  names(icv) <- icn
  ic <- c(ic, icv)
  gff[,inblkOG := ic[arrv]]
  gff[,arrv := NULL]
  return(gff)
}

#' @title run_ofFromObj
#' @description
#' \code{run_ofFromObj} run_ofFromObj
#' @rdname run_orthofinder
#' @import data.table
#' @importFrom Biostrings writeXStringSet
#' @export
run_ofFromObj <- function(blast00,
                          blast01,
                          blast10,
                          blast11,
                          pep0,
                          pep1,
                          writeDir){

  id <- Orthogroup <- ofID <- og <- NULL
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
  di0 <- names(id0)
  di1 <- names(id1)
  names(di0) <- id0
  names(di1) <- id1

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
  ofID1 <- ofID2 <- NULL
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

#' @title read_hits4of
#' @description
#' \code{read_hits4of} read_hits4of
#' @rdname run_orthofinder
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
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
