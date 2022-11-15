#' @title Set synteny parameters and parse required data
#' @description
#' \code{set_syntenyParams} Generate all data needed to run synteny. This
#' includes the synteny parameters, combined/annotated bed file and annotated
#' blast files.
#' @name set_syntenyParams
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param resultsDir file.path to the results directory
#' @param orthofinderDir file.path to the raw orthofinder directory
#'
#' \cr
#' If called, \code{set_syntenyParams} returns its own arguments.
#'
#' @title find orthofinder results files in the results directory
#' @description
#' \code{set_syntenyParams} Locate the required orthofinder results for genespace,
#' which are stored in the /results directory. Typically for internal use. If
#' called directly, use with caution.
#' @rdname set_syntenyParams
#' @import data.table
#' @importFrom utils combn
#' @export
set_syntenyParams <- function(gsParam){

  query <- target <- tmp1 <- genNum1 <- genNum2 <- tmp2 <- query1 <-
    fs1 <- fs2 <- queryBlast <- synHits <- synBuff <- NULL

  resultsDir <- gsParam$paths$results
  synHitsDir <- gsParam$paths$syntenicHits
  genomeIDs  <- gsParam$genomeIDs
  ploidy     <- gsParam$ploidy

  if(!is.na(gsParam$outgroup)){
    gids <- c(genomeIDs, gsParam$outgroup)
    gids <- gids[!duplicated(gids)]
  }else{
    gids <- genomeIDs
  }

  ##############################################################################
  # 1. check that the results directory exists
  runPass <- TRUE
  if(!dir.exists(resultsDir))
    runPass <- FALSE

  ##############################################################################
  # 2. check that all the required files exist and that
  if(runPass){
    fs <- file.path(resultsDir, c(
      "SpeciesIDs.txt", "SequenceIDs.txt", "Orthogroups.tsv",
      "N0.tsv", "SpeciesTree_rooted.txt"))
    names(fs) <- c("SpeciesIDs", "SequenceIDs", "ogs", "hogs", "speciesTree")
    if(!all(file.exists(fs[1:3])))
      runPass <- FALSE
  }

  ##############################################################################
  # 3. check that the genomeIDs are all represented in the speciesIDs.txt
  if(runPass){
    sids <- read_orthofinderSpeciesIDs(fs["SpeciesIDs"])
    chk <- all(gids %in% names(sids)) && all(names(sids) %in% gids)
    if(!chk)
      runPass <- FALSE
  }

  ##############################################################################
  # 4. check that all the genomes exist in blast files
  if(runPass){
    # -- make all combinations of genomes
    blMd <- data.table(CJ(query = genomeIDs, target = genomeIDs))
    blMd[,`:=`(queryPloidy = ploidy[query],
               targetPloidy = ploidy[target],
               genNum1 = sids[query],
               genNum2 = sids[target])]

    # -- get file paths to the two blast files
    blMd[,tmp1 := file.path(resultsDir, sprintf(
      "Blast%s_%s.txt.gz", genNum1, genNum2))]
    blMd[,tmp2 := file.path(resultsDir, sprintf(
      "Blast%s_%s.txt.gz", genNum2, genNum1))]

    # -- check that files exist
    blMd[,`:=`(ex1 = file.exists(tmp1), ex2 = file.exists(tmp2))]

    # -- if only one of the files exists, copy to other
    blMd$tmp1[with(blMd, !ex1 & ex2)] <- blMd$tmp2[with(blMd, !ex1 & ex2)]
    blMd$tmp2[with(blMd, !ex2 & ex1)] <- blMd$tmp1[with(blMd, !ex2 & ex1)]

    # -- get file sizes for blast files
    blMd[,`:=`(fs1 = file.size(tmp1), fs2 = file.size(tmp2))]

    # -- if files missing, file.size returns NA, so drop missing files
    blMd <- subset(blMd, complete.cases(blMd))

    # -- ensure that all genome combinations are in there
    cmb <- combn(genomeIDs, 2, simplify = F)
    u <- with(blMd, paste(query, target))
    chk <- all(sapply(cmb, function(x)
      paste(x[1], x[2]) %in% u || paste(x[2], x[1]) %in% u))
    if(!chk)
      runPass <- FALSE
  }

  ##############################################################################
  # 5. Get the query and target genomes figured out
  if(runPass){
    # -- query genome is the one with the larger file as the query
    blMd[,query1 := fs1 >= fs2]
    blMd[,`:=`(queryBlast = ifelse(query1, tmp1, tmp2),
               targetBlast = ifelse(target == query, NA,
                                    ifelse(query1, tmp2, tmp1)))]
    blMd <- subset(blMd, queryBlast == tmp1)
    blMd <- blMd[,c("query", "target", "queryPloidy", "targetPloidy",
                    "queryBlast", "targetBlast")]

    cmb <- combn(genomeIDs, 2, simplify = F)
    u <- with(blMd, paste(query, target))
    chk <- all(sapply(cmb, function(x)
      paste(x[1], x[2]) %in% u || paste(x[2], x[1]) %in% u))
    if(!chk)
      runPass <- FALSE
  }

  ##############################################################################
  # 6. Check and add the ortholog files
  if(runPass){
    blMd[,`:=`(
      queryOrthologs = file.path(
        resultsDir, sprintf("%s__v__%s.tsv", query, target)),
      targetOrthologs = file.path(
        resultsDir, sprintf("%s__v__%s.tsv", target, query)))]
    blMd$queryOrthologs[!file.exists(blMd$queryOrthologs)] <- NA
    blMd$targetOrthologs[!file.exists(blMd$targetOrthologs)] <- NA
    blMd$targetOrthologs[with(blMd, queryOrthologs == targetOrthologs)] <- NA
  }

  ##############################################################################
  # 7. Add the synHits file paths
  if(runPass)
    blMd[,synHits := file.path(
      synHitsDir, sprintf("%s_vs_%s.synBlast.txt.gz", query, target))]

  ##############################################################################
  # 8. Add the parameters to the blMd
  if(runPass){
    pnames <- c(
      "orthofinderInBlk", "blkSize", "blkRadius", "nGaps", "synBuff",
      "onlyOgAnchors", "nSecondaryHits", "blkSizeSecond", "blkRadiusSecond",
      "nGapsSecond", "onlyOgAnchorsSecond")
    for(i in pnames)
      blMd[[i]] <- gsParam$params[[i]]

    blMd <- data.table(blMd)
    blMd[, `:=`(synRad = sqrt(2) * synBuff,
                inBufferRadius = synBuff/2)]
  }

  ##############################################################################
  # 9 return the list
  if(runPass)
    gsParam$synteny <- list(
      combBed = file.path(resultsDir, "combBed.txt"),
      SpeciesIDs = file.path(resultsDir, "SpeciesIDs.txt"),
      SequenceIDs = file.path(resultsDir, "SequenceIDs.txt"),
      ogs = file.path(resultsDir, "Orthogroups.tsv"),
      hogs = file.path(resultsDir, "N0.tsv"),
      speciesTree = file.path(resultsDir, "SpeciesTree_rooted.txt"),
      blast = blMd)
  if(!file.exists(file.path(resultsDir, "N0.tsv")) & runPass)
    gsParam$synteny$hogs <- NA
  if(!file.exists(file.path(resultsDir, "SpeciesTree_rooted.txt")) & runPass)
    gsParam$synteny$speciesTree <- NA
  return(gsParam)
}

#' @title find orthofinder results files
#' @description
#' \code{find_ofFiles} Find raw orthofinder results files. Typically for
#' internal use with copy_of2results. If called directly, use with caution.
#' @rdname set_syntenyParams
#' @import data.table
#' @export
find_ofFiles <- function(orthofinderDir){

  ##############################################################################
  # 1. check to ensure that a valid orthofinder run exists
  ofDir <- check_filePathParam(orthofinderDir)
  if(is.na(ofDir))
    stop("could not coerce", ofDir, "to a valid file.path\n")

  ofDir <- list.files(
    path = ofDir, pattern = "Results_", full.names = TRUE)

  # -- stop if there is no run
  if(length(ofDir) > 1)
    stop("found multiple orthofinder runs in ",
         orthofinderDir, "\n\tremove all but one.\n")
  if(length(ofDir) != 1){
    hasRun <- FALSE
  }else{
    hasRun <- TRUE
  }

  ofPaths <- list()

  ##############################################################################
  # 2. check that the working directory is in ofDir
  ofWd <- file.path(ofDir, "WorkingDirectory")
  if(!dir.exists(ofWd)){
    hasRun <- FALSE
  }

  ##############################################################################
  # 3. species IDs
  if(hasRun){
    tmp <- file.path(ofWd, "SpeciesIDs.txt")
    if(file.exists(tmp)){
      ofPaths$SpeciesIDs <- tmp
      sids <- read_orthofinderSpeciesIDs(ofPaths$SpeciesIDs)
    }else{
      hasRun <- FALSE
    }
  }

  ##############################################################################
  # 4. sequence IDs
  if(hasRun){
    tmp <- file.path(ofWd, "SequenceIDs.txt")
    if(file.exists(tmp)){
      ofPaths$SequenceIDs <- tmp
    }else{
      hasRun <- FALSE
    }
  }

  ##############################################################################
  # 5. OGs
  if(hasRun){
    tmp <- file.path(ofDir, "Orthogroups", "Orthogroups.tsv")
    if(file.exists(tmp)){
      ofPaths$ogs <- tmp
    }else{
      hasRun <- FALSE
    }
  }

  ##############################################################################
  # 6. HOGs
  if(hasRun){
    tmp <- file.path(ofDir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
    ofPaths$hogs <- tmp
  }

  ##############################################################################
  # 7. species tree
  if(hasRun){
    tmp <- file.path(ofDir, "Species_Tree", "SpeciesTree_rooted.txt")
    if(file.exists(tmp)){
      ofPaths$speciesTree <- tmp
    }else{
      ofPaths$speciesTree <- NA
    }
  }

  ##############################################################################
  # 8. blast files
  if(hasRun){
    if(!is.null(sids)){

      genome1 <- genome2 <- blastFile <- genNum1 <- genNum2 <-
        orthologFile <-NULL

      # -- get combination of genomes
      genomeIDs <- names(sids)
      blMd <- data.table(CJ(genome1 = genomeIDs, genome2 = genomeIDs))

      # -- get species ID numbers
      blMd[,`:=`(genNum1 = sids[genome1], genNum2 = sids[genome2])]

      # -- get expected blast file locations
      blMd[,blastFile := file.path(ofWd, sprintf(
        "Blast%s_%s.txt.gz", genNum1, genNum2))]

      # -- get ortholog files
      blMd[,orthologFile := file.path(
        ofDir,
        "Orthologues",
        sprintf("Orthologues_%s", genome1),
        sprintf("%s__v__%s.tsv", genome1, genome2))]

      # -- set missing files to NA
      blMd$blastFile[!file.exists(blMd$blastFile)] <- NA
      blMd$orthologFile[!file.exists(blMd$orthologFile)] <- NA

      # -- save to ofPaths to return
      ofPaths$blast <- blMd
    }
  }

  if(!hasRun)
    ofPaths <- NA
  return(ofPaths)
}


#' @title copy raw orthofinder results to the genespace directory
#' @description
#' \code{copy_of2results} copy files from an orthofinder run to the genespace
#' results directory. Typically for internal use. If called directly, use with
#' caution.
#' @rdname set_syntenyParams
#' @import data.table
#' @export
copy_of2results <- function(orthofinderDir, resultsDir){
  ofDir   <- check_filePathParam(orthofinderDir)
  resDir  <- check_filePathParam(resultsDir)
  ofFiles <- find_ofFiles(ofDir)

  # -- copy the single files to the parent results dir
  if(file.exists(ofFiles$hogs)){
    tmp <- unlist(ofFiles[c("SpeciesIDs", "SequenceIDs", "ogs", "hogs", "speciesTree")])
  }else{
    tmp <- unlist(ofFiles[c("SpeciesIDs", "SequenceIDs", "ogs", "speciesTree")])
  }

  tmp <- tmp[!is.na(tmp)]
  tmp <- tmp[file.exists(tmp)]
  nu  <- file.copy(unlist(tmp), resDir)

  # -- copy the blast results to the parent dir
  tmp <- ofFiles$blast$blastFile
  tmp <- tmp[!is.na(tmp)]
  tmp <- tmp[file.exists(tmp)]
  nu  <- file.copy(tmp, resDir)

  # -- copy the ortholog results to the parent dir
  tmp <- ofFiles$blast$orthologFile
  tmp <- tmp[!is.na(tmp)]
  tmp <- tmp[file.exists(tmp)]
  nu  <- file.copy(tmp, resDir)
}

