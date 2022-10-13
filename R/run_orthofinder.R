#' @title Accessory functions for OrthoFinder
#' @description
#' \code{run_orthofinder} Methods to call OrthoFinder from R, check installation
#' parse results and copy files needed for GENESPACE
#' @name run_orthofinder
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param genomeIDs character vector of genomeIDs to consider
#' @param resultsDir file.path to the genespace /results directory
#' @param orthofinderDir file.path to the orthofinder run. When used with
#' copy_of2results, can move files from a separate orthofinder run to
#' the genespace directory without copying the thousands of orthofinder files.
#' @param overwrite deprecated, retained for backward compatibility. If an
#' orthofinder run exists, stops the run and returns a warning.
#' @param path2orthofinder file path or string the the $PATH to the OrthoFinder
#' program.
#' @param path2diamond file path or string the the $PATH to the DIAMOND2
#' program.
#' \cr
#' If called, \code{run_orthofinder} returns its own arguments.
#'
#'

#' @title Get the version of the orthofinder install
#' @description
#' \code{get_orthofinderVersion} Checks that orthofinder is installed and if so,
#' returns the installed version.
#' @rdname run_orthofinder
#' @export
get_orthofinderVersion <- function(path2orthofinder){
  # -- check if orthofinder is callable
  path <- path.expand(path2orthofinder)
  wh <- Sys.which(as.character(path))
  isThere <- basename(wh) == "orthofinder"
  # -- if orthofinder is installed, check version
  if(isThere){
    ver <- system2(path, "-h", stdout = TRUE)[2]
    ver <- strsplit(ver, " ")[[1]][3]
    vern <- strsplit(ver, ".", fixed = T)[[1]]
    vern <- as.numeric(sprintf("%s.%s%s", vern[1], vern[2], vern[3]))
    return(vern)
  }else{
    return(NA)
  }
}

#' @title Get the version of the DIAMOND install
#' @description
#' \code{get_diamondVersion} Checks that DIAMOND is installed and if so,
#' returns the installed version.
#' @rdname run_orthofinder
#' @export
get_diamondVersion <- function(path2diamond){
  path <- path.expand(path2diamond)
  chk <- tryCatch(
    {
      system2(path, "help",
              stdout = TRUE, stderr = FALSE)[1]
    },
    error = function(err) {
      return(NA)
    }
  )
  if(!is.na(chk)){
    ver <- strsplit(gsub(" |diamond|v", "", chk),"(", fixed = T)[[1]][1]
    vern <- strsplit(ver, ".", fixed = T)[[1]]
    vern <- as.numeric(sprintf("%s.%s%s", vern[1], vern[2], vern[3]))
    chk <- vern
  }
  return(chk)
}

#' @title Call OrthoFinder from R
#' @description
#' \code{run_orthofinder} Takes the results of init_genespace and calls
#' OrthoFinder.
#' @rdname run_orthofinder
#' @export
run_orthofinder <- function(gsParam,
                            genomeIDs = NULL,
                            overwrite = FALSE){

  ##############################################################################
  # 1. handle an existing run
  ofDir <- gsParam$paths$orthofinder
  resDir <- gsParam$paths$results
  ofRun <- length(list.files(ofDir)) > 0

  # -- 1.1 Handle overwrite = TRUE
  if(ofRun & overwrite)
    stop(strwrap(sprintf(
      "overwrite = TRUE is now deprecated ... to re-run orthofinder, you need to
      remove the directory %s and all of its contents, then re-call
      `run_orthofinder()`",
      ofDir)))
  if(ofRun){
    cat(strwrap(sprintf(
      "An existing orthofinder run exists in %s. To overwrite this result,
      remove this directory and all of its contents, then re-call
      `run_orthofinder()`. **NOT RE-RUNNING ORTHOFINDER** Instead just checking
      to make sure all files are OK ...",
      ofDir)))
  }

  # -- if no existing run, remove the orthofinder directory if it exists
  if(!ofRun){

    onewayBlast <- gsParam$params$onewayBlast
    diamondUltraSens <- gsParam$params$diamondUltraSens
    path2orthofinder <- gsParam$shellCalls$orthofinder
    nCores <- gsParam$params$nCores

    if(dir.exists(ofDir))
      unlink(ofDir, recursive = T)

    ############################################################################
    # 2. set up the directory structure and make sure things look good
    # -- make the tmp directory
    tmpDir <- gsParam$paths$tmp
    cat(strwrap(sprintf(
      "Copying files over to the temporary directory: %s",
      tmpDir), indent = 8, exdent = 16), sep = "\n")
    if(dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    dir.create(tmpDir)

    # -- copy peptides over to tmp directory
    # -- this allows for more genomeIDs in /peptide than just those in gsParam
    pepDir <- gsParam$paths$peptide
    genomeIDs <- gsParam$genomeIDs
    pepf <- file.path(pepDir, sprintf("%s.fa", genomeIDs))
    if(!all(file.exists(pepf)))
      stop(sprintf(
        "something is wrong with the peptide files. could not find: %s",
        paste(pepf[!file.exists(pepf)], collapse = "\n")))
    nu <- file.copy(pepf, tmpDir)

    ############################################################################
    # 3. Get the orthofinder command
    ofComm <- sprintf(
      "-f %s -t %s -a %s %s %s -X -o %s",
      tmpDir, nCores, nCores,
      ifelse(onewayBlast, "-1", ""),
      ifelse(diamondUltraSens, "-S diamond_ultra_sens", ""),
      ofDir)
    ofComm <- gsub("  ", " ", gsub("  ", " ", ofComm))
    if(!is.na(path2orthofinder)){
      cat(strwrap(sprintf(
        "Running the following command in the shell: `%s %s`.This can take a
        while. To check the progress, look in the `WorkingDirectory` in the
        output (-o) directory", path2orthofinder, ofComm),
        indent = 8, exdent = 16), sep = "\n")

      outp <- system2(
        path2orthofinder,
        ofComm,
        stdout = TRUE, stderr = TRUE)
      cat(paste(c("\t", outp), collapse = "\n\t"))
      gsParam$ofFiles <- find_ofFiles(orthofinderDir = ofDir)
    }else{
      cat(strwrap(
        "Could not find a valid path to the orthofinder program. To run
        orthofinder, ensure that the program is in the $PATH, the call the
        function from the shell using: \n", indent = 0, exdent = 8), sep = "\n")
      cat(sprintf("orthofinder %s", ofComm))
    }
  }else{
    gsParam$ofFiles <- find_ofFiles(orthofinderDir = ofDir)
  }
  return(gsParam)
}

#' @title find orthofinder results files
#' @description
#' \code{find_ofFiles} Find raw orthofinder results files. Typically for
#' internal use with copy_of2results. If called directly, use with caution.
#' @rdname run_orthofinder
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
  if(length(ofDir) != 1)
    stop("could not find a valid orthofinder run in ",
         orthofinderDir, "\n")

  # -- make a skeleton list
  ofPaths <- list(
    SpeciesIDs = NA, SequenceIDs = NA, ogs = NA, hogs = NA, blast = NA)

  ##############################################################################
  # 2. get paths to all expected files
  # -- paths to single files (ogs, hogs, speciesIDs, sequenceIDs)
  ofWd <- file.path(ofDir, "WorkingDirectory")

  # -- species IDs
  tmp <- file.path(ofWd, "SpeciesIDs.txt")
  if(file.exists(tmp)){
    ofPaths$SpeciesIDs <- tmp
    sids <- read_orthofinderSpeciesIDs(ofPaths$SpeciesIDs)
  }else{
    sids <- NULL
  }

  # -- sequence IDs
  tmp <- file.path(ofWd, "SequenceIDs.txt")
  if(file.exists(tmp))
    ofPaths$SequenceIDs <- tmp

  # -- OGs
  tmp <- file.path(ofDir, "Orthogroups", "Orthogroups.tsv")
  if(file.exists(tmp))
    ofPaths$ogs <- tmp

  # -- HOGs
  tmp <- file.path(ofDir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
  if(file.exists(tmp))
    ofPaths$hogs <- tmp

  if(!is.null(sids)){
    # -- blast files
    genome1 <- genome2 <- blastFile <- genNum1 <- genNum2 <-
      orthologFile <-NULL
    genomeIDs <- names(sids)
    blMd <- data.table(CJ(genome1 = genomeIDs, genome2 = genomeIDs))

    blMd[,`:=`(genNum1 = sids[genome1], genNum2 = sids[genome2])]
    blMd[,blastFile := file.path(ofWd, sprintf(
      "Blast%s_%s.txt.gz", genNum1, genNum2))]
    blMd$blastFile[!file.exists(blMd$blastFile)] <- NA

    # -- ortholog files
    blMd[,orthologFile := file.path(
      ofDir,
      "Orthologues",
      sprintf("Orthologues_%s", genome1),
      sprintf("%s__v__%s.tsv", genome1, genome2))]
    blMd$orthologFile[!file.exists(blMd$orthologFile)] <- NA
    ofPaths$blast <- blMd
  }

  return(ofPaths)
}

#' @title find orthofinder results files in the results directory
#' @description
#' \code{find_gsResults} Locate the required orthofinder results for genespace,
#' which are stored in the /results directory. Typically for internal use. If
#' called directly, use with caution.
#' @rdname run_orthofinder
#' @import data.table
#' @export
find_gsResults <- function(genomeIDs = NULL,
                           resultsDir = NULL,
                           gsParam = NULL){

  if(!is.null(gsParam)){
    genomeIDs <- gsParam$genomeIDs
    resultsDir <- gsParam$paths$results
  }

  if(is.null(genomeIDs))
    stop("genomeIDs or gsParam must be specified\n")
  if(is.null(resultsDir))
    stop("resultsDir or gsParam must be specified\n")

  ofPaths <- list(
    SpeciesIDs = NA, SequenceIDs = NA, ogs = NA, hogs = NA, blast = NA)

  ##############################################################################
  # 1. Check that all the four main files are there
  fs <- file.path(
    resultsDir,
    c("SpeciesIDs.txt", "SequenceIDs.txt", "Orthogroups.tsv", "N0.tsv"))
  names(fs) <- c("SpeciesIDs", "SequenceIDs", "ogs", "hogs")

  if(!all(file.exists(fs))){
    cat(strwrap(sprintf(
    "**NOTE** OrthoFinder files: SpeciesIDs.txt, SequenceIDs.txt,
    Orthogroups.tsv and N0.tsv must all be in the %s directory", resultsDir),
    indent = 0, exdent = 8), sep = "\n")
  }else{
    ############################################################################
    # 2. Make sure that genomeIDs match speciesIDs
    sids <- read_orthofinderSpeciesIDs(file.path(resultsDir, "SpeciesIDs.txt"))
    chk <- all(genomeIDs %in% names(sids)) && all(names(sids) %in% genomeIDs)
    if(!chk){
      cat(strwrap(sprintf(
        "**NOTE** There appears to be a problem with the orthofinder results in
        %s. GenomeIDs used in this GENESPACE run are `%s`. However, the species
        IDs in the orthofinder run are `%s`. These need to match exactly.",
        resultsDir,
      paste(genomeIDs, collapse = ", "),
      paste(names(sids), collapse = ", ")), indent = 0, exdent = 8), sep = "\n")
    }else{

      ##########################################################################
      # 3. Get paths to blast and ortholog files
      genNum1 <- genNum2 <- orthologFile <- genome1 <- genome2 <-
        blastFile <- NULL
      blMd <- data.table(CJ(genome1 = genomeIDs, genome2 = genomeIDs))
      blMd[,`:=`(genNum1 = sids[genome1], genNum2 = sids[genome2])]
      blMd[,blastFile := file.path(resultsDir, sprintf(
        "Blast%s_%s.txt.gz", genNum1, genNum2))]
      blMd$blastFile[!file.exists(blMd$blastFile)] <- NA

      blMd[,orthologFile := file.path(
        resultsDir,
        sprintf("%s__v__%s.tsv", genome1, genome2))]
      blMd$orthologFile[!file.exists(blMd$orthologFile)] <- NA

      ofPaths$SpeciesIDs <- file.path(resultsDir, "SpeciesIDs.txt")
      ofPaths$SequenceIDs <- file.path(resultsDir, "SequenceIDs.txt")
      ofPaths$ogs <- file.path(resultsDir, "Orthogroups.tsv")
      ofPaths$hogs <- file.path(resultsDir, "N0.tsv")
      ofPaths$blast <- blMd
    }
  }
  return(ofPaths)
}

#' @title copy raw orthofinder results to the genespace directory
#' @description
#' \code{copy_of2results} copy files from an orthofinder run to the genespace
#' results directory. Typically for internal use. If called directly, use with
#' caution.
#' @rdname run_orthofinder
#' @import data.table
#' @export
copy_of2results <- function(orthofinderDir, resultsDir, genomeIDs){
  ofDir <- check_filePathParam(orthofinderDir)
  resDir <- check_filePathParam(resultsDir)
  ofFiles <- find_ofFiles(ofDir)

  # -- copy the single files to the parent results dir
  tmp <- unlist(ofFiles[c("SpeciesIDs", "SequenceIDs", "ogs", "hogs")])
  tmp <- tmp[!is.na(tmp)]
  tmp <- tmp[file.exists(tmp)]
  nu <- file.copy(unlist(tmp), resDir)

  # -- copy the blast results to the parent dir
  tmp <- ofFiles$blast$blastFile
  tmp <- tmp[!is.na(tmp)]
  tmp <- tmp[file.exists(tmp)]
  nu <- file.copy(tmp, resDir)

  # -- copy the ortholog results to the parent dir
  tmp <- ofFiles$blast$orthologFile
  tmp <- tmp[!is.na(tmp)]
  tmp <- tmp[file.exists(tmp)]
  nu <- file.copy(tmp, resDir)

  gsFiles <- find_gsResults(resultsDir = resultsDir, genomeIDs = genomeIDs)
  return(gsFiles)
}
