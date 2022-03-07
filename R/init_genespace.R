#' @title Find files and directories for a GENESPACE run
#'
#' @description
#' \code{init_genespace} Searches for desired genome files in the
#' raw genome repo director.
#'
#' @param genomeIDs character vector of length > 1, matching length
#' of speciesIDs, versions and ploidy. Specifies the name to assign
#' to each genome. This vector must be unique and can be any string
#' that begins with a letter (a-z, A-Z). An X will be appended to
#' IDs starting with special characters or numbers.
#' @param outgroup character string matching one of the genomeIDs that will
#' be used in the orthofinder -og run but not in the synteny search. Suggested
#' to ensure that there is an outgroup that predates any WGD that the user
#' would like to study.
#' @param orthofinderMethod character string either 'fast' or 'default'. See
#' build_OFdb for details.
#' @param orthofinderInBlk logical, should orthofinder be re-run within
#' syntenic regions? Highly recommended for polyploids
#' @param speciesIDs file path character vector. This is the subdirectory in
#' rawGenomeDir containing the files for each genomeID.
#' @param versionIDs file path character vector. This is the subdirectory in
#' each speciesID with the genome version that contains the files for
#' each genomeID.
#' @param rawGenomeDir single file path, pointing to the directory that
#' contains all the genome annotation raw files.
#' @param gffString character string of length 1, coercable to a
#' regular expression. File names in the 'annotation' subdirectories are
#' screened for this string. If a single, unique file is found, it is returned
#' as the gff source file. Otherwise an error is produced.
#' @param pepString character string of length 1, coercable to a
#' regular expression. See gffSearchString.
#' @param ploidy integer string specifying ploidy of genome assemblies. This is
#' usually half of the actual ploidy, that is an inbred diploid usually is
#' represented by a haploid genome assembly.
#' @param path2orthofinder character string coercible to a file path that points
#' to the orthofinder executable. If orthofinder is in the path, specify with
#' "orthofinder"
#' @param path2mcscanx see path2orthofinder, except to the mcscanx directory.
#' This must contain the MCScanX_h folder.
#' @param verbose logical length 1, should updates be printed to the console?
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#' @param wd file.path where the analysis will be run
#' @param minPepLen integer, the shortest peptide to be analyzed
#' @param overwrite logical, should existing directories be overwritten?
#' @param diamondMode character string to match one of the diamond search modes
#' @param path2diamond character string coercible to a file path that points
#' to the diamond executable. If diamond is in the path, specify with
#' "diamond"
#' @details Simple directory parser to find and check the paths to all
#' annotation and assembly files.
#'
#' @return A list containing paths to the raw files. If a file is not found,
#' path is returned as null and a warning is printed.
#'
#' @examples
#' \dontrun{
#' runwd <- file.path(getwd(), "testGenespace")
#' make_exampleDataDir(writeDir = runwd)
#'
#' #############################################################################
#' ###### Default method, can run orthofinder outside of R on a cluster etc. ###
#' gpar <- init_genespace(
#'   genomeIDs = c("human","chimp","rhesus"),
#'   speciesIDs = c("human","chimp","rhesus"),
#'   versionIDs = c("human","chimp","rhesus"),
#'   ploidy = rep(1,3),
#'   wd = runwd,
#'   nCores = 4,
#'   gffString = "gff",
#'   pepString = "pep",
#'   path2orthofinder = "orthofinder",
#'   path2mcscanx = "~/MCScanX",
#'   rawGenomeDir = file.path(runwd, "rawGenomes"))
#'
#' #############################################################################
#' ###### Slower but more accurate inference  ##################################
#' ###### especially for distantly related genomes or ploidy > 1 ###############
#' # add in orthofinder in block orthogroup ...
#' gpar <- init_genespace(
#'   genomeIDs = c("human","chimp","rhesus"),
#'   speciesIDs = c("human","chimp","rhesus"),
#'   versionIDs = c("human","chimp","rhesus"),
#'   ploidy = rep(1,3),
#'   wd = runwd,
#'   nCores = 4,
#'   orthofinderInBlk = TRUE,
#'   gffString = "gff",
#'   pepString = "pep",
#'   path2orthofinder = "orthofinder",
#'   path2mcscanx = "~/MCScanX",
#'   rawGenomeDir = file.path(runwd, "rawGenomes"))
#'
#' #############################################################################
#' ###### Fast visualization of closely related genome #########################
#' # just a quick visualization ... not recommended for a full run, although
#' # probably not much less sensitive with closely related genomes
#' # **NOTE** this method can only be done when orthofinder is installed in the
#' # path for this R session (e.g. opening R from a conda env with orthofinder)
#' # **NOTE** increasing minPepLen results in a loss of fidelity, but an
#' # increase in speed.
#'
#' # This is the parameterization used in the help files.
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
#' }
#' @import R.utils
#' @import data.table
#' @export
init_genespace <- function(genomeIDs,
                           wd,
                           outgroup = NULL,
                           speciesIDs,
                           versionIDs,
                           rawGenomeDir,
                           ploidy,
                           diamondMode = "more-sensitive",
                           orthofinderMethod = "default",
                           orthofinderInBlk = any(ploidy > 1),
                           nCores = detectCores()/2,
                           minPepLen = 20,
                           overwrite = T,
                           path2orthofinder = "orthofinder",
                           path2mcscanx = "MCScanX",
                           path2diamond = "diamond",
                           gffString = "gff",
                           pepString = "pep|prot",
                           verbose = TRUE){

  ##############################################################################
  # -- ad hoc function to build the skeleton parameters
  skeleton_params <- function(){
    list(genomes = list(
      genomeIDs = NULL, outgroup = NULL, ploidy = NULL),
      params = list(
        wd = NULL, orthofinderMethod = NULL, nCores = NULL,
        verbose = NULL, synteny = NULL),
      paths = list(
        rawGff = NULL, rawPeptide = NULL, gff = NULL, peptide = NULL,
        orthofinder = NULL, results = NULL, orthofinderCall = NULL,
        mcscanxCall = NULL))
  }

  ##############################################################################
  # -- ad hoc function to make sure the genomeIDs conform
  check_genomeIDs <- function(x){
    lets <- c(LETTERS, letters)
    x <- as.character(x)
    firstChar <- substr(x, 1, 1)
    if(!all(firstChar %in% lets))
      stop("all genomeIDs must start with a letter (a-z A-Z) character\n")
    if(any(grepl("[^a-zA-Z0-9]", x)))
      stop("non-alphanumeric characters are not allowed in genomeIDs\n")
    if(any(is.na(x)) || any(is.null(x)))
      stop("problem with genomeIDs - could not coerce to character string\n")
    if(any(duplicated(x)))
      stop("problem with genomeIDs - not all are unique\n")
    if(length(x) < 2)
      stop("problem with genomeIDs - need 2 or more genomeIDs\n")
    return(x)
  }

  ##############################################################################
  # -- ad hoc function to make sure the outgroup conforms
  check_outgroup <- function(x, genomeIDs){
    if(any(is.null(x))){
      return(NA)
    }else{
      x <- as.character(x)
      lets <- c(LETTERS, letters)
      x <- as.character(x)
      firstChar <- substr(x, 1, 1)
      if(!all(firstChar %in% lets))
        stop("outgroup ID must start with a letter (a-z A-Z) character\n")
      if(grepl("[^a-zA-Z0-9]", x))
        stop("non-alphanumeric characters are not allowed in outgroup ID\n")
      if(any(is.na(x)) || any(is.null(x)))
        stop("problem with outgroup ID - could not coerce to character string\n")
      if(any(duplicated(x)))
        stop("problem with outgroup ID - not all are unique\n")
      if(!all(x %in% genomeIDs))
        stop("problem with outgroup - needs to be listed in genomeIDs\n")
      return(x)
    }
  }

  ##############################################################################
  # -- ad hoc function to check if ploidy specs are ok
  check_ploidy <- function(x, genomeIDs){
    x <- as.integer(x)
    if(any(is.na(x)) || any(is.null(x)))
      stop("problem with ploidy - could not coerce all to integers\n")
    if(length(x) == 1)
      x <- rep(x, length(genomeIDs))
    if(length(x) != length(genomeIDs))
      stop("problem with ploidy - length not same as genomeIDs (or 1)\n")
    names(x) <- genomeIDs
    return(x)
  }

  ##############################################################################
  # -- ad hoc function to make sure the wd is ok and created.
  check_wd <- function(x){
    if(!dir.exists(dirname(x)))
      stop(sprintf(
        "parent directory of wd does not exist\n\tChoose a different wd or create %s first",
        dirname(x)))
    if(!dir.exists(x))
      dir.create(x, recursive = T)
    return(x)
  }

  ##############################################################################
  # -- ad hoc function to make sure raw files are present
  check_rawFiles <- function(path,
                             genomeIDs,
                             speciesIDs,
                             versionIDs,
                             pattern){
    if(!dir.exists(path))
      stop("can't find the directory for raw genomes:", path)

    if(length(genomeIDs) != length(speciesIDs) ||
       length(genomeIDs) != length(versionIDs))
      stop("genomeIDs, speciesIDs, and versionIDs must all be the same length\n")

    pattern <- as.character(pattern)
    if(is.na(pattern) || is.null(pattern) || length(pattern) != 1)
      stop("search pattern must be a single character string")

    pv <- file.path(path, speciesIDs, versionIDs, "annotation")
    if(any(!dir.exists(pv)))
      stop(sprintf(
        "problem with the raw genome repo directory structure:
        Cant find the following directories where the raw genome info should be stored: \n%s\n",
        paste(pv[!dir.exists(pv)], collapse = "\n")))

    names(pv) <- genomeIDs
    ps <- sapply(pv, USE.NAMES = T, simplify = F, function(x)
      list.files(
        path = x,
        pattern = pattern,
        full.names = T))
    vs <- sapply(ps, length)
    if(min(vs) == 0)
      stop(sprintf("problems with the raw genome files:
                     \tCould not find %s in: %s\n",
                   pattern, paste(pv[which(vs == 0)], collapse = "; ")))

    if(max(vs) > 1)
      stop(sprintf("problems with the raw genome files:
                     \tPattern (%s) not unique for: %s\n",
                   pattern, paste(pv[which(vs > 1)], collapse = "; ")))
    ps0 <- unlist(ps)
    names(ps0) <- names(ps)
    return(ps0)
  }

  ##############################################################################
  # -- ad hoc function to check min peptide length parameter
  check_minPepLen <- function(x){
    x <- as.integer(x)[1]
    if(is.na(x) || is.null(x))
      stop("minPepLen must be an integer\n")
    if(x < 0)
      x <- 0
    return(x)
  }

  ##############################################################################
  # -- ad hoc function to make the parsed annotation file paths
  make_parsedAnnotPaths <- function(path, genomeIDs){
    gp <- file.path(path, "gff")
    if(!dir.exists(gp))
      dir.create(gp)

    pp <- file.path(path, "peptide")
    if(!dir.exists(pp))
      dir.create(pp)

    gpo <- file.path(gp, sprintf("%s.gff.gz", genomeIDs))
    ppo <- file.path(pp, sprintf("%s.fa", genomeIDs))
    names(gpo) <- names(ppo) <- genomeIDs
    return(list(gff = gpo, pep = ppo))
  }

  ##############################################################################
  # -- ad hoc function to make genespace directories
  make_genespaceDirs <- function(path, overwrite = FALSE){
    of <- file.path(path, "orthofinder")
    re <- file.path(path, "results")
    if(!dir.exists(of) & overwrite)
      unlink(of, recursive = T)
    if(!dir.exists(re))
      dir.create(re)
    return(c(of, re))
  }

  ##############################################################################
  # -- ad hoc function to check the mcscanx install
  check_MCScanXhInstall <- function(path){
    path <- as.character(path)
    if(is.na(path) || is.null(path)){
      pth <- NA
    }else{
      path <- path.expand(path)
      pth <- file.path(path, "MCScanX_h")
      if(!dir.exists(path))
        stop("path to MCScanX is not valid, check that its specified correctly\n")
      if(!file.exists(pth))
        stop("found MCScanX folder, but not executable. Has it been installed with make?\n")

      chk <- tryCatch(
        suppressWarnings(system2(pth, "-h", stdout = TRUE, stderr = FALSE)),
        error = function(err) NA)
      if(!grepl("prefix_fn", chk[1]) || is.na(chk))
        pth <- NA
    }
    if(is.na(pth))
      stop("cannot find MCScanX_h executable in", pth,"\n")
    return(pth)
  }

  ##############################################################################
  # -- ad hoc function to check diamond installation
  check_diamondInstall <- function(path){
    path <- as.character(path)
    if(is.na(path) || is.null(path)){
      path <- NA
    }else{
      path <- path.expand(path)
      wh <- Sys.which(as.character(path))
      ex <- file.exists(path)
      if(!ex && wh == ""){
        path <- NA
      }else{
        chk <- tryCatch(
          suppressWarnings(system2(path, "help", stdout = TRUE, stderr = TRUE)),
          error = function(err) NA)
        if(!grepl("diamond v2", chk[1]))
          path <- NA
      }
    }
    return(path)
  }

  ##############################################################################
  # -- ad hoc function to check orthofinder install
  check_orthofinderInstall <- function(path, verbose = FALSE){
    if(is.na(path) || is.null(path)){
      path <- NA
    }else{
      path <- path.expand(path)
      wh <- Sys.which(as.character(path))
      ex <- file.exists(path)
      if(!ex && wh == ""){
        path <- NA
      }else{
        ver <- tryCatch(
          suppressWarnings(system2(path, "-h", stdout = TRUE, stderr = TRUE)),
          error = function(err) NA)
        if(!grepl("OrthoFinder", ver[2]) || is.na(ver) || length(ver) == 1){
          path <- NA
        }else{
          ver <- strsplit(ver[2], " ")[[1]][3]
          vern <- strsplit(ver, ".", fixed = T)[[1]]
          vern <- as.numeric(sprintf("%s.%s%s", vern[1], vern[2], vern[3]))
          if(vern < 2.52)
              warning(sprintf("Orthofinder >= 2.5.2 must be installed (path is to version %s)\n\tAssuming orthofinder will be run outside of R with v2.5.2 or later\n", ver))
        }
      }
    }
    return(path)
  }

  ##############################################################################
  # -- ad hoc function to check the nCores parameter
  check_nCores <- function(nCores){
    if(length(nCores > 1)) nCores <- nCores[1]
    nCores <- as.integer(nCores)
    if(is.null(nCores)) nCores <- detectCores()/2
    if(is.na(nCores)) nCores <- detectCores()/2
    if(nCores < 1) nCores <- 1
    if(nCores > detectCores()){
      warning(sprintf(
        "user specified %s cores, but only %s available\n\tSetting nCores to %s",
        nCores, detectCores(), detectCores()))
      nCores <- detectCores()
    }
    return(nCores)
  }

  ##############################################################################
  # -- ad hoc function to check which orthofinder method can be used
  choose_ofMethod <- function(path2orthofinder, orthofinderMethod = "default"){
    ofm <- match.arg(orthofinderMethod, choices = c("fast", "default"))
    tmp <- check_orthofinderInstall(path2orthofinder)
    isInstall <- !is.na(tmp)
    if(ofm == "fast" && isInstall){
      return("fast")
    }else{
      if(isInstall){
        return("default inside R")
      }else{
        return("default")
      }
    }
  }
  ##############################################################################
  # -- ad hoc function to check if orthofinder in block can be run
  choose_ofInBlkMethod <- function(path2orthofinder, orthofinderInBlk = FALSE){
    orthofinderInBlk <- check_logicalArg(orthofinderInBlk)
    return(orthofinderInBlk && check_orthofinderInstall(path2orthofinder))
  }

  ##############################################################################
  # -- ad hoc function to check diamond mode parameter
  choose_diamondMode <- function(path2diamond, diamondMode){
    path2diamond <- check_diamondInstall(path2diamond)
    arg <- match.arg(diamondMode, choices = c(
      "fast", "mid-sensitive", "sensitive", "more-sensitive",
      "very-sensitive", "ultra-sensitive"))

    return(sprintf("--%s", arg))
  }




  ##############################################################################
  # 1. make the skeleton and do basic checks
  ##############################################################################
  # make sure genomeIDs look good
  if(verbose)
    cat("Initializing GENESPACE run\n\tchecking genomeIDs ... ")
  p <- skeleton_params()
  p$genomes$genomeIDs <- check_genomeIDs(genomeIDs)
  if(verbose & all(outgroup %in% genomeIDs))
    cat(sprintf("PASS (%s)\n\tchecking outgroup ... ",
                paste(p$genomes$genomeIDs, collapse = ", ")))
  p$genomes$outgroup <- check_outgroup(outgroup, genomeIDs = genomeIDs)
  if(is.na(p$genomes$outgroup) & verbose)
    cat("PASS (no outgroup)\n")
  if(!is.na(p$genomes$outgroup) & verbose)
    cat(sprintf("PASS (%s)\n",
                paste(p$genomes$outgroup, collapse = ", ")))
  if(verbose & all(outgroup %in% genomeIDs))
    cat("\tchecking ploidy ... ")
  p$genomes$ploidy <- check_ploidy(ploidy, genomeIDs = genomeIDs)
  if(verbose & all(outgroup %in% genomeIDs))
    cat(sprintf("PASS (%s)\n\tchecking the number of parallel processes ... ",
                paste(p$genomes$ploidy, collapse = ", ")))

  p$params$nCores <- check_nCores(nCores)
  if(verbose)
    cat(sprintf("PASS (%s)\n\tVerbosity ... ", p$params$nCores))
  p$params$verbose <- check_logicalArg(verbose)
  if(verbose)
    cat(sprintf("PASS (%s)\n\tminPepLen ... ", p$params$verbose))
  p$params$minPepLen <- check_minPepLen(minPepLen)
  if(verbose)
    cat(sprintf("PASS (%s)\n\tchecking working directory ... ", p$params$minPepLen))

  # -- check the working directory
  p$params$wd <- check_wd(wd)
  setwd(p$params$wd)
  if(verbose)
    cat(sprintf("PASS (%s)\n", p$params$wd))

  # -- check for parsed gff annotation files
  if(verbose)
    cat("\tchecking parsed gff files ... ")
  tmp <- make_parsedAnnotPaths(
    path = p$params$wd,
    genomeIDs = p$genomes$genomeIDs)
  p$paths$gff <- tmp$gff
  foundGff <- FALSE
  if(!all(file.exists(p$paths$gff))){
    if(verbose)
      cat("not found - need to run parse_annotations (etc.)\n")
  }else{
    foundGff <- TRUE
    if(verbose)
      cat(sprintf("PASS (%s)\n", dirname(p$paths$gff)[1]))
  }

  # -- check for parsed peptide annotation files
  if(verbose)
    cat("\tchecking parsed peptide files ... ")
  p$paths$peptide <- tmp$pep
  foundPep <- FALSE
  if(!all(file.exists(p$paths$peptide))){
    if(verbose)
      cat("not found - need to run parse_annotations (etc.)\n")
  }else{
    foundPep <- TRUE
    if(verbose)
      cat(sprintf("PASS (%s)\n", dirname(p$paths$peptide)[1]))
  }

  # -- check for raw gff and peptide annotation files
  if(!foundPep || !foundGff){
    if(verbose)
      cat("Since some parsed gff and peptide files are missing, checking raw files\n")
    p$paths$rawGff <- check_rawFiles(
      path = rawGenomeDir,
      genomeIDs = p$genomes$genomeIDs,
      speciesIDs = speciesIDs,
      versionIDs = versionIDs,
      pattern = gffString)

    p$paths$rawPeptide <- check_rawFiles(
      path = rawGenomeDir,
      genomeIDs = p$genomes$genomeIDs,
      speciesIDs = speciesIDs,
      versionIDs = versionIDs,
      pattern = pepString)

    cat(sprintf("\tfound raw gff files:\n\t\t%s\n",
        paste(basename(p$paths$rawGff), collapse = "\n\t\t")))
    cat(sprintf("\tfound raw peptide files:\n\t\t%s\n",
        paste(basename(p$paths$rawPeptide),collapse = "\n\t\t")))
  }

  ##############################################################################
  # 2. check dependencies
  ##############################################################################
  # -- mcscanx
  if(verbose)
    cat("Checking dependencies and 3rd party installations\n\tMCScanX installation ... ")
  p$paths$mcscanxCall <- check_MCScanXhInstall(path2mcscanx)
  if(verbose)
    cat(sprintf("PASS (%s)\n\tOrthofinder installation ... ",
                p$paths$mcscanxCall))

  # -- orthofinder
  p$paths$orthofinderCall <- check_orthofinderInstall(
    path2orthofinder, verbose = verbose)
  if(is.na(p$paths$orthofinderCall)){
    if(verbose)
      cat(sprintf("no installation found at %s\n",
                  path2orthofinder))
  }else{
    if(verbose)
      cat(sprintf("PASS (%s)\n",
                  p$paths$orthofinderCall))
  }

  # -- diamond
  if(verbose)
    cat("\tOrthoFinder method ... ")
  p$params$orthofinderMethod <- choose_ofMethod(
    orthofinderMethod = orthofinderMethod,
    path2orthofinder = path2orthofinder)
  if(verbose){
    if(p$params$orthofinderMethod == "default" & orthofinderMethod == "fast" & verbose){
        cat(sprintf("switched to %s\n", p$params$orthofinderMethod))
      message(
        "\n\t######***NOTE***########\n",
        "\tCould not find orthofinder in system path, but fast method specified\n",
        "\tFor this run, setting orthofinder method = default (slower but more accurate)\n",
        "\tIf this isn't right, open R from a terminal with orthofinder in the path\n")
    }else{
      cat(sprintf("PASS - (%s)\n", p$params$orthofinderMethod))
    }
  }
  if(is.na(p$paths$orthofinderCall) & verbose)
    message("\n\t######***NOTE***########\n",
        "\tYou will need to conduct an orthofinder run outside of R prior to using GENESPACE\n",
        "\tFollow instructions printed from run_orthofinder\n")

  # -- orthofinder in block method checks
  if(verbose)
    cat("\tOrthofinder in block method ... ")
  p$params$orthofinderInBlk <- choose_ofInBlkMethod(
    path2orthofinder = path2orthofinder,
    orthofinderInBlk = orthofinderInBlk)
  if(verbose){
    if(!p$params$orthofinderInBlk & orthofinderInBlk){
      cat(
        "Could not find orthofinder in system path, but orthofinderInBlk specified\n",
        "For this run, setting orthofinderInBlk to FALSE (may not be appropriate for polyploids)\n",
        "If this isn't right, open R from a terminal with orthofinder in the path\n")
    }else{
      cat(sprintf("PASS (%s)\n", p$params$orthofinderInBlk))
    }
  }

  # -- diamond calls for just fast method
  if(p$params$orthofinderMethod == "fast"){
    if(verbose)
      cat("\tdiamond mode ... ")
    p$paths$diamondCall <- check_diamondInstall(path2diamond)
    p$params$diamondMode <- choose_diamondMode(
      path2diamond = path2diamond,
      diamondMode = diamondMode)
    if(verbose)
      cat(sprintf("PASS (%s)\n", p$params$diamondMode))
  }

  # -- final initialization
  if(verbose)
    cat("GENESPACE run successfully initialized\n")
  tmp <- make_genespaceDirs(p$params$wd, overwrite = overwrite)
  p$paths$orthofinder <- tmp[1]
  p$paths$results <- tmp[2]
  if(verbose)
    cat(sprintf("\tresults will be stored in %s\n", p$paths$results))

  p$paths$blastDir <- NA
  p$paths$orthogroupsDir <- NA
  p$paths$orthologuesDir <- NA

  return(p)
}
