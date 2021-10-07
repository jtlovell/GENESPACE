#' @title Build orthofinder database
#'
#' @description
#' \code{run_orthofinder} Simplified blast database construction for and orthogroup
#' construction in orthofinder.
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param overwrite logical, should results be overwritten?
#' @param onlyCheckRun logical, should nothing be done but see if there is a run
#'
#' @details Build pairwise blast database without running both sides of
#' pairwise blasts. The genome with more genes is set as the query, then to
#' complete the blast database, a mirrored blast file is written with
#' the query genome as the target. Self blast runs are done as usual.
#'
#' @return A data.table that is needed for many downstream analyses. This has
#' metadata for each pair of genomes with the following columns:
#' \enumerate{
#' \item genome1/2: the genome IDs for each pairwise run
#' \item uGenes1/2: the number of unique genes for each genome
#' \item query/target: the identity of query/target genomes in each blast run
#' \item run/mirrorBlast: logical whether the blast was run or mirrored
#' \item gn1/2: orthofinder genomeID numbers for each genome
#' \item db/fa/1/2/blFile: diamond database, fasta files and blast files.
#' }
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @note \code{run_orthofinder} is a generic name for the functions documented.
#' \cr
#' If called, \code{run_orthofinder} returns its own arguments.
#'

#' @title fast_ofDb
#' @description
#' \code{fast_ofDb} fast_ofDb
#' @rdname run_orthofinder
#' @import data.table
#' @export
run_orthofinder <- function(gsParam, overwrite = FALSE){
  if(is.null(gsParam$params$synteny))
    stop("must run set_syntenyParams first\n")
  beenRun <- find_orthofinderResults(gsParam, onlyCheckRun = T)
  if(beenRun & !overwrite){
    warning("orthofinder run exists & !overwrite, so not running")
    gsParam <- find_orthofinderResults(gsParam, onlyCheckRun = F)
  }else{
    if(is.na(gsParam$paths$orthofinderCall) || !gsParam$paths$orthofinderCall){
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
    gsParam$params$orthofinderCall <- com
  }
  return(gsParam)
}

#' @title drop_unusedPeptides
#' @description
#' \code{drop_unusedPeptides} drop_unusedPeptides
#' @rdname run_orthofinder
#' @import data.table
#' @export
drop_unusedPeptides <- function(gsParam){
  f <- list.files(path = dirname(gsParam$paths$peptide[1]), full.names = F)
  fi <- basename(gsParam$paths$peptide)
  if(any(!f %in% fi)){
    fo <- f[!f %in% fi]
    for(i in fo)
      file.remove(file.path(dirname(gsParam$paths$peptide[1]), i))
  }
}

#' @title prep_ofDbFromPeptide
#' @description
#' \code{prep_ofDbFromPeptide} prep_ofDbFromPeptide
#' @rdname run_orthofinder
#' @import data.table
#' @export
prep_ofDbFromPeptide <- function(gsParam){

  # clean up leftover peptides if necessary
  drop_unusedPeptides(gsParam)

  # check for and remove orthofinder directory if needed
  if(dir.exists(gsParam$paths$orthofinder)){
    unlink(gsParam$paths$orthofinder, recursive = T)
  }

  # convert to orthofinder
  com <- sprintf(
    "%s -f %s -t %s -op -o %s 1>/dev/null 2>&1",
    gsParam$paths$orthofinderCall,
    dirname(gsParam$paths$peptide[1]),
    gsParam$params$nCores,
    gsParam$paths$orthofinder)
  system(com)
  return(com)
}

#' @title fast_ofDb
#' @description
#' \code{fast_ofDb} fast_ofDb
#' @rdname run_orthofinder
#' @import data.table
#' @export
fast_ofDb <- function(gsParam){
  runBlast <- genome1 <- genome2 <- db2 <- fa1 <- blFile <- NULL
  ##############################################################################
  # Ad hoc internal functions
  ##############################################################################
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
  ##############################################################################
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
  ##############################################################################
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
    p[,com := sprintf(
      "%s blastp %s --quiet -e %s -p %s --compress 1 -d %s -q %s -o %s",
      "diamond", .1, gsParam$params$nCores, db2, fa1, blFile)]
    return(p)
  }

  ##############################################################################
  # 1. Remove existing orthofinder directory if it exists
  if(dir.exists(gsParam$paths$orthofinder)){
    unlink(gsParam$paths$orthofinder, recursive = T)
  }

  ##############################################################################
  # 2. convert to orthofinder
  com <- sprintf(
    "%s -f %s -t %s -a 1 -op -o %s 1>/dev/null 2>&1",
    gsParam$paths$orthofinderCall,
    dirname(gsParam$paths$peptide[1]),
    gsParam$params$nCores,
    gsParam$paths$orthofinder)
  system(com)

  ##############################################################################
  # 3. place orthofinder input files in paths$orthofinder
  reorg_ofInput(gsParam$paths$orthofinder)
  si <- read_orthofinderSpeciesIDs(gsParam$paths$orthofinder)

  ##############################################################################
  # 4. get blast parameters
  p <- add_blastInfo2syn(gsParam = gsParam, ofSpeciesIDs = si)

  ##############################################################################
  # 5. Run blasts
  pwp <- subset(p, runBlast)
  for(i in 1:nrow(pwp)){
    if(gsParam$params$verbose)
      with(pwp[i,], cat(sprintf("\t\tRunning %s/%s (%s vs. %s)\n",
                                i, nrow(pwp), genome1, genome2)))
    system(pwp$com[i])
  }

  ##############################################################################
  # 6. Invert blasts if necessary
  if(gsParam$params$verbose)
    cat("\t\tDone!\n\tInverting intergenomic files ... ")
  ip <- subset(p, !runBlast & genome1 != genome2)
  for(i in 1:nrow(ip)){
    invert_blast(fileIn = ip$invertFile[i], fileOut = ip$blFile[i])
  }

  ##############################################################################
  # 7. Run orthofinder
  if(gsParam$params$verbose)
    cat("Done!\n\tRunning full orthofinder on pre-computed blast:\n")
  com <- with(gsParam, sprintf(
    "%s -b %s -t %s -a 1 -X",
    paths$orthofinderCall, paths$orthofinder, params$nCores, params$nCores))

  system(com)
  return(com)
}

#' @title prep_ofDbFromPeptide
#' @description
#' \code{prep_ofDbFromPeptide} prep_ofDbFromPeptide
#' @rdname run_orthofinder
#' @import data.table
#' @importFrom R.utils gunzip
#' @export
default_ofDb <- function(gsParam){
  if(all(is.na(gsParam$params$synteny)))
    stop("must run set_syntenyParams first\n")

  if(gsParam$params$verbose)
    cat("\tCleaning out orthofinder directory and prepping run\n")
  ##############################################################################
  # 1. clean out peptide directory of unused fastas, if necessary
  drop_unusedPeptides(gsParam)

  ##############################################################################
  # 2. Remove existing orthofinder directory if it exists
  if(dir.exists(gsParam$paths$orthofinder)){
    unlink(gsParam$paths$orthofinder, recursive = T)
  }

  ##############################################################################
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

  ##############################################################################
  # 4. run it
  if(dontRun){
    cat("\tCould not find valid orthofinder executable in the path\n",
        "\tRun the following command outside of R (assuming orthofinder is in the path):",
        "\n################\n",com,"\n################\n", sep = "")
  }else{
    system(com)
  }

  return(com)
}

#' @title find_orthofinderResults
#' @description
#' \code{find_orthofinderResults} find_orthofinderResults
#' @rdname run_orthofinder
#' @import data.table
#' @importFrom R.utils gunzip
#' @export
find_orthofinderResults <- function(gsParam, onlyCheckRun = F){
  ogsFile <- order_filesByMtime(
    path = gsParam$paths$orthofinder,
    pattern = "Orthogroups.tsv",
    recursive = T)
  if(onlyCheckRun){
    return(length(ogsFile) > 0)
  }else{
    if(length(ogsFile) > 1)
      warning("Found multiple orthofinder runs, only using the most recent\n")
    if(length(ogsFile) == 0)
      stop("Can't find the 'orthogroups.tsv' file\n\tHave you run orthofinder yet?\n")

    ofResDir <- dirname(dirname(ogsFile[1]))

    pfile <- file.path(ofResDir, "Gene_Duplication_Events")
    if(!dir.exists(pfile)){
      paralogsDir <- NA
    }else{
      paralogsDir <- pfile
    }

    orthfile <- file.path(ofResDir, "Orthologues")
    if(!dir.exists(orthfile)){
      orthologuesDir <- NA
    }else{
      orthologuesDir <- orthfile
    }

    blsFile <- order_filesByMtime(
      path = gsParam$paths$orthofinder,
      pattern = "diamondDBSpecies0.dmnd",
      recursive = T)

    blastDir <- dirname(blsFile[1])

    gsParam$paths$blastDir <- blastDir
    gsParam$paths$orthogroupsDir <- dirname(ogsFile[1])
    gsParam$paths$paralogsDir <- paralogsDir
    gsParam$paths$orthologuesDir <- orthologuesDir

    return(gsParam)
  }
}
