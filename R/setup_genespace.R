#' @title Find files and directories for a GENESPACE run
#'
#' @description
#' \code{setup_genespace} Searches for desired genome files in the
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
#' @param speciesIDs file path character vector. This is the subdirectory in
#' rawGenomeDir containing the files for each genomeID.
#' @param versionIDs file path character vector. This is the subdirectory in
#' each speciesID with the genome version that contains the files for
#' each genomeID.
#' @param rawGenomeDir single file path, pointing to the directory that
#' contains all the genome annotation raw files.
#' @param genomeDir character string coercible to a file.path. Specifying the
#' directory to store the parsed genome information (gff, peptide. )
#' @param blastDir character string coercible to a file.path. Specifying the
#' directory to store the raw blast files
#' @param syntenyDir character string coercible to a file.path. Specifying the
#' directory to store the synteny-constrained blast files
#' @param resultsDir character string coercible to a file.path. Specifying the
#' directory to store the results
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
#' @param path2diamond see path2orthofinder, except to the diamond executable
#' @param path2mcscanx see path2orthofinder, except to the mcscanx directory.
#' This must contain the MCScanX_h folder.
#' @param verbose logical length 1, should updates be printed to the console?
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#'
#' @details Simple directory parser to find and check the paths to all
#' annotation and assembly files.
#'
#' @return A list containing paths to the raw files. If a file is not found,
#' path is returned as null and a warning is printed.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @import R.utils
#' @import data.table
#' @export
setup_genespace <- function(
  genomeIDs,
  orthofinderMethod = "fast",
  outgroup = NULL,
  speciesIDs,
  versionIDs,
  rawGenomeDir,
  ploidy,
  nCores = detectCores()/2,
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = "MCScanX",
  genomeDir = file.path(getwd(), "genome"),
  blastDir = file.path(getwd(), "blast"),
  syntenyDir = file.path(getwd(), "synteny"),
  resultsDir = file.path(getwd(), "results"),
  gffString = "gene.gff",
  pepString = "pep|prot",
  verbose = TRUE){

  ##############################################################################
  # -- 1. Argument checking
  ##############################################################################

  # Check logical arguments
  verbose <- check_logicalArg(verbose)
  orthofinderMethod <- match.arg(orthofinderMethod, choices = c("fast", "default"))

  if(missing(outgroup))
    outgroup <- NULL
  if(!is.null(outgroup)){
    outgroup <- outgroup[outgroup %in% genomeIDs]
  }else{
    outgroup <- NA
  }

  # check that the number of cores is within what the machine can handle
  nCores <- check_nCores(nCores)

  # check that paths to dependencies are valid
  if(!check_MCScanXhInstall(file.path(path2mcscanx, "MCScanX_h")))
    stop("Cant find MCScanX_h in", path2mcscanx, "check install\n")
  if (!check_orthofinderInstall(path2orthofinder))
    stop("cannot call orthofinder with", path2orthofinder,"\n")
  if (!check_diamondInstall(path2diamond))
    stop("cannot call diamond with", path2diamond,"\n")

  # check that the length of species, version and genome IDs match
  if(length(genomeIDs) != length(speciesIDs) |
     length(genomeIDs) != length(versionIDs))
    stop("genomeIDs, speciesIDs and versionIDs must all be the same length\n")

  # make sure ploidy is specified correctly.
  if(length(ploidy) == 1)
    ploidy <- rep(ploidy, length(genomeIDs))
  if(length(ploidy) != length(genomeIDs))
    stop("ploidy must be a numeric vector of the same length as genomeIDs\n")
  names(ploidy) <- genomeIDs

  # ensure that the species and other parameters are character vectors
  argInfo <- c(as.list(environment()), list())
  argClasses <- sapply(argInfo, class)
  argLens <- sapply(argInfo, length)
  charArgs <- c("genomeIDs", "speciesIDs", "versionIDs", "rawGenomeDir",
                "gffString","pepString")

  for(x in charArgs){
    if(argClasses[x] != "character")
      stop(x, "must be of class character\n")
  }

  # ensure that the rawGenomeDir exists
  if (!dir.exists(rawGenomeDir))
    stop(rawGenomeDir, "Does not exist\n")

  ##############################################################################
  # -- 2. Make output directories
  ##############################################################################

  # Genome directory and the peptide and gff subdirectories
  if(length(genomeDir) != 1)
    stop("genomeDir must be a single char. string coercible to a file.path\n")
  if(!dir.exists(genomeDir))
    dir.create(genomeDir, recursive = T)
  gffDir <- file.path(genomeDir, "gff")
  if(!dir.exists(gffDir))
    dir.create(gffDir)
  pepDir <- file.path(genomeDir, "peptide")
  if(!dir.exists(pepDir))
    dir.create(pepDir)

  # synteny, results and blast directories to store intermediate data
  if(length(blastDir) != 1)
    stop("blastDir must be a single char. string coercible to a file.path\n")
  if(!dir.exists(blastDir))
    dir.create(blastDir, recursive = T)
  if(length(syntenyDir) != 1)
    stop("syntenyDir must be a single char. string coercible to a file.path\n")
  if(!dir.exists(syntenyDir))
    dir.create(syntenyDir, recursive = T)
  if(length(resultsDir) != 1)
    stop("resultsDir must be a single char. string coercible to a file.path\n")
  if(!dir.exists(resultsDir))
    dir.create(resultsDir, recursive = T)
  if(verbose)
    cat("Created the following output directories:\n",
        "\tParsed gff:",gffDir,"\n",
        "\tParsed peptide fasta:",pepDir,"\n",
        "\tRaw BLASTs:",blastDir,"\n",
        "\tSyntenic BLASTs:",syntenyDir,"\n",
        "\tResults:",resultsDir,"\n")


  ##############################################################################
  # -- 3. Check that rawGenome files exist
  ##############################################################################

  # make paths to the expected locations of the raw files
  genomePaths <- file.path(rawGenomeDir, speciesIDs, versionIDs)
  names(genomePaths) <- genomeIDs

  # ensure that the search strings for the peptides and gffs look right
  strList <- list(pep = pepString, gff = gffString)
  if (length(unique(sapply(strList, length)))!=1)
    stop("strings to find files must all be the same length\n")
  if (!(all(sapply(strList, length) == 1) |
        all(sapply(strList, length) == length(genomeIDs))))
    stop("strings to find files must have a length=1 or the # of genomes\n")

  # rep the search string for each genome
  if(length(pepString) == 1)
    strList <- lapply(strList, function(x){
      x <- rep(x, length(genomeIDs))
      names(x) <- genomeIDs
      return(x)
    })

  # for each genomeID ...
  pathsOut <- list()
  for(i in genomeIDs){
    pathsOut[[i]] <- list()

    # get the expected path to the peptide and gff files
    annPath <- file.path(genomePaths[i], "annotation")

    pepPath <- list.files(
      path = annPath,
      pattern = strList$pep[i],
      full.names = T)
    gffPath <- list.files(
      path = annPath,
      pattern = strList$gff[i],
      full.names = T)

    # check that each has exactly one entry / genome
    if(length(pepPath) == 0)
      stop("Could not find annotation with ", strList$pep[i],
           " string in peptide file name for ",i, " in:", annPath,"\n")
    if(length(pepPath) > 1)
      stop("Search string, (",strList$pep[i],
           ") for peptide annotation is not unique. Found:", pepPath,"\n")

    if(length(gffPath) == 0)
      stop("Could not find annotation with ", strList$gff[i],
           " string in gff file name for ",i, " in:", annPath,"\n")
    if(length(gffPath) > 1)
      stop("Search string, (",strList$gff[i],
           ") for gff annotation is not unique. Found:", gffPath,"\n")

    # if all looks good, return the path to the raw annotation files
    pathsOut[[i]]$pepPath <- pepPath
    pathsOut[[i]]$gffPath <- gffPath
  }


  ##############################################################################
  # -- 4. Format parameters and return as a list
  ##############################################################################

  outList <- list(
    nCores = nCores,
    orthofinderMethod =  orthofinderMethod,
    peptideRaw = sapply(pathsOut, function(x) x$pepPath),
    gffRaw = sapply(pathsOut, function(x) x$gffPath),
    gff = gffDir,
    peptide = pepDir,
    blast = blastDir,
    synteny = syntenyDir,
    results = resultsDir,
    genomeIDs = genomeIDs,
    ploidy = ploidy,
    path2orthofinder = path2orthofinder,
    path2diamond = path2diamond,
    path2mcscanx = path2mcscanx,
    verbose = verbose,
    outgroup = outgroup)

  if(verbose)
    cat("\nFound the following gff3-formatted annotation files:\n\t")
  if(verbose)
    cat(outList$gffRaw, sep = "\n\t")

  if(verbose)
    cat("\nFound the following peptide fasta files:\n\t")
  if(verbose)
    cat(outList$peptideRaw, sep = "\n\t")

  return(outList)
}
