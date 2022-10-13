#' @title Find files and directories for a GENESPACE run
#'
#' @description
#' \code{init_genespace} Searches for desired genome files in the
#' raw genome repo director.
#'
#' @param genomeIDs character vector of length > 1, matching length
#' of speciesIDs, versions and ploidy. Specifies the name to assign
#' to each genome. This vector must be unique and can be any string
#' that begins with a letter (a-z, A-Z) and is alphanumeric. '.' and '_' are
#' allowed as long as they are not the first character.
#' @param ignoreTheseGenomes character string matching one of the genomeIDs that
#' will be used in the orthofinder -og run but not in the synteny search.
#' Suggested to ensure that there is an outgroup that predates any WGD that the
#' user would like to study.
#' @param orthofinderInBlk logical, should orthofinder be re-run within
#' syntenic regions? Highly recommended for polyploids. When called, HOGs within
#' blocks replace global HOGs or OGs. See useHOGs for more information.
#' @param ploidy integer string specifying ploidy of genome assemblies. This is
#' usually half of the actual ploidy, that is an inbred diploid usually is
#' represented by a haploid genome assembly.
#' @param path2orthofinder character string coercible to a file path that points
#' to the orthofinder executable. If orthofinder is in the path, specify with
#' "orthofinder"
#' @param path2diamond character string coercible to a file path that points
#' to the diamond executable. If diamond is in the path, specify with
#' "diamond"
#' @param path2mcscanx see path2orthofinder, except to the mcscanx directory.
#' This must contain the MCScanX_h folder.
#' @param onewayBlast logical of length 1, specifying whether one-way blasts
#' should be run via `orthofinder -1 ...`. This replaces orthofinderMethod =
#' "fast", but  uses `diamond2 --more-sensitive` whereas the previous method
#' used --fast specification. Substantial speed improvements in large runs with
#' little loss of fidelity.
#' @param diamondUltraSens logical of length 1, specifying whether the diamond
#' mode run within orthofinder should be --more-sensitive (default, FALSE) or
#' --ultra-sensitive.
#' @param useHOGs logical of length 1 or NA, specifying whether to use
#' phylogenetically hierarchical orthogroups (HOGs) or raw orthogroups. By
#' default (NA), this is decided internally by `annotate_bed`, where the
#' orthogroup type with members that best match the genome ploidy is used. In
#' general, HOGs should be used for any run where all genomes are haploid, since
#' they have been shown to have ~20% better precision than raw orthogroups.
#' However, in cases where we want both homeologs, HOGs may be problematic and
#' probably should not be used for syntenic region calculations. That said,
#' HOGs are always used for within-block orthofinder, which is also the default
#' when any genomes have ploidy > 1. So, the only way to use the deprecated
#' orthogroups.tsv for pan-genome calculation is to set useHOGs = FALSE AND
#' orthofinderInBlk = FALSE.
#' @param rawOrthofinderDir file.path of length 1, specifying the location of
#' an existing raw orthofinder run. Defaults to the $wd/orthofinder, but can
#' be any path point to a valid orthofinder run. If not a valid path, this is
#' ignored.
#' @param nGaps integer of length 1, specifying the -m param to mcscanx
#' for the primary MCScanX run. This acts on the results from the initial
#' MCScanX run.
#' @param blkSize integer of length 1, specifying the -s param to mcscanx
#' @param nSecondHits integer of length 1, specifying the number of blast
#' hits to include after masking.
#' @param synBuff Numeric > 0, specifying the distance from an anchor
#' to consider a hit syntenic. This parameter is also used to limit the search
#' radius in dbscan-based blk calculation. Larger values will return larger
#' tandem arrays but also may permit inclusion of spurious non-syntenic networks
#' @param nGapsSecond see nGaps, but passed to secondary hits after masking
#' primary hits.
#' @param blkSizeSecond see blkSize, but passed to the secondary scan if
#' nSecondaryHits > 0.
#' @param synBuffSecond see syntenyBuffer. Applied only to synteny
#' construction of secondary hits.
#' @param onlyOgAnchors logical, should only hits in orthogroups be considered
#' for anchors?
#' @param onlyOgAnchorsSecond logical should only hits in orthogroups be
#' considered for anchors in secondary blocks?
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#' @param wd file.path where the analysis will be run
#' @param maxOgPlaces xxxx
#' @param arrayJump xxxx
#' @param nSecondaryHits xxxx
#' @param minPepLen xxxx
#'
#' @param outgroup deprecated in V1. See ignoreTheseGenomes.
#' @param orthofinderMethod deprecated in V1. See onewayBlast.
#' @param speciesIDs deprecated in V1. See `parse_annotations`.
#' @param versionIDs deprecated in V1. See `parse_annotations`.
#' @param rawGenomeDir deprecated in V1. See `parse_annotations`.
#' @param gffString deprecated in V1. See `parse_annotations`.
#' @param pepString deprecated in V1. See `parse_annotations`.
#' @param diamondMode deprecated in V1. 'fast' mode is no longer available.
#' --ultra-sensitive is available via diamondUltraSens.
#' @param overwrite deprecated in V1. Results are never over-written.
#' @param verbose deprecated in V1. All updates are printed to the console
#'
#' @details Simple directory parser to find and check the paths to all
#' annotation and assembly files.
#'
#' @return A list containing paths to the raw files. If a file is not found,
#' path is returned as null and a warning is printed.
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @import data.table
#' @importFrom parallel detectCores
#' @export
init_genespace <- function(wd,

                           genomeIDs = NULL,
                           ploidy = 1,
                           ignoreTheseGenomes = NULL,

                           path2orthofinder = "orthofinder",
                           path2diamond = "diamond",
                           path2mcscanx = "MCScanX",

                           onewayBlast = FALSE,
                           orthofinderInBlk = any(ploidy > 1),
                           useHOGs = NA,
                           rawOrthofinderDir = NA,
                           diamondUltraSens = FALSE,

                           nCores = min(c(detectCores()/2, 16)),
                           maxOgPlaces = ploidy * 8,
                           blkSize = 5,
                           nGaps = 5,
                           synBuff = 100,
                           arrayJump = ceiling(synBuff/2),
                           onlyOgAnchors = TRUE,

                           nSecondaryHits = 0,
                           nGapsSecond = nGaps * 2,
                           blkSizeSecond = blkSize,
                           onlyOgAnchorsSecond = FALSE,

                           # -- deprecated arguments here for backwards compat.
                           outgroup = ignoreTheseGenomes,
                           nSecondHits = nSecondaryHits,
                           synBuffSecond = NULL,
                           orthofinderMethod = NULL,
                           speciesIDs = NULL,
                           minPepLen = NULL,
                           versionIDs = NULL,
                           rawGenomeDir = NULL,
                           diamondMode = NULL,
                           overwrite = NULL,
                           gffString = NULL,
                           pepString = NULL,
                           verbose = NULL){

  ##############################################################################
  ##############################################################################
  # -- internal functions w/o documentation
  ##############################################################################


  ##############################################################################
  check_outgroup <- function(outgroup, genomeIDs){
    x <- check_character(outgroup)
    x <- x[!duplicated(x)]
    x <- x[!is.null(x) & !is.na(x)]
    x <- x[x %in% genomeIDs]
    if(length(x) == 0)
      x <- NA
    if(sum(!genomeIDs %in% outgroup) < 2)
      stop(sprintf("%s genomes available after excluding outgroup. Needs to be > 1",
                   sum(!genomeIDs %in% outgroup)))
    return(x)
  }

  ##############################################################################
  check_ploidy <- function(ploidy, genomeIDs, outgroup){
    x <- check_integer(ploidy)
    if(any(is.na(x)) || any(is.null(x)))
      stop("problem with ploidy - could not coerce all to integers\n")
    if(length(x) == 1)
      x <- rep(x, length(genomeIDs))
    if(length(x) != length(genomeIDs))
      stop("problem with ploidy - length not same as genomeIDs (or 1)\n")
    names(x) <- genomeIDs
    if(!is.null(outgroup) || !is.na(outgroup))
      x <- x[!names(x) %in% outgroup]
    return(x)
  }

  ##############################################################################
  check_ogPlaces <- function(maxOgPlaces, genomeIDs, outgroup){
    x <- check_integer(maxOgPlaces)
    if(any(is.na(x)) || any(is.null(x)))
      stop("problem with maxOgPlaces - could not coerce all to integers\n")
    if(length(x) == 1)
      x <- rep(x, length(genomeIDs))
    if(length(x) != length(genomeIDs))
      stop("problem with maxOgPlaces - length not same as genomeIDs (or 1)\n")
    names(x) <- genomeIDs
    if(!is.null(outgroup) || !is.na(outgroup))
      x <- x[!names(x) %in% outgroup]
    return(x)
  }

  ##############################################################################
  check_genomeIDs <- function(genomeIDs){
    x <- check_character(genomeIDs)
    if(any(is.na(x)) || any(is.null(x)))
      stop("could not coerce to genomeIDs to character string\n")
    if(any(duplicated(x)))
      stop("genomeIDs are not all unique\n")
    if(length(x) < 2)
      stop("need 2 or more genomeIDs\n")

    okChars <- c(LETTERS, letters, 0:9, "_", ".")

    y <- lapply(x, function(y) strsplit(y, "")[[1]])
    allGood <- sapply(y, function(z) all(z %in% okChars))
    noDotFirst <- sapply(y, function(z) z[1] != ".")

    if(!all(allGood))
      stop("some genomeIDs include non-alphanumeric characters aside from _ and .\n")
    if(!all(noDotFirst))
      stop("some genomeIDs start with .\n")

    return(x)
  }

  ##############################################################################
  check_wd <- function(wd){
    wdChk <- path.expand(check_character(x = wd, onlySingleValue = T))
    if(is.na(wdChk)){
      stop(sprintf("could not coerce wd: `%s` to a file.path", wd))
    }else{
      if(!dir.exists(wdChk)){
        stop(sprintf("wd: `%s` does not exist", wd))
      }else{
        if(!dir.exists(file.path(wdChk, "peptide"))){
          stop(sprintf(
            "could not find `peptide` directory in wd: %s.
            This should contain the parsed peptide fasta files",
            wd))
        }else{
          if(!dir.exists(file.path(wdChk, "bed"))){
            stop(sprintf(
              "could not find `bed` directory in wd: %s.
            This should contain the parsed bed files",
              wd))
          }else{
            cat(strwrap(sprintf("PASS: `%s`", wd), exdent = 8), sep = "\n")
          }
        }
      }
    }
    return(wdChk)
  }

  ##############################################################################
  ##############################################################################

  ##############################################################################
  ##############################################################################
  # 0. Check if old parameters are specified and return notes if so
  ##############################################################################
  # -- 0.1 parameters that are now in parse_annotations only
  toParse <- c(
    speciesIDs, versionIDs, rawGenomeDir, gffString, pepString, minPepLen)
  if(any(!is.null(toParse))){
    depr <- paste(toParse[!is.null(toParse)], collapse = ", ")
    cat(strwrap(sprintf(
      "**NOTE** argument(s) '%s' passed to `init_genespace` are now deprecated.
      Raw gff and peptide fasta parsing must be completed before calling
      `init_genespace`. See `parse_annotations()` for details",
      depr), indent = 0, exdent = 8), sep = "\n")
  }

  ##############################################################################
  # -- 0.2 "fast" parameters(orthofinderMethod, diamondMode)
  if(!is.null(orthofinderMethod))
    cat(strwrap(
      "**NOTE** argument 'orthofinderMethod' is now deprecated. Previous 'fast'
      method is now approximated by `onewayOF = TRUE`, which uses
      `orthofinder -1`, where reciprocal blast hits are not calculated for
      intergenomic searches.",
      indent = 0, exdent = 8), sep = "\n")

  if(!is.null(diamondMode))
    cat(strwrap(
      "**NOTE** argument 'diamondMode' is now deprecated. Two modes of diamond
      searches are now available. The default is `--more-sensitive`. To set
      `--ultra-sensitive`, use `diamondUltraSens = TRUE`, which sets
      `orthofinder -S diamond_ultra_sens`. The previously provided `--fast`
      method is no longer available.",
      indent = 0, exdent = 8), sep = "\n")

  ##############################################################################
  # -- 0.3 other parameters that are no longer used.
  if(!is.null(overwrite))
    cat(strwrap(
      "**NOTE** argument 'overwrite' is now deprecated. Any case where an
      existing directory exists and files would previously be overwritten
      result in a warning and the function stops with instructions about
      manually removing the files in question. This is to avoid inadvertant
      file deletions.",
      indent = 0, exdent = 8), sep = "\n")

  if(!is.null(verbose))
    cat(strwrap(
      "**NOTE** argument 'verbose' is now deprecated. GENESPACE now always
      provides updates to the console when run interactively in R.",
      indent = 0, exdent = 8), sep = "\n")

  ##############################################################################
  ##############################################################################
  # 1. Parameter checking
  ##############################################################################
  # -- 1.1 genome labels etc.
  cat("Checking Working Directory ... ")
  wd <- check_wd(wd)

  # -- make sure that the genome IDs are good
  cat("Checking user-defined parameters ...\n\tGenome IDs & ploidy ... ")
  if(is.null(genomeIDs)){
    pepDir <- file.path(wd, "peptide")
    pepFiles <- list.files(pepDir, pattern = ".fa$")
    if(length(pepFiles) == 0 || !dir.exists(pepDir))
      stop(sprintf(
        "if genomeIDs are not given, there must be peptide annotations in %s",
        wd))
    genomeIDs <- gsub(".fa$", "", pepFiles)
    pepFiles <- pepDir <- NULL
  }
  genomeIDs <- check_genomeIDs(genomeIDs)

  # -- make sure the outgroup is unique and one of the genomeIDs
  if(!identical(outgroup, ignoreTheseGenomes) && !is.null(ignoreTheseGenomes)){
    cat(strwrap(sprintf(
      "**NOTE** parameter `outgroup` is now deprecated and replaced with
      `ignoreTheseGenomes`. Setting ignoreTheseGenomes to %s",
      outgroup), indent = 8, exdent = 16), sep = "\n")
    outgroup <- ignoreTheseGenomes
  }
  outgroup <- check_outgroup(outgroup = outgroup, genomeIDs = genomeIDs)

  # -- check that ploidy is good
  ploidy <- check_ploidy(
    ploidy = ploidy, genomeIDs = genomeIDs, outgroup = outgroup)
  maxOgPlaces <- check_ogPlaces(
    maxOgPlaces = maxOgPlaces, genomeIDs = genomeIDs, outgroup = outgroup)
  cat(sprintf(
    "\n\t\t%s\n\tOutgroup ... ",
    paste(apply(cbind(align_charLeft(genomeIDs), ploidy), 1, function(x)
      paste(x, collapse = ": ")), collapse = "\n\t\t")))
  if(is.na(outgroup)){
    cat("NONE\n")
  }else{
    cat(strwrap(sprintf("%s", paste(outgroup, collapse = ", ")),
                indent = 0, exdent = 16))
  }

  ##############################################################################
  # -- 1.2 synteny parameters
  # -- number of cores
  cat("\tn. parallel processes ... ")
  nCores <- check_integer(
    x = nCores, min = 1, max = Inf, default = min(c(detectCores()/2, 16)))

  # -- block size
  cat(sprintf("%s\n\tcollinear block size ... ", nCores))
  blkSize <- check_integer(x = blkSize, min = 1, max = Inf, default = 5)

  # -- n. gaps
  cat(sprintf("%s\n\tn gaps in collinear block ... ", blkSize))
  nGaps <- check_integer(x = nGaps, min = 1, max = Inf, default = 5)

  # -- synteny buffer
  cat(sprintf("%s\n\tsynteny buffer size... ", nGaps))
  synBuff <- check_integer(x = synBuff, min = 1, max = Inf, default = 100)

  # -- og anchors
  cat(sprintf("%s\n\tonly orthogroups hits as anchors ... ", synBuff))
  onlyOgAnchors <- check_logical(onlyOgAnchors)

  # -- n secondary hits
  cat(sprintf("%s\n\tn secondary hits ... ", onlyOgAnchors))
  nSecondaryHits <- check_integer(
    x = nSecondaryHits, min = 0, max = Inf, default = 0)

  # -- if nSecondary hits > 0, parameterize those, otherwise set to NA
  if(nSecondaryHits > 0){
    cat(sprintf("%s\n\tn gaps in secondary blocks ... ", nSecondaryHits))
    nGapsSecond <- check_integer(
      x = nGapsSecond, min = 1, max = Inf, default = 5)
    cat(sprintf("%s\n\tsecondardy block size ... ", nGapsSecond))
    blkSizeSecond <- check_integer(
      x = blkSizeSecond, min = 1, max = Inf, default = 5)
    cat(sprintf("%s\n\tonly orthogroup hits in secondary block ... ",
                blkSizeSecond))
    onlyOgAnchorsSecond <- check_logical(onlyOgAnchorsSecond)
    cat(sprintf("%s\n ... ", onlyOgAnchorsSecond))
  }else{
    cat(sprintf("%s\n", nSecondaryHits))
    nGapsSecond <- NA
    blkSizeSecond <- NA
    onlyOgAnchorsSecond <- NA
  }

  ##############################################################################
  # -- 1.3 basic argument checking (just do the checking, no reporting)
  # -- check that the arguments are specified correctly
  path2orthofinder <- check_character(
    x = path2orthofinder, default = "orthofinder", onlySingleValue = T, na.rm = T)
  path2mcscanx <- check_character(
    x = path2mcscanx, default = "MCScanX", onlySingleValue = T, na.rm = T)
  diamondUltraSens <- check_logical(diamondUltraSens, onlySingleValue = T)
  onewayBlast <- check_logical(onewayBlast, onlySingleValue = T)

  if(!is.null(diamondMode))
    cat(strwrap(
      "**NOTE** 'diamondMode' specification is deprecated in GENESPACE V1. If
      you wish to increase sensitivity, set `diamond_ultra_sens` to TRUE. Other
      diamond modes (e.g. --fast) are no longer supported in GENESPACE\n", indent = 8, exdent = 8),
      sep = "\n")

  orthofinderInBlk <- check_logical(orthofinderInBlk, onlySingleValue = T)
  if(!is.na(useHOGs))
    useHOGs <- check_logical(useHOGs, onlySingleValue = T)
  rawOrthofinderDir <- check_filePathParam(check_character(
    x = rawOrthofinderDir, onlySingleValue = T))

  params <- list(
    useHOGs = useHOGs, nCores = nCores,
    orthofinderInBlk = orthofinderInBlk, blkSize = blkSize, nGaps = nGaps,
    synBuff = synBuff, onlyOgAnchors = onlyOgAnchors,
    nSecondaryHits = nSecondaryHits, blkSizeSecond = blkSizeSecond,
    nGapsSecond = nGapsSecond, onlyOgAnchorsSecond = onlyOgAnchorsSecond)

  ##############################################################################
  ##############################################################################
  # 2. Check annotation files
  cat("Checking annotation files (.bed and peptide .fa):\n")
  inFiles <- check_annotFiles(path = wd, genomeIDs = genomeIDs)

  ##############################################################################
  ##############################################################################
  # 3 Check dependencies (Orthofinder / diamond2)
  ##############################################################################
  # -- 3.1 define the paths and versions
  cat("Checking dependencies ...\n")
  orthoInstall <- get_orthofinderVersion(path2orthofinder)
  diamondInstall <- get_diamondVersion(path2diamond)

  ##############################################################################
  # -- 3.2 if an old install return a warning
  if(orthoInstall < 2.54 && !is.na(orthoInstall)){
    cat(strwrap(sprintf(
      "**WARNING!!** OrthoFinder version %s is installed but GENESPACE needs
      >= v2.54. Setting path2orthofinder as NA. Install an up-to-date version or
      run OrthoFinder in an environment with a current install.",
      orthoInstall),
      indent = 6, exdent = 6),
      sep = "\n")
    orthoInstall <- NA
    path2orthofinder <- NA
  }

  ##############################################################################
  # -- 3.3 if not DIAMOND2 return a warning
  if(diamondInstall < 2 && !is.na(orthoInstall)){
    cat(strwrap(sprintf(
      "**WARNING!!** DIAMOND version %s is installed but OrthoFinder needs
      DIAMOND2. Setting path2diamond and path2orthofinder as NA. Install an
      up-to-date version or run OrthoFinder in an environment with DIAMOND2.",
      diamondInstall),
      indent = 6, exdent = 6),
      sep = "\n")
    orthoInstall <- NA
    diamondInstall <- NA
    path2orthofinder <- NA
    path2diamond <- NA
  }

  ##############################################################################
  # -- 3.4 if no orthofinder install and orthofinderInBlk, set to FALSE
  if(is.na(orthoInstall) && orthofinderInBlk){
    cat(strwrap(
      "**WARNING!!**, orthofinderInBlk was set to TRUE; however, this
    routine requires a valid path to the OrthoFinder program to be specified
    within R. Therefore, orthofinderInBlk has been set to FALSE. If this is not
    desired, please re-run build_params with a valid path2orthofinder.\n",
      indent = 6, exdent = 6),
      sep = "\n")
    orthofinderInBlk <- FALSE
  }

  ##############################################################################
  # -- 3.5 if all looks good, tell the user so.
  if(!is.na(orthoInstall))
    cat(sprintf("\tFound valid path to OrthoFinder v%s: `%s`\n",
                orthoInstall, path2orthofinder))
  if(!is.na(diamondInstall))
    cat(sprintf("\tFound valid path to DIAMOND2 v%s: `%s`\n",
                diamondInstall, path2diamond))
  orthoFinderCall <- ifelse(
    is.na(orthoInstall) || is.na(path2orthofinder), NA, path2orthofinder)
  diamondCall <- ifelse(
    is.na(diamondInstall) || is.na(path2diamond), NA, path2diamond)

  ##############################################################################
  # -- 3.6 check for MCScanX_h
  tmp <- check_MCScanXhInstall(path2mcscanx)

  if(is.na(tmp)){
    cat(strwrap(sprintf(
      "**WARNING!!** Can't find valid path to the MCScanX_h executable in `%s`.
      Only plotting and query GENESPACE functions will be functional. If you
      want to run the main `synteny` function, you will need to specify a path
      to a valid MCScanX installation (see README).",
      path2mcscanx), exdent = 8),
      sep = "\n")
    mcscanxhInstall <- NA
  }else{
    mcscanxhInstall <- file.path(path.expand(path2mcscanx), "MCScanX_h")
    cat(sprintf("\tFound valid MCScanX_h executable: `%s`\n",
                mcscanxhInstall))
  }
  MCScanX_hCall <- ifelse(
    is.na(mcscanxhInstall) || is.na(path2mcscanx), NA, mcscanxhInstall)

  ##############################################################################
  # -- 3.7 combine 3rd party calls into a list
  shellCalls <- list(
    orthofinder = orthoFinderCall,
    diamond = diamondCall,
    mcscanx_h = MCScanX_hCall)

  ##############################################################################
  ##############################################################################
  # 4. Build the parameter output
  ##############################################################################
  # -- 4.1 check orthofinder run parameters
  orthofinderInBlk <- check_logical(orthofinderInBlk, onlySingleValue = T)
  useHOGs <- check_logical(useHOGs, onlySingleValue = T)
  ofParams <- list(diamondUltraSens, onewayBlast, orthofinderInBlk, useHOGs)
  names(ofParams) <- c("diamondUltraSens", "onewayBlast", "ofInBlk", "useHOGs")

  ##############################################################################
  # -- 4.2 get the final paths
  paths <- list(
    wd = wd,
    peptide = file.path(wd, "peptide"),
    bed = file.path(wd, "bed"),
    results = file.path(wd, "results"),
    syntenicHits = file.path(wd, "syntenicHits"),
    orthofinder = file.path(wd, "orthofinder"),
    dotplots = file.path(wd, "dotplots"),
    riparian = file.path(wd, "riparian"),
    pangenome = file.path(wd, "pangenome"),
    tmp = file.path(wd, "tmp"))

  ##############################################################################
  # -- 4.3 make directories if they do not exist
  tmp <- sapply(paths[names(paths) != "orthofinder"], function(x){
    if(!dir.exists(x))
      dir.create(x)
  })

  ##############################################################################
  # -- 5.4 return the list
  p <- list(
    genomeIDs = genomeIDs,
    outgroup = outgroup,
    ploidy = ploidy,
    maxOgPlaces = maxOgPlaces,
    shellCalls = shellCalls,
    paths = paths,
    params = c(params, ofParams)[!duplicated(names(c(params, ofParams)))])
}


