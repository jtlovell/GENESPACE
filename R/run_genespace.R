#' @title The GENESPACE pipeline
#'
#' @description
#' \code{run_genespace} Run the entire GENESPACE pipeline, from begining to end,
#' with one function call.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param overwriteBed logical, should the bed file be re-created and
#' overwritten?
#'
#' @details The full genespace pipeline is as follows.
#' \enumerate{
#' \item Details coming soon.
#' }
#'
#' @return a gsParam list.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#' @export
run_genespace <- function(gsParam,
                          overwriteBed = FALSE){

  ##############################################################################
  # ad hoc function to make sure the combined bed file is ok
  check_combBedFile <- function(bedFile, checkOGs = FALSE){
    bedPass <- FALSE
    if(file.exists(bedFile)){

      # -- check that the bed file has all the columns needed
      bedNames <- c("chr", "start", "end", "id", "ofID", "pepLen",
                    "ord", "genome", "arrayID", "isArrayRep", "globOG", "globHOG",
                    "synOG", "inblkOG", "noAnchor", "og")
      bedHed <- strsplit(readLines(bedFile, 1), "\t")[[1]]

      # -- check that all the genomes are present in the bed file
      if(all(bedNames %in% bedHed)){
        bed <- read_combBed(bedFile)
        bedPass <- all(!is.na(bed$og)) && all(gsParam$genomeIDs %in% bed$genome)
        if(checkOGs){
          bedPass <- bedPass & all(!is.na(bed$synOG))
          if(gsParam$params$orthofinderInBlk & bedPass)
            bedPass <- all(!is.na(bed$inblkOG))
        }
      }
    }
    return(bedPass)
  }

  ##############################################################################
  # 1. Run orthofinder ...
  # -- The function checks if an orthofinder run has been completed and if so,
  # (optionally) moves the files and returns a new gsParam object with the
  # updated paths
  cat("\n############################", strwrap(
    "1. Running orthofinder (or parsing existing results)", indent = 0, exdent = 8), sep = "\n")
  gsParam <- run_orthofinder(gsParam = gsParam, verbose = TRUE)

  # -- get the files in order if the run is complete
  gsParam <- run_orthofinder(gsParam = gsParam, verbose = FALSE)

  # -- if the species tree exists, re-order the genomeIDs
  tmp <- gsParam$ofFiles$speciesTree

  if(requireNamespace("ape", quietly = T)){
    if(!is.na(tmp) && !is.null(tmp)){
      if(file.exists(tmp) && length(gsParam$genomeIDs) > 2){
        treLabs <- ape::ladderize(ape::read.tree(tmp))$tip.label
        cat(strwrap(sprintf(
          "re-ordering genomeIDs by the species tree: %s",
          paste(treLabs, collapse = ", ")), indent = 8, exdent = 16),
          sep = "\n")
        gsParam$genomeIDs <- treLabs[treLabs %in% gsParam$genomeIDs]
      }
    }
  }

  ##############################################################################
  # 2. Annotate the bed file ...
  # -- This first checks if an annotated bed file (if overwriteBed = T) exists
  # and, if not, concatenates all the bed files and builds a new combBed.txt
  # file including array information
  hasBed <- check_combBedFile(
    bedFile = file.path(gsParam$paths$results, "combBed.txt"))
  if(hasBed && overwriteBed)
    hasBed <- FALSE

  if(!hasBed){
    cat("\n############################", strwrap(
      "2. Combining and annotating the bed files with orthogroups and tandem
      array information ... ", indent = 0, exdent = 8), sep = "\n")
    bed <- annotate_bed(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "2. Annotated/concatenated bed file exists", indent = 0, exdent = 8), sep = "\n")
  }

  ##############################################################################
  # 3. Annotate the blast files ...
  # -- First make sure that the blast files are all there, then go through
  # and annotate them with the combined bed file
  # -- This also makes the first round of dotplots
  gsf <- find_gsResults(
    resultsDir = gsParam$paths$results,
    genomeIDs = gsParam$genomeIDs, verbose = FALSE)
  gids <- unique(unlist(gsf$blast[,1:2]))
  if(!(all(!is.na(unlist(gsf[1:4]))) && all(gsParam$genomeIDs %in% gids)))
    stop("could not find complete blast hits")

  cat("\n############################", strwrap(
    "3. Combining and annotating the blast files with orthogroup info ...",
    indent = 0, exdent = 8), sep = "\n")
  gsParam <- annotate_blast(gsParam = gsParam)

  ##############################################################################
  # 4. Run synteny
  # -- goes through each pair of genomes and pulls syntenic anchors and the hits
  # nearby. This is the main engine of genespace
  cat("\n############################", strwrap(
    "4. Flagging synteny for each pair of genomes ...",
    indent = 0, exdent = 8), sep = "\n")
  gsParam <- synteny(gsParam = gsParam)

  ##############################################################################
  # 5. Build syntenic orthogroups
  # -- in the case of a bunch of haploid genomes, this just aggregates
  # orthogroups and splits them by membership in syntenic regions
  # -- in the case of polyploid genomes, this also runs orthofinder in blocks,
  # then re-runs synteny and re-aggregates blocks,.
  cat("\n############################", strwrap(
    "5. Building synteny-constrained orthogroups ...",
    indent = 0, exdent = 8), sep = "\n")
  gsParam <- build_synOGs(gsParam)

  ##############################################################################
  # 6. Integrate syntenic positions across genomes
  # -- the main purpose here is to find the interpolated syntenic positions
  # across all genomes
  # -- this also calculates syntenic block coordinates and phased reference-
  # specific coordinates to feed directly into plot_riparian
  cat("\n############################", strwrap(
    "6. Integrating syntenic positions across genomes ...",
    indent = 0, exdent = 8), sep = "\n")
  gsParam <- integrate_synteny(gsParam)

  ##############################################################################
  # 7. Build the pan-genome annotations and riparian plots
  # -- loops through each genome and uses that as the reference
  glab <- align_charLeft(sprintf("%s: ",c("genome", gsParam$genomeIDs)))
  names(glab) <- c("head", gsParam$genomeIDs)
  cat("\n############################", strwrap(
    "7. Building pan-genome annotations and riparian plots",
    indent = 0, exdent = 8),
    sprintf("\t...%sn pos. || n array || n NS ortho", glab["head"]), sep = "\n")
  for(refi in gsParam$genomeIDs){
    tmp <- plot_riparian(
      gsParam = gsParam, verbose = FALSE, refGenome = refi)
    ripSourceData <- tmp$sourceData
    save(ripSourceData, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_geneOrder_riparianSourceData.rda", refi)))
    ripPlotObj <- tmp$ggplotObj
    save(ripPlotObj, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_geneOrder_riparianGgplotObj.rda", refi)))
    tmp <- plot_riparian(
      gsParam = gsParam, verbose = FALSE, refGenome = refi, useOrder = FALSE)
    ripSourceData <- tmp$sourceData
    save(ripSourceData, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_bp_riparianSourceData.rda", refi)))
    ripPlotObj <- tmp$ggplotObj
    save(ripPlotObj, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_bp_riparianGgplotObj.rda", refi)))
    pg <- pangenome(
      gsParam = gsParam, refGenome = refi, verbose = F)
    pgout <- fread(file.path(gsParam$paths$pangenome,
                             sprintf("%s_refPangenomeAnnot.txt", refi)))
    with(pgout, cat(sprintf(
      "\t...%s%s || %s || %s\n",
      glab[refi], uniqueN(pgID), sum(!isArrayRep & !isNSOrtho), sum(isNSOrtho))))
  }

  gpFile <- file.path(gsParam$paths$results, "gsParams.rda")
  cat("\n############################", strwrap(sprintf(
    "GENESPACE run complete!\n All results are stored in %s in the following subdirectories:",
    gsParam$paths$wd), indent = 0, exdent = 0),
    "\traw dotplots           : /dotplots (...rawHits.pdf)",
    "\tsyntenic block dotplots: /dotplots (...synHits.pdf)",
    "\tannotated blast files  : /syntenicHits",
    "\tannotated/combined bed : /results/combBed.txt",
    "\tsyntenic block coords. : /results/blkCoords.txt",
    "\tsyn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt",
    "\tpan-genome annotations : /pangenome (...PangenomeAnnot.txt)",
    "\triparian plots         : /riparian",
    "\tgenespace param. list  : /results/gsParams.rda",
    "############################",
    strwrap(sprintf(
      "**NOTE** the genespace parameter object is returned or can be loaded
      into R via `load('%s', verbose = TRUE)`. Then you can customize your
      riparian plots by calling `plot_riparian(gsParam = gsParam, ...)`. The
      source data and ggplot2 objects are also stored in the /riparian
      directory and can also be accessed by `load(...)`. ",
      gpFile), indent = 0, exdent = 8),
    strwrap(
      "**NOTE** To query genespace results by position or gene,
      use `query_genespace(...)`. See specifications in ?query_genespace for
      details.",  indent = 0, exdent = 8),
    "############################",
    sep = "\n")
  save(gsParam, file = gpFile)
  return(gsParam)
}

