#' @title The GENESPACE pipeline
#'
#' @description
#' \code{run_genespace} Run the entire GENESPACE pipeline, from begining to end,
#' with one function call.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param overwriteBed logical, should the bed file be re-created and
#' overwritten?
#' @param overwriteSynHits logial, should the annotated blast files be
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
                          overwrite = FALSE,
                          overwriteBed = overwrite,
                          overwriteSynHits = overwrite){
  gsParam$paths$rawOrthofinder <- gsParam$paths$orthofinder
  ##############################################################################
  # 1. Run orthofinder ...
  cat("\n############################", strwrap(
    "1. Running orthofinder (or parsing existing results)",
    indent = 0, exdent = 8), sep = "\n")

  ##############################################################################
  # -- 1.1 Check for existing parsed orthofinder results
  cat("\tChecking for existing orthofinder results ...\n")
  gsParam <- set_syntenyParams(gsParam)
  if("synteny" %in% names(gsParam)){
    noResults <- is.na(gsParam$synteny$SpeciesIDs)
  }else{
    noResults <- TRUE
  }

  ##############################################################################
  # -- 1.2 If no results exist, check for raw orthofinder run
  if(noResults){
    if(dir.exists(gsParam$paths$rawOrthofinder)){
      chkOf <- find_ofFiles(gsParam$paths$rawOrthofinder)
      noOrthofinder <- is.na(chkOf[[1]])
    }else{
      noOrthofinder <- TRUE
    }
  }else{
    noOrthofinder <- FALSE
  }

  ##############################################################################
  # -- 1.3 If raw results exist, copy them over
  if(!noOrthofinder && noResults){
    with(gsParam, copy_of2results(
      orthofinderDir = paths$rawOrthofinder, resultsDir = paths$results))

    if(dir.exists(gsParam$paths$rawOrthofinder)){
      chkOf <- find_ofFiles(gsParam$paths$rawOrthofinder)
      noOrthofinder <- is.na(chkOf[[1]])
      noResults <- noOrthofinder
    }else{
      noOrthofinder <- TRUE
    }
  }

  # noResults <- is.na(gsParam$synteny$SpeciesIDs)
  if(!noResults)
    cat("\t... found existing run, not re-running orthofinder\n")

  ##############################################################################
  # -- 1.4 if no orthofinder run, make one
  if(noResults)
    tmp <- run_orthofinder(gsParam = gsParam, verbose = TRUE)

  gsParam <- set_syntenyParams(gsParam)

  ##############################################################################
  # -- 1.5 get the files in order if the run is complete
  if(noResults){
    chkOf <- find_ofFiles(gsParam$paths$orthofinder)
    noOrthofinder <- is.na(chkOf[[1]])
    if(noOrthofinder)
      stop("could not find orthofinder files!")
    with(gsParam, copy_of2results(
      orthofinderDir = paths$orthofinder,
      resultsDir = paths$results))
  }
  gsParam <- set_syntenyParams(gsParam)

  ##############################################################################
  # -- 1.6 if the species tree exists, re-order the genomeIDs
  tmp <- gsParam$synteny$speciesTree

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

  # -- 1.7 if useHOGs, check if the N0.tsv file exists
  useHOGs <- gsParam$params$useHOGs
  if(useHOGs){
    if(is.na(gsParam$synteny$hogs)){
      useHOGs <- FALSE
    }else{
      if(!file.exists(gsParam$synteny$hogs)){
        useHOGs <- FALSE
      }
    }
  }
  gsParam$params$useHOGs <- useHOGs
  useHOGs <- NULL

  ##############################################################################
  # 2. Get the data ready for synteny
  hasBed <- FALSE
  bedf <- gsParam$synteny$combBed
  if(file.exists(bedf))
    hasBed <- is.data.table(read_combBed(bedf))
  if(overwriteBed)
    hasBed <- FALSE

  if(!hasBed){
    cat("\n############################", strwrap(
      "2. Combining and annotating bed files w/ OGs and tandem array info ... ",
      indent = 0, exdent = 8), sep = "\n")
    bed <- annotate_bed(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "2. Annotated/concatenated bed file exists", indent = 0, exdent = 8),
      sep = "\n")
  }

  ##############################################################################
  # 3. Annotate the blast files ...
  # -- First make sure that the blast files are all there, then go through
  # and annotate them with the combined bed file
  # -- This also makes the first round of dotplots
  # -- 3.1 check if all the synHits files exist. If so, and !overwriteSynHits
  # don't re-annotate
  hasHits <- FALSE
  if(all(file.exists(gsParam$synteny$blast$synHits)))
    if(!overwriteSynHits)
      hasHits <- TRUE

  # -- 3.2 iterate through and annotate all synHits files
  if(!hasHits){
    cat("\n############################", strwrap(
      "3. Combining and annotating the blast files with orthogroup info ...",
      indent = 0, exdent = 8), sep = "\n")
    gsParam <- annotate_blast(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "3. Annotated/blast files exists", indent = 0, exdent = 8),
      sep = "\n")
  }

  ##############################################################################
  # 4. Run synteny
  # -- goes through each pair of genomes and pulls syntenic anchors and the hits
  # nearby. This is the main engine of genespace

  cat("\n############################", strwrap(
    "4. Flagging synteny for each pair of genomes ...",
    indent = 0, exdent = 8), sep = "\n")
  gsParam <- synteny(gsParam = gsParam, overwrite = overwriteSynHits)
  gparsv <<- gsParam

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

