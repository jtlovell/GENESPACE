#' @title The GENESPACE pipeline
#'
#' @description
#' \code{run_genespace} Run the entire GENESPACE pipeline from beginning to end
#' with one function call.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param overwrite logical, should all raw files be overwritten except
#' orthofinder results
#' @param overwriteBed logical, should the bed file be re-created and
#' overwritten?
#' @param overwriteSynHits logial, should the annotated blast files be
#' overwritten?
#'
#' @details The function calls required to run the full genespace pipeline are
#' printed below. See each function for detailed descriptions. Also, see
#' `init_genespace`for details on parameter specifications.
#'
#' \enumerate{
#' \item `run_orthofinder` runs orthofinder or finds and copies over data from
#' a previous run.
#' \item `set_syntenyParams` converts parameters in the gsParam list into a
#' matrix of file paths and parameters for each pairwise combination of query
#' and target genomes
#' \item `annotate_bed` reads in all of the bed files, concatenates them and
#' adds some important additional information, including gene rank order,
#' orthofinder IDs, orthogroup information, tandem array identity etc.
#' \item `annotate_blast` reads in all the blast files and adds information from
#' the annotated/combined bed file
#' \item `synteny` is the main engine for genespace. this flags syntenic blocks
#' and make dotplots
#' \item `build_synOGs` integrates syntenic orthogroups across all blast files
#' \item `run_orthofinderInBlk` optionally re-runs orthofinder within each
#' syntenic block, returning phylogenetically hierarchical orthogroups (HOGs)
#' \item `integrate_synteny` interpolates syntenic position of all genes across
#' all genomes
#' \item `pangenome` combines positional and orthogroup information into a
#' single matrix anchored to the gene order coordinates of a single reference
#' \item `plot_riparian` is the primary genespace plotting routine, which stacks
#' the genomes and connects syntenic regions to color-coded reference
#' chromosomes
#' }
#'
#' @return a gsParam list.
#'
#' @examples
#' \dontrun{
#' ###############################################
#' # -- change paths to those valid on your system
#' genomeRepo <- "~/path/to/store/rawGenomes"
#' wd <- "~/path/to/genespace/workingDirectory"
#' path2mcscanx <- "~/path/to/MCScanX/"
#' ###############################################
#'
#' dir.create(genomeRepo)
#' dir.create(wd)
#' rawFiles <- download_exampleData(filepath = genomeRepo)
#'
#' parsedPaths <- parse_annotations(
#'   rawGenomeRepo = genomeRepo,
#'   genomeDirs = c("human", "chicken"),
#'   genomeIDs = c("human", "chicken"),
#'   presets = "ncbi",
#'   genespaceWd = wd)
#'
#' gpar <- init_genespace(
#'   wd = wd, nCores = 4,
#'   path2mcscanx = path2mcscanx)
#'
#' out <- run_genespace(gpar)
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
  gsParam <- synteny(gsParam = gsParam)

  cat("\n############################", strwrap(
    "Generating dotplots ... ",
    indent = 0, exdent = 8), sep = "\n")
  nu <- plot_hits(gsParam = gsParam)

  ##############################################################################
  # 5. Build syntenic orthogroups

  cat("\n############################", strwrap(
    "5. Building synteny-constrained orthogroups ... ",
    indent = 0, exdent = 8), sep = "\n")

  # -- in the case of a bunch of haploid genomes, this just aggregates
  # orthogroups and splits them by membership in syntenic regions
  cat("\tAggregating syntenic orthogroups ... ")
  gsParam <- build_synOGs(gsParam)
  cat("Done\n")

  # -- in the case of polyploid genomes, this also runs orthofinder in blocks,
  # then re-runs synteny and re-aggregates blocks,.
  if(gsParam$params$orthofinderInBlk){

    # -- returns the gsparam obj and overwrites the bed file with a new og col
    cat("\t##############\n\tRunning Orthofinder within syntenic regions\n")
    gsParam <- run_orthofinderInBlk(
      gsParam = gsParam, overwrite = overwriteSynHits, onlyInBuffer = TRUE)

    # -- takes the new og column and refreshes the sameOG column in blast files
    cat("\tDone!\n\tRe-annotating blast files ...\n")
    gsParam <- annotate_blast(gsParam = gsParam)

    # -- re-run synteny with new sameOG column
    cat("\tDone!\n\tRe-running synteny ...\n")
    gsParam <- synteny(gsParam = gsParam)
  }

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

  # -- for each reference genome
  for(refi in gsParam$genomeIDs){

    # -- gene-order riparian
    tmp <- plot_riparian(
      gsParam = gsParam, verbose = FALSE, refGenome = refi)
    ripSourceData <- tmp$sourceData
    save(ripSourceData, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_geneOrder_riparianSourceData.rda", refi)))
    ripPlotObj <- tmp$ggplotObj
    save(ripPlotObj, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_geneOrder_riparianGgplotObj.rda", refi)))

    # -- bp position riparian
    tmp <- plot_riparian(
      gsParam = gsParam, verbose = FALSE, refGenome = refi, useOrder = FALSE)
    ripSourceData <- tmp$sourceData
    save(ripSourceData, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_bp_riparianSourceData.rda", refi)))
    ripPlotObj <- tmp$ggplotObj
    save(ripPlotObj, file = file.path(gsParam$paths$riparian, sprintf(
      "%s_bp_riparianGgplotObj.rda", refi)))

    # -- pangenome
    pg <- pangenome(
      gsParam = gsParam, refGenome = refi, verbose = F)
    pgout <- read_pangenome(
      file.path(gsParam$paths$pangenome,
                sprintf("%s_refPangenomeAnnot.txt", refi)),
      which = "long")
    with(pgout, cat(sprintf(
      "\t...%s%s || %s || %s\n",
      glab[refi], uniqueN(pgID), sum(!isArrayRep & !isNSOrtho), sum(isNSOrtho))))
  }

  ##############################################################################
  # 8. Print summaries and return
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

