#' @title Parsing of pairwise hits into synteny
#' @description
#' \code{synteny} The main GENESPACE engine to constrain genome-wide homology
#' hits to synteny.
#'
#' @name synteny
#'
#' @param gsParam a list containing all parameters for a GENESPACE run. See
#' init_genespace
#' @param genomeIDs an optional vector of genomeIDs to consider. If not
#' specified (default) taken from gsParam$genomeIDs$genomeIDs
#' @param hits data.table of hits
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
#' @param arrayBuffer Numeric length 1, specifying the maximum gene-rank order
#' position of two genes to be considered members of the same collinear array.
#' Set to synBuff + 1 if smaller than the maximum synBuff or sqrt(2*(synBuff^2))
#' if not specified.
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
#' @param selfRegionMask integer specifying the size of the region (radius, gene
#' rank order) surrounding self hits to mask for secondary/homeologous hits
#' @param dropInterleavesSmallerThan integer specifying the smallest block to
#' keep in the split interleaves synteny step
#' @param minRbhScore integer specifying the minimum blast bit score to allow
#' for a RBH to be included
#' @param path2mcscanx character string file.path pointing to the install
#' directory for the MCScanX program. In particular, there must be an
#' executable in this directory for MCScanX_h.
#' @param overwrite logical, should existing directories be overwritten?
#' @param ... additional arguments passed to set_syntenyParam().
#'
#' @details The main engine for GENESPACE synteny searching. This
#' finds syntenic 'anchors' that are high-confidence synteny- and homology-
#' constrained hits, then pulls nearby hits within a specified buffer around
#' 'anchor' hits between two genomes. Combined, this provides a framework to
#' both analyze syntenic duplicates (e.g. tandem arrays) and have high
#' confidence that only the desired types of hits (orthologs, homoeologs, etc.)
#' are considered.
#'
#' The basic premise is that we can find synteny in haploid genome comparisons (most
#' diploid species have haploid genome representations) best by first
#' constraining the initial search to the single best scoring pairwise hits for
#' each gene. Then, if desired further subsetting this set to only gene pairs
#' that fall into the same orthogroups. This approach effectively removes
#' relic syntenic regions from ancient whole genome duplications and produces
#' a set of high-confidence synteny "anchors" which serve as known positions of
#' syntenic (ortho/para/homeolo)gous regions. We then search in a fixed-radius
#' for nearest neighbors within a gene-rank order buffer around the anchors. For
#' intra-genomic hits, the self hits are the anchors and a buffer is calculated
#' via euclidean distance. For intra-genomic hits in polyploids, the self-hit
#' regions are masked and the synteny search is re-run to more effectively find
#' homeologs.
#'
#' It is important to note that this does NOT produce finalized block
#' coordinates, but only large-scale regions that are syntenic. These results
#' are usually passed to an additional orthofinder run, either globally, or for
#' polyploids or searches with paralogs, within-block pairwise runs. See
#' rerun_orthofinder.
#'
#' Parameterization of this function is complex and varies by the type of
#' contrast desired. To simplify parameterization, we have build a convenience
#' function that infers various parameters based on basic
#' genespace parameters and ploidy. Set set_syntenyParams for more information
#' about the input synParam data.table.
#'
#' @return updated gsParam list
#'
#' \cr
#' If called, \code{synteny} returns its own arguments.
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
#' gpar <- run_orthofinder(gsParam = gpar, overwrite = F)
#'
#' gpar <- synteny(gsParam = gpar) # use defaults
#'
#' # -- run again (need to set overwrite = T) with blkSize = 10
#' gpar$params$synteny <- NULL
#' gpar <- synteny(gsParam = gpar, blkSize = 10, overwrite = T)
#'
#' # -- run again with custom specs
#' # **NOTE** if params$synteny is a data.table, synteny will respect it. So,
#' # this can be generated and modified before running. If it is NULL (as above)
#' #synteny will calculate parameters internally with additional arguments.
#'
#' # here, increase the synteny buffer and
#' # make the blocks need to be bigger between human and chimp
#' gpar <- set_syntenyParams(gpar, synBuff = 200)
#' wh <- with(gpar$params$synteny, which(
#'   genome1 %in% c("humnan", "chimp") & genome2 %in% c("humnan", "chimp")))
#' gpar$params$synteny$blkSize[wh] <- 10
#' gpar <- synteny(gsParam = gpar, overwrite = T)
#' }
#'
#'
#' @rdname synteny
#' @import R.utils
#' @import data.table
#' @export
synteny <- function(gsParam,
                    genomeIDs = NULL,
                    overwrite = F, ...){

  ##############################################################################
  # 1. Set up the environment and parameter checking
  # -- genomeIDs
  setDTthreads(1)
  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(!any(genomeIDs %in% gsParam$genomes$genomeIDs))
      stop("specified genomeIDs dont look right\n")
  genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]

  # -- shortcuts and output files
  verbose <- gsParam$params$verbose
  writeTo <- gsParam$paths$results
  nCores <- gsParam$params$nCores
  gffFile <- file.path(writeTo, "gffWithOgs.txt.gz")
  blksFile <- file.path(writeTo, "syntenicBlocks.txt.gz")

  # -- set the synteny parameters
  if(is.data.table(gsParam$params$synteny))
    if(!all(genomeIDs %in% gsParam$params$synteny$genome1))
      gsParam$params$synteny <- NULL
  if(!is.data.table(gsParam$params$synteny)){
    tmp <- set_syntenyParams(gsParam, ...)
    gsParam <- set_syntenyParams(gsParam)
    if(identical(tmp$params$synteny, gsParam$params$synteny)){
      cat("Synteny Parameters have not been set! Setting to defaults\n")
    }else{
      cat("Synteny Parameters have not been set! Using user-defined settings\n")
      gsParam$params$synteny <- tmp$params$synteny
    }
  }
  synBuff <- max(gsParam$params$synteny$synBuff)

  # -- find the orthogroups
  if(is.na(gsParam$paths$orthogroupsDir)){
    if(verbose)
      cat("\tIndexing location of orthofinder results ... ")
    gsParam <- find_orthofinderResults(gsParam)
    if(verbose)
      cat("Done!\n")
  }

  ##############################################################################
  # 2. load and parse the gff
  # -- add in global syntenic orthogroup arrays, orthofinder IDs, etc
  if(verbose)
    cat("Parsing the gff files ... \n\tReading the gffs and adding orthofinder IDs ... ")
  gff <- annotate_gff(
    gsParam = gsParam,
    genomeIDs = genomeIDs)

  ##############################################################################
  # -- add arrays. This takes the specifed column in the gff and parses array
  # and array representatives. We do this again below if inblk OG is used.
  # With verbose = T, prints the counts of reps and total members / genome.
  if(verbose)
    cat("Done!\nDefining collinear orthogroup arrays ... \n")
  gff <- add_arrayReps2gff(
    gff = gff,
    synBuff = synBuff,
    ogColumn = "globOG",
    verbose = verbose)


  ##############################################################################
  # 3. Run the initial synteny builder
  # -- This is the full pipeline: split up synParams into chunks each with
  # nCores entries, read in blast results, parse to collinear hits, cull to
  # syntenic regions, form large regions, build syntenic blocks within regions.
  gsParam <- pipe_synteny(
    gsParam = gsParam,
    gff = gff,
    overwrite = T)

  ##############################################################################
  # 4. Build syntenic orthogroups from syntenic hits
  # -- This will also run orthofinder within the syntenic regions if desired and
  # combine syntenic and inblk orthogroups.
  gff <- pull_synOGs(
    gsParam = gsParam)

  ##############################################################################
  # 5. If method == "orthofinderInBlk", re-calculate synArrays and synteny
  recall_synteny <- function()



  fwrite(blks, file = blksFile, sep = "\t", quote = F, showProgress = F)
  if(verbose)
    cat(sprintf(
      "\tSynteny constraints - Done!\n\tSyntenic block coordinates written to /results/%s\n",
        basename(blksFile)))

  # -- get syntenic orthologs
  gff <- pull_synOGs(gsParam = gsParam)
  gffo <- combine_inblkSynOG(
    genomeIDs = genomeIDs,
    gff = gff,
    gsParam = gsParam)
  fwrite(gffo, file = gffFile, sep = "\t", quote = F, showProgress = F)
  if(verbose)
    cat("\tWrote gff to file: /results/gffWithOgs.txt.gz\n\tDone!\n")
  return(gsParam)
}

#' @title Set synteny parameters
#' @description
#' \code{set_syntenyParams} Calculate paramters for pairwise synteny search.
#' Other synteny functions require this as input.
#' @rdname synteny
#' @importFrom Biostrings readAAStringSet
#' @import data.table
#' @export
set_syntenyParams <- function(gsParam,
                              onlyOgAnchors = TRUE,
                              onlyOgAnchorsSecond = TRUE,
                              blkSize = 5,
                              blkSizeSecond = blkSize,
                              nGaps = 5,
                              nGapsSecond = nGaps,
                              nSecondHits = 0,
                              synBuff = 100,
                              synBuffSecond = synBuff,
                              selfRegionMask = synBuff * 2,
                              dropInterleavesSmallerThan = 2,
                              minRbhScore = 50,
                              arrayBuffer = sqrt(2*(synBuff^2))){

  ##############################################################################
  # 1. setup environment
  ##############################################################################

  genome1 <- genome2 <- tiebreak <- nGenes1 <- nGenes2 <- query <- target <- NULL
  nSecondHits1 <- nSecondHits2 <- nhits1 <- nhits2 <- NULL

  # -- check the parameters
  onlyOgAnchors <- as.logical(onlyOgAnchors[1])
  if(is.null(onlyOgAnchors) || is.na(onlyOgAnchors) || length(onlyOgAnchors) == 0)
    stop("cannot coerce onlyOgAnchors to logical\n")
  onlyOgAnchorsSecond <- as.logical(onlyOgAnchorsSecond[1])
  if(is.null(onlyOgAnchorsSecond) || is.na(onlyOgAnchorsSecond) || length(onlyOgAnchorsSecond) == 0)
    stop("cannot coerce onlyOgAnchors to logical\n")

  blkSize <- as.integer(blkSize)[1]
  if(is.null(blkSize) || is.na(blkSize) || length(blkSize) == 0)
    stop("cannot coerce blkSize to integer\n")
  blkSizeSecond <- as.integer(blkSizeSecond)[1]
  if(is.null(blkSizeSecond) || length(blkSizeSecond) == 0)
    stop("cannot coerce blkSizeSecond to integer\n")

  nGaps <- as.integer(nGaps)[1]
  if(is.null(nGaps) || is.na(nGaps) || length(nGaps) == 0)
    stop("cannot coerce nGaps to integer\n")
  nGapsSecond <- as.integer(nGapsSecond)[1]
  if(is.null(nGapsSecond) || length(nGapsSecond) == 0)
    stop("cannot coerce nGapsSecond to integer\n")

  synBuff <- as.integer(synBuff)[1]
  if(is.null(synBuff) || is.na(synBuff) || length(synBuff) == 0)
    stop("cannot coerce synBuff to integer\n")
  synBuffSecond <- as.integer(synBuffSecond)[1]
  if(is.null(synBuffSecond) || is.na(synBuffSecond) || length(synBuffSecond) == 0)
    stop("cannot coerce synBuffSecond to integer\n")

  selfRegionMask <- as.integer(selfRegionMask)[1]
  if(is.null(selfRegionMask) || is.na(selfRegionMask) || length(selfRegionMask) == 0)
    stop("cannot coerce selfRegionMask to integer\n")
  nSecondHits <- as.integer(nSecondHits)[1]
  if(is.null(nSecondHits) || is.na(nSecondHits) || length(nSecondHits) == 0)
    stop("cannot coerce nSecondHits to integer\n")
  dropInterleavesSmallerThan <- as.integer(dropInterleavesSmallerThan)[1]
  if(is.null(dropInterleavesSmallerThan) || is.na(dropInterleavesSmallerThan) || length(dropInterleavesSmallerThan) == 0)
    stop("cannot coerce dropInterleavesSmallerThan to integer\n")
  minRbhScore <- as.integer(minRbhScore)[1]
  if(is.null(minRbhScore) || is.na(minRbhScore) || length(minRbhScore) == 0)
    stop("cannot coerce minRbhScore to integer\n")

  # -- shorten some names
  genomeIDs <- gsParam$genomes$genomeIDs
  ploidy <- gsParam$genomes$ploidy[genomeIDs]
  blastDir <- gsParam$paths$orthofinder
  pepFiles <- gsParam$paths$peptide

  # -- Check the peptide files
  if(!all(file.exists(pepFiles))){
    wh <- which(!file.exists(pepFiles))
    stop(sprintf("the following genome annotations have not been correctly parsed:\n\t%s\n\tRerun parse_annotations\n",
                 paste(names(pepFiles)[wh], collapse = "\n\t")))
  }

  # -- Calculate the number of genes for each genome
  nGenes <- sapply(pepFiles, function(x) length(readAAStringSet(x)))
  if(min(nGenes) == 0){
    wh <- which(nGenes == 0)
    stop(sprintf("the following genome annotations have not been correctly parsed:\n\t%s\n\tRerun parse_annotations\n",
                 paste(names(pepFiles)[wh], collapse = "\n\t")))
  }

  ##############################################################################
  # 2. build the pairwise parameter data.table scaffold
  ##############################################################################

  # make the database of unique combinations of genomes
  cmb <- data.table(
    genome1 = names(ploidy),
    genome2 = names(ploidy))
  cmb <- cmb[,CJ(genome1, genome2)]

  # -- add ploidy, orthofinder IDs, n genes
  speciesIDs <- (1:length(genomeIDs))-1
  names(speciesIDs) <- genomeIDs[order(genomeIDs)]
  cmb[,`:=`(u = paste(genome1, genome2),
            ploidy1 = ploidy[genome1],
            ploidy2 = ploidy[genome2],
            ofID1 = speciesIDs[genome1],
            ofID2 = speciesIDs[genome2],
            nGenes1 = nGenes[genome1],
            nGenes2 = nGenes[genome2],
            nhits1 = ploidy[genome2],
            nhits2 = ploidy[genome1])]

  # -- Choose query and target genomes, based on the total number of genes
  cmb[,tiebreak := as.numeric(factor(genome1, levels = genomeIDs)) -
        as.numeric(factor(genome2, levels = genomeIDs))]
  cmb[,`:=`(query = ifelse(nGenes1 > nGenes2, genome1,
                           ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                  genome1, genome2)),
            target = ifelse(nGenes1 > nGenes2, genome2,
                            ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                   genome2, genome1))),
      by = "u"]
  cmb[,`:=`(runBlast = genome1 == query,
            mirrorBlast = genome1 == target & genome1 != query,
            tiebreak = NULL)]

  # -- choose the number of secondary hits
  cmb[,nSecondHits1 :=
        ifelse(nhits1 == 1 & nhits2 == 1,
               nSecondHits,
               ifelse(genome1 == genome2 & nhits2 > 1 & nSecondHits > 0,
                      (nSecondHits * nhits2/2) + (nhits2/2),
                      ifelse(nhits2 > 1 & nSecondHits > 0,
                             nSecondHits * nhits2,
                             nSecondHits)))]
  cmb[,nSecondHits2 :=
        ifelse(nhits1 == 1 & nhits2 == 1,
               nSecondHits,
               ifelse(genome1 == genome2 & nhits1 > 1 & nSecondHits > 0,
                      (nSecondHits * nhits1/2) + (nhits1/2),
                      ifelse(nhits1 > 1 & nSecondHits > 0,
                             nSecondHits * nhits1,
                             nSecondHits)))]

  ##############################################################################
  # 3. Add all the parameters into the output data.table.
  ##############################################################################
  params <- c("onlyOgAnchors","onlyOgAnchorsSecond","blkSize", "blkSizeSecond",
              "nGaps", "nGapsSecond",
              "synBuff", "synBuffSecond", "selfRegionMask",
              "nSecondHits", "dropInterleavesSmallerThan", "minRbhScore")
  for(i in params)
    cmb[[i]] <- get(i)

  cmb <- cmb[with(cmb, order(-(nGenes1 + nGenes2))),]
  gsParam$params$synteny <- cmb
  mxBuff <- max(gsParam$params$synteny$synBuff)
  if(is.null(arrayBuffer))
    arrayBuffer <- mxBuff
  arrayBuffer <- as.numeric(arrayBuffer[1])
  if(is.na(arrayBuffer) || arrayBuffer < mxBuff)
    arrayBuffer <- mxBuff
  gsParam$params$arrayBuffer <- arrayBuffer
  return(gsParam)
}


#' @title Run MCScanX from R
#' @description
#' \code{run_mcscanx} Internal GENESPACE function to run MCScanX on hits when
#' searching for synteny.
#' @rdname synteny
#' @import R.utils
#' @import data.table
#' @export
run_mcscanx <- function(hits,
                        blkSize,
                        nGaps,
                        gsParam,
                        path2mcscanx){
  ##############################################################################
  # parameter argument checking
  gn2 <- mcsID1 <- mcsID2 <- blkID <- ord2 <- ofID1 <- ofID2 <- ord1 <- NULL
  runBlast <- chr1 <- chr2 <- blkID <- NULL
  hitCols <- c("ofID1","ofID2","chr1","chr2","ord1","ord2","score")
  if(!all(hitCols %in% colnames(hits)))
    stop("Looks like the hits data.table is misformed. It must be a data.table",
         "comparing two genomes with the columns:", hitCols,"\n")
  if(length(blkSize) != 1)
    stop("blkSize must be a integer of length 1\n")
  if(length(nGaps) != 1)
    stop("blkSize must be a integer of length 1\n")
  if(!is.integer(blkSize)){
    blkSize <- as.integer(blkSize)
    if(is.na(blkSize))
      stop("Cannot coerce blkSize to an integer\n")
  }
  if(!is.integer(nGaps)){
    nGaps <- as.integer(nGaps)
    if(is.na(nGaps))
      stop("Cannot coerce nGaps to an integer\n")
  }
  MCScanX_h <- file.path(path2mcscanx, "MCScanX_h")
  if(!file.exists(MCScanX_h))
    stop("Cannot find MCScanX_h executable in", path2mcscanx,"\n")

  ##############################################################################
  # set tmp directory
  tmpDir <- file.path(
    gsParam$paths$results,
    paste0("tmp_",
           paste(sample(c(letters,LETTERS), 20, replace = T), collapse = "")))
  if (dir.exists(tmpDir))
    unlink(tmpDir, recursive = T)
  dir.create(tmpDir)
  on.exit(expr = unlink(tmpDir, recursive = T))

  ##############################################################################
  # convert gene locations to 'gff' like mcscanx file
  u1 <- subset(hits, !duplicated(ofID1))
  u2 <- subset(hits, !duplicated(ofID2))
  setkey(u1, chr1, ord1)
  setkey(u2, chr2, ord2)
  u1[,mcsID1 := paste0("aa", as.numeric(factor(chr1, levels = unique(chr1))))]
  u2[,mcsID2 := paste0("bb", as.numeric(factor(chr2, levels = unique(chr2))))]

  mcs1 <- u1[,c("mcsID1", "ofID1","ord1","ord1")]
  setnames(mcs1, c("chr", "id","start","end"))

  mcs2 <- u2[,c("mcsID2", "ofID2","ord2","ord2")]
  setnames(mcs2, c("chr", "id","start","end"))
  mcsGffIn <- rbind(mcs1, mcs2)

  ##############################################################################
  # convert hits to mcscanx_h blast input
  mcsBlsIn <- hits[,c("ofID1", "ofID2", "score")]
  # mcsBlsIn[,score := 1]

  mcsBlsIn[,ofID2 := paste0(ofID2,"xxxx")]
  mcsGffIn$id[grepl("bb", mcsGffIn$chr)] <- paste0(mcsGffIn$id[grepl("bb", mcsGffIn$chr)],"xxxx")

  blFile <- file.path(tmpDir, "mcs.homology")
  gfFile <- file.path(tmpDir, "mcs.gff")
  colFile <- file.path(tmpDir, "mcs.collinearity")

  fwrite(
    mcsGffIn,
    file = gfFile,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    showProgress = FALSE,
    verbose = FALSE)
  fwrite(
    mcsBlsIn,
    file = blFile,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    showProgress = FALSE,
    verbose = FALSE)

  ##############################################################################
  # run mcscanx_h
  mcsCom <- sprintf(
    "%s -a -b 2 -c 2 -m %s -s %s %s 1>/dev/null 2>&1",
    MCScanX_h, nGaps, blkSize, file.path(tmpDir,"mcs"))
  system(mcsCom)
  ##############################################################################
  idg <- strsplit(as.character(hits$ofID1[1]), "_")[[1]][1]
  # parse collinearity file
  suppressWarnings(collin <-  fread(
    cmd = sprintf("cat %s | grep %s_ | grep :", colFile, idg),
    col.names = c("blkID","gn1","gn2"),
    select = 1:3,
    header = F))
  if(nrow(collin) > 1){
    collin[,blkID := as.numeric(sapply(blkID, function(x)
      strsplit(x, "-")[[1]][1])) + 1]
    collin[,gn2 := gsub("xxxx$", "", gn2)]
    mcsb <- collin$blkID
    names(mcsb) <- with(collin, paste(gn1, gn2))
    return(mcsb)
  }
}

#' @title pull_synOGs
#' @description
#' \code{pull_synOGs} pull_synOGs
#' @rdname synteny
#' @import data.table
#' @importFrom grDevices pdf dev.off
#' @export
pull_synOGs <- function(gsParam,
                        genomeIDs = NULL){

  ##############################################################################
  # -- function to get syntenic ogs from blast hits
  add_synOg2gff <- function(gff,
                            hits = NULL,
                            gsParam,
                            genomeIDs,
                            allowRBHinOg,
                            useBlks){
    blkBuffer <- isSelf <- ofID1 <- ofID2 <- og <- ofID <- synOG <- blkID <- NULL
    # -- find hits files
    if(is.null(hits)){
      eg <- CJ(genomeIDs, genomeIDs)
      fs <- file.path(gsParam$paths$results,
                      sprintf("%s_%s_synHits.txt.gz",eg[[1]], eg[[2]]))
      fs <- fs[file.exists(fs)]

      # -- read all hits
      nCores <- gsParam$params$nCores
      hts <- rbindlist(mclapply(fs, mc.cores = nCores, function(i){
        if(useBlks){
          x <- fread(
            i, select = c("ofID1","ofID2","og","blkBuffer","blkAnchor","blkID"),
            na.strings = c("","NA"))
        }else{
          x <- fread(
            i, select = c("ofID1","ofID2","og","regBuffer","regAnchor","regID"),
            na.strings = c("","NA"),
            col.names = c("ofID1","ofID2","og","blkBuffer","blkAnchor","blkID"))
        }
        x <- subset(x, blkBuffer)
        x[,isSelf := any(ofID1 == ofID2), by = "blkID"]
        x <- subset(x, !isSelf & !is.na(og))
        if(!allowRBHinOg)
          x <- subset(x, !grepl("RBH", og))
        return(x[,c("ofID1", "ofID2", "blkAnchor", "blkID")])
      }))
    }else{
      hts <- data.table(hits)
      hts <- subset(hts, !is.na(og))
      if(!allowRBHinOg)
        hts <- subset(hts, !grepl("RBH", og))
      if(useBlks){
        hts <- hts[,c("ofID1", "ofID2", "blkBuffer", "blkAnchor", "blkID")]
      }else{
        hts <- with(hts, data.table(
          ofID1 = ofID1, ofID2 = ofID2,
          blkBuffer = regBuffer, blkAnchor = regAnchor, blkID = regID))
      }
      hts <- subset(hts, blkBuffer)
      hts[,blkBuffer := NULL]
    }

    # -- convert to syntenic orthogroups
    ic <- with(subset(hts, !is.na(blkID)), clus_igraph(
      id1 = c(ofID1, ofID2), id2 = c(ofID2, ofID1)))
    gff[,synOG := ic[ofID]]
    nmis <- sum(is.na(gff$synOG))
    mol <- max(gff$synOG, na.rm = T)
    gff$synOG[is.na(gff$synOG)] <- (mol + 1):(mol + nmis)
    gff[,synOG := as.integer(factor(synOG, levels = unique(synOG)))]
    return(gff)
  }

  ##############################################################################
  # -- function to combine the inblk and syntenic OGs
  combine_inblkSynOG <- function(genomeIDs,
                                 gff,
                                 gsParam){

    ofID <- ofID1 <- ofID2 <- clus <- og <- combOG <- inBlkOG <- synOG <- NULL
    if(gsParam$params$verbose)
      cat("Combining synteny-constrained and inblock orthogroups ...\n")

    if(gsParam$params$verbose)
      cat(sprintf("\tsyn OGs: %s, inblk OGs: %s",
                  uniqueN(gff$synOG, na.rm = T), uniqueN(gff$inBlkOG, na.rm = T)))
    noInBlk <- all(is.na(gff$inBlkOG))
    if(noInBlk)
      gff[,inBlkOG := synOG]

    inblk <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "inBlkOG"]
    syn <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "synOG"]
    u <- with(inblk, paste(ofID1, ofID2))
    syn <- subset(syn, !paste(ofID1, ofID2) %in% u & ofID1 != ofID2)
    tmp <- rbind(syn[,c("ofID1", "ofID2")], inblk[,c("ofID1", "ofID2")])
    tmp[,clus := clus_igraph(ofID1, ofID2)]
    ov <- with(tmp, c(clus, clus)); names(ov) <- with(tmp, c(ofID1, ofID2))
    gff[,combOG := ov[ofID]]
    nmis <- sum(is.na(gff$combOG))
    mol <- max(gff$combOG, na.rm = T)
    gff$combOG[is.na(gff$combOG)] <- (mol + 1):(mol + nmis)
    gff[,combOG := as.integer(factor(combOG, levels = unique(combOG)))]

    if(gsParam$params$verbose)
      cat(sprintf(", combined OGs: %s\n",
                  uniqueN(gff$combOG)))
    gff[,`:=`(og = combOG, combOG = NULL)]
    if(noInBlk)
      gff[,inBlkOG := NA]
    return(gff)
  }

  isArrayRep <- ogInblk <- ofID1 <- ofID2 <- inBlkOG <- ofID <- NULL
  og <- synOG <- globOG <- u <- arrayID <- NULL
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, na.strings = c("", "NA"), showProgress = F)
  verbose <- gsParam$params$verbose

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs

  # -- pull global syntenic orthogroups
  if(verbose)
    cat(sprintf(
      "Checking synteny-constrained global orthogroups for synOGs\n\tn. global OGs = %s\n",
      uniqueN(gff$globOG)))
  synGff <- add_synOg2gff(
    gff = subset(gff, isArrayRep),
    useBlks = F,
    gsParam = gsParam,
    genomeIDs = genomeIDs,
    allowRBHinOg = T)
  synGff[,og := globOG]
  # -- add the syntenic OGs to the original gff
  gff[,u := ifelse(is.na(arrayID), ofID, arrayID)]
  sog <- synGff$synOG
  names(sog) <- with(synGff, ifelse(is.na(arrayID), ofID, arrayID))
  gff[,`:=`(synOG = sog[u], u = NULL)]

  if(verbose)
    cat(sprintf("\tn. syntenic OGs = %s\n", uniqueN(gff$synOG)))

  # -- if desired, pull orthogroups within blocks
  if(gsParam$params$orthofinderInBlk){
    if(verbose)
      cat("Adding syntenic orthogroups from pairwise w/in-region hits\n")
    ofh <- blkwise_orthofinder(
      gsParam = gsParam,
      gff = synGff)
    ofh[, ogInblk := clus_igraph(ofID1, ofID2)]
    ofv <- ofh$ogInblk; names(ofv) <- ofh$ofID1
    synGff[,inBlkOG := ofv[ofID]]
    mol <- max(synGff$inBlkOG, na.rm = T)
    nmis <- sum(is.na(synGff$inBlkOG))
    synGff$inBlkOG[is.na(synGff$inBlkOG)] <- (mol + 1):(mol + nmis)
    synGff[,inBlkOG := as.integer(factor(inBlkOG, levels = unique(inBlkOG)))]

    gff[,u := ifelse(is.na(arrayID), ofID, arrayID)]
    bog <- synGff$inBlkOG
    names(bog) <- with(synGff, ifelse(is.na(arrayID), ofID, arrayID))
    gff[,`:=`(inBlkOG = bog[u], u = NULL)]
    gff[,og := inBlkOG]

    # -- combine inblk and syntenic OGs
    gff <- combine_inblkSynOG(
      genomeIDs = genomeIDs,
      gff = gff,
      gsParam = gsParam)
  }else{
    gff[,inBlkOG := NA]
    gff[,og := synOG]
  }
  return(gff)
}

#' @title add_arrayReps2gff
#' @description
#' \code{add_arrayReps2gff} add_arrayReps2gff
#' @rdname synteny
#' @import data.table
#' @importFrom grDevices pdf dev.off
#' @export
add_arrayReps2gff <- function(gff,
                             synBuff,
                             ogColumn,
                             verbose){

  pull_synArray <- function(gff, synBuff, ogColumn){
    # -- make global arrays from orthogroups
    gff[,arrayID := gff[[ogColumn]]]
    gff[,arrayID := sprintf("%s_%s_%s", genome, chr, arrayID)]
    gff[,n := .N, by = "arrayID"]

    # -- combine 1x ogs with ogs in regions < synBuff
    g1 <- subset(gff, n == 1)
    g2 <- subset(gff, n > 1)
    if(nrow(g2) > 1){
      g2[,rng := diff(range(ord)),  by = "arrayID"]
      g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
      g2 <- subset(g2, rng > synBuff)

      # -- combine above with ogs without max gap < synBuff
      g2[,rng := max(diff(ord[order(ord)])), by = "arrayID"]
      g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
      g2 <- subset(g2, rng > synBuff)

      # -- split ogs with gaps
      g2[,clus := dbscan(frNN(cbind(ord, ord), eps = synBuff), minPts = 1)$cluster,
         by = "arrayID"]
      g2[,arrayID := sprintf("%s_%s", arrayID, clus)]
      gff <- rbind(g1, g2[,colnames(g1), with = F])
    }
    gff[,n := NULL]
    return(gff)
  }

  pull_arrayRep <- function(gff){
    if(!"arrayID" %in% colnames(gff))
      stop("Cannot find a column named arrayID ... can't add arrayReps\n")

    if(!all(c("start", "ord", "pepLen") %in% colnames(gff)))
      stop("cannot find start, ord and/or pepLen in data.table columns\n")

    # -- split single and multi-member arrays
    gff[,n := .N, by = "arrayID"]
    out <- subset(gff, n == 1)
    out[,isArrayRep := TRUE]
    tmp <- subset(gff, n > 1)
    if(nrow(tmp) > 1){
      setkey(tmp, genome, ord)

      # -- calulate the distance to the median for each gene
      tmp[,`:=`(med = as.numeric(median(ord)),
                medbp = as.numeric(median(start)))]
      tmp[,`:=`(dist2med = abs(med - ord),
                dist2bp = abs(medbp - start))]

      # -- rank and choose representatives
      setorder(tmp, dist2med, dist2bp, -pepLen)
      tmp[,rnk := 1:.N, by = "arrayID"]
      tmp[,`:=`(isArrayRep = rnk == 1, rnk = NULL, dist2med = NULL, dist2bp = NULL,
                med = NULL, medbp = NULL)]

      # -- combine and return
      out <- rbind(out, tmp)
      out[,n := NULL]
      setkey(out, genome, ord)
    }
    return(out)
  }
  ovf <- gff$ord; names(ovf) <- gff$ofID

  gff <- pull_synArray(gff = gff, synBuff = synBuff, ogColumn = ogColumn)
  # -- add array representatives
  gff <- pull_arrayRep(gff)
  tmp <- data.table(gff)
  no <- 0; ni <- 1
  for(i in 1:10){
    if(ni > no){
      # -- get vector of gene order for array reps and convert order column
      ni <- uniqueN(tmp$arrayID)
      tmps <- subset(tmp, isArrayRep)
      setkey(tmps, genome, ord)
      tmps[,ord := 1:.N, by = "genome"]
      ov <- tmps$ord; names(ov) <- tmps$arrayID
      tmp[,ord := ov[arrayID]]

      # -- re-calculate arrayIDs with new order
      tmp <- pull_synArray(gff = tmp, synBuff = synBuff, ogColumn = ogColumn)
      tmp <- pull_arrayRep(tmp)
      no <- uniqueN(tmp$arrayID)
    }
  }
  gff <- data.table(tmp)

  if(verbose){
    cat("\tFound the following counts of arrays / genome:\n")
    spl <- subset(gff, arrayID %in% subset(gff, duplicated(arrayID))$arrayID)
    nu <- lapply(split(spl, by = "genome"), function(x)
      cat(sprintf("\t%s: %s genes in %s collinear arrays\n",
                  x$genome[1], nrow(x), uniqueN(x$arrayID))))
  }

  gff[,ord := ovf[ofID]]
  return(gff)
}

#' @title flag_synteny
#' @description
#' \code{flag_synteny} flag_synteny
#' @rdname synteny
#' @export
flag_synteny <- function(hits,
                         gff,
                         blkSize,
                         nGaps,
                         synBuff,
                         nhits1,
                         nhits2,
                         selfOnly = NULL,
                         onlyOgAnchors,
                         twoStep = F){

  ##############################################################################
  # -- function to pull hits within a buffer to any flagged with "isAnchor"
  find_hitsInBuff <- function(hits, synBuff){
    inBuffer <- ord1 <- ord2 <- isAnchor <- NULL
    if(!"isAnchor" %in% colnames(hits))
      stop("column name isAnchor must be in hits\n")
    splHit <- split(hits, by = c("gen1","gen2","chr1","chr2"))
    out <- rbindlist(lapply(splHit,  function(x){
      if(any(x$isAnchor)){
        nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synBuff))
        wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
        x[,inBuffer := 1:.N %in% wh]
        return(x)
      }
    }), fill = T)
    if(!"inBuffer" %in% colnames(out)){
      return(NULL)
    }else{
      return(subset(out, inBuffer))
    }
  }

  ##############################################################################
  # -- function to pull hits within block coordinates
  get_hitsInBlks <- function(blks, hits){
    ord1 <- ord2 <- blkID <- NULL
    splh <- split(hits, by = c("chr1","chr2","gen1","gen2"))
    splb <- split(blks, by = c("chr1","chr2","gen1","gen2"))
    splh <- splh[names(splb)]
    hits <- rbindlist(lapply(1:nrow(blks), function(i){
      y <- blks[i,]
      x <- splh[[with(y, paste(chr1, chr2, gen1, gen2, sep = "."))]]
      x <- subset(x, ord1 >= y$startOrd1 & ord1 <= y$endOrd1 &
                    ord2 >= y$minOrd2 & ord2 <= y$maxOrd2)
      x[,blkID := y$blkID]
      return(x)
    }))
    return(hits)
  }

  ##############################################################################
  # -- function to re-rank the scores of hits on a subset of the data
  rerank_score <- function(hits, byBlk = F){
    hits <- data.table(hits)
    if(byBlk){
      setorder(hits, ofID1, -score)
      hits[,scrRank1 := 1:.N, by = c("blkID","ofID1")]
      setorder(hits, ofID2, -score)
      hits[,scrRank2 := 1:.N, by = c("blkID","ofID2")]
    }else{
      setorder(hits, ofID1, -score)
      hits[,scrRank1 := 1:.N, by = "ofID1"]
      setorder(hits, ofID2, -score)
      hits[,scrRank2 := 1:.N, by = "ofID2"]
    }
    return(hits)
  }

  ##############################################################################
  # -- function to flag hits that are syntenic anchors
  flag_inBuffAnchor <- function(anchorhits,
                                allhits,
                                gsParam,
                                blkSize,
                                synBuff,
                                nGaps){
    h <- data.table(allhits)
    a <- data.table(anchorhits)
    # -- initial mcscan
    anchu <- run_mcscanx(
      hits = a,
      gsParam = gsParam,
      blkSize = blkSize,
      nGaps = nGaps,
      path2mcscanx = gsParam$paths$mcscanxCall)
    h[,isAnchor := paste(ofID1, ofID2) %in% names(anchu)]

    # -- pull hits around anchors
    a <- find_hitsInBuff(h, synBuff = synBuff)

    # -- initial anchor clustering
    a <- clus_dbscan(hits = a, radius = synBuff, blkSize = blkSize)
    anchv <- a$blkID; names(anchv) <- with(a, paste(ofID1, ofID2))
    h[,blkID := anchv[paste(ofID1, ofID2)]]
    h[,inBuffer := !is.na(blkID)]
    h[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2,
                           as.numeric(as.factor(blkID)))]
    h$blkID[with(h, !inBuffer)] <- NA
    return(h)
  }

  ##############################################################################
  # -- function to flag hits in synteny with self hits
  flag_selfSyn <- function(synBuff,
                           hits){
    tmp <- subset(hits, chr1 == chr2)
    tmp[,isAnchor := ofID1 == ofID2]
    tmp <- find_hitsInBuff(tmp, synBuff = synBuff)
    anchoru <- with(subset(tmp, isAnchor), paste(ofID1, ofID2))
    bufferu <- with(subset(tmp, inBuffer), paste(ofID1, ofID2))
    tmp[,`:=`(
      regID = sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, 1),
      regBuffer = paste(ofID1, ofID2) %in% bufferu,
      regAnchor = paste(ofID1, ofID2) %in% anchoru,
      blkID = sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, 1),
      blkBuffer = paste(ofID1, ofID2) %in% bufferu,
      blkAnchor = paste(ofID1, ofID2) %in% anchoru,
      lgRegID = sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, 1))]
    tmp$regID[with(tmp, !paste(ofID1, ofID2) %in% bufferu)] <- NA
    tmp$blkID[with(tmp, !paste(ofID1, ofID2) %in% bufferu)] <- NA
    tmp$lgRegID[with(tmp, !paste(ofID1, ofID2) %in% bufferu)] <- NA
    return(tmp)
  }

  ##############################################################################
  # -- function to flag anchors within regions
  flag_blksInReg <- function(gsParam,
                             hits,
                             blkSize,
                             nGaps){
    if(!all(c("regAnchor","regID") %in% colnames(hits)))
      stop("must have regions defined prior to running\n")

    # -- re-order within regions
    a <- subset(h, regAnchor)
    a[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense")),
      by = "regID"]

    # -- choose collinear hits in each region
    spl <- split(a, by = "regID")
    a <- rbindlist(lapply(spl, function(x){
      suppressWarnings(x[,isCollin := !is.na(run_mcscanx(
        hits = x, gsParam = gsParam, blkSize = blkSize, nGaps = nGaps,
        path2mcscanx = gsParam$paths$mcscanxCall)[paste(ofID1, ofID2)])])
      if(!"isCollin" %in% colnames(x))
        x[,isCollin := FALSE]
      return(x)
    }), use.names=TRUE)
    a <- subset(a, isCollin)
    a[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense")),
      by = "regID"]

    # -- cluster within regions to blocks
    spl <- split(a, by = "regID")
    iblks <- rbindlist(lapply(spl, function(x){
      x[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = sqrt((blkSize^2)*2)),
                         minPts = blkSize)$cluster]
      return(x)
    }))
    iblks <- subset(iblks, blkID != 0)

    # -- clustering within blks, again dropping blocks too small
    iblks[,blkID := sprintf("%s_%s_%s_%s",gen1, gen2, regID, blkID)]
    iblks[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                ord2 = frank(ord2, ties.method = "dense")),
          by = "blkID"]
    iblks[,nblk := dbscan(frNN(cbind(ord1, ord2), eps = sqrt(2)+.1),
                          minPts = blkSize)$cluster,
          by = "blkID"]
    fblks <- subset(iblks, blkID != 0)

    # -- final clustering within blks
    fblks[,blkID := sprintf("%s_%s", blkID, nblk)]
    fblks[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                ord2 = frank(ord2, ties.method = "dense")),
          by = "blkID"]
    fblks[,nblk := dbscan(frNN(cbind(ord1, ord2), eps = sqrt(2)+.1),
                          minPts = 0)$cluster,
          by = "blkID"]
    blv <- fblks$blkID; names(blv) <- with(fblks, paste(ofID1, ofID2))
    blAnc <- with(fblks, paste(ofID1, ofID2))

    # --
    hits[,blkID := blv[paste(ofID1, ofID2)]]
    hits[,blkID := ifelse(
      is.na(blkID), NA,
            sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2,
                    as.numeric(as.factor(blkID))))]
    return(hits)
  }

  ##############################################################################
  # -- function to find large regions
  flag_lgRegionv <- function(hits, synBuff){
    tmp <- subset(hits, regBuffer)
    tmp <- clus_dbscan(
      hits = tmp, radius = synBuff, blkSize = 1)
    lgv <- tmp$blkID; names(lgv) <- with(tmp, paste(ofID1, ofID2))
    hits[,lgRegID := lgv[paste(ofID1, ofID2)]]
    hits[,lgRegID := ifelse(
      is.na(lgRegID), NA,
      sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2,
              as.numeric(as.factor(lgRegID))))]
    return(hits)
  }

  ##############################################################################
  # -- function to pull hits that are in a buffer
  pull_anchBuffv <- function(hits, synBuff){
    tmp <- data.table(hits)
    bc <- calc_blkCoords(subset(tmp, !is.na(blkID)))
    tmp <- get_hitsInBlks(blks = bc, hits = tmp)
    tmp <- find_hitsInBuff(hits = tmp, synBuff = synBuff)
    return(tmp)
  }

  ##############################################################################
  ##############################################################################
  # -- Pipeline here   #########################################################
  # -- *NOTES*: 'h' is short for hits (all), 'a' is short for anchor hits
  # -- 'tmp' is reserved for intermediate calculatation that are then removed

  ##############################################################################
  # 1. Subset the data to hits that will always be used for analysis.
  # -- drop non-standard hits columns
  cols <- c("ofID1", "ofID2", "score", "gen1", "gen2", "start1", "start2",
            "end1", "end2", "chr1", "chr2", "ord1", "ord2", "scrRank1",
            "scrRank2", "isOg", "nChr1", "nChr2")
  hits <- hits[,cols, with = F]

  # -- keep only array reps
  tmp <- with(hits, unique(c(gen1, gen2)))
  genome1 <- tmp[1]
  genome2 <- tmp[length(tmp)]
  if(genome1 == genome2 && is.null(selfOnly)){
    selfOnly <- TRUE
  }else{
    selfOnly <- FALSE
  }

  # -- only ofIDs in the hits object and array reps
  gffRep <- subset(gff, isArrayRep)
  # -- make vectors of orders
  globOrdVect <- gff$ord; names(globOrdVect) <- gff$ofID
  repOrdVect <- gffRep$ord; names(repOrdVect) <- gffRep$ofID
  isArrayRepv <- gff$ofID[gff$isArrayRep]
  h <- subset(hits, ofID1 %in% isArrayRepv & ofID2 %in% isArrayRepv)

  # -- make sure all chrs are big enuf
  h <- subset(h, nChr1 >= blkSize & nChr2 >= blkSize)
  iv <- with(h, unique(c(ofID1, ofID2)))
  gffRep <- subset(gffRep, ofID %in% iv)

  # -- get order of genes in subset hits
  setkey(gffRep, genome, ord)
  gffRep[,ord := 1:.N, by = "genome"]

  # -- make vectors of orders
  globOrdVect <- gff$ord; names(globOrdVect) <- gff$ofID
  repOrdVect <- gffRep$ord; names(repOrdVect) <- gffRep$ofID
  isArrayRepv <- gff$ofID[gff$isArrayRep]

  ##############################################################################
  # 2. If self, parse accordingly
  if(genome1 == genome2 & selfOnly){
    h <- flag_selfSyn(
      hits = h, synBuff = synBuff)
  }else{
    ############################################################################
    # 3. parse hits as needed

    # -- overall re-ranking of positions
    h[,`:=`(ord1 = frank(repOrdVect[ofID1], ties.method = "dense"),
            ord2 = frank(repOrdVect[ofID2], ties.method = "dense"))]
    a <- data.table(h)

    # -- subset to orthogroups if desired
    if(onlyOgAnchors)
      a <- subset(a, isOg)

    # -- re-call score rank and subset accordingly
    a <- rerank_score(a)
    a <- subset(a, scrRank1 <= nhits1 & scrRank2 <= nhits1)

    # -- re-call order by array rep order
    a[,`:=`(ord1 = frank(repOrdVect[ofID1], ties.method = "dense"),
            ord2 = frank(repOrdVect[ofID2], ties.method = "dense"))]

    ############################################################################
    # 4. Pull initial anchors
    h <- flag_inBuffAnchor(
      anchorhits = a, allhits = h, gsParam = gsParam,
      blkSize = blkSize, nGaps = nGaps, synBuff = synBuff)

    # -- do again if 2step (add rbhs to ogs)
    if(twoStep){
      a <- subset(h, inBuffer)
      a <- rerank_score(hits = a, byBlk = T)
      a <- subset(a, isOg | (scrRank1 == 1 & scrRank2 == 1))

      a[,`:=`(ord1 = frank(repOrdVect[ofID1], ties.method = "dense"),
              ord2 = frank(repOrdVect[ofID2], ties.method = "dense"))]
      h <- flag_inBuffAnchor(
        anchorhits = a, allhits = h, gsParam = gsParam,
        blkSize = blkSize, nGaps = nGaps, synBuff = synBuff)
    }

    ############################################################################
    # 5. Finalize regions
    # -- Split overlapping regions
    tmp <- split_overlaps(
      hits = subset(h, isAnchor), dropInterleavesSmallerThan = 2,
      minPropDup = .05, verbose = F)
    bv <- tmp$blkID; names(bv) <- with(tmp, paste(ofID1, ofID2))
    h[,blkID := bv[paste(ofID1, ofID2)]]
    tmp <- NULL

    # -- final region assignment
    a <- pull_anchBuffv(hits = h, synBuff = synBuff)
    regv <- a$blkID; names(regv) <- with(a, paste(ofID1, ofID2))
    h[,`:=`(regID = regv[paste(ofID1, ofID2)],
            regAnchor = paste(ofID1, ofID2) %in% names(bv),
            regBuffer = paste(ofID1, ofID2) %in% names(regv),
            isAnchor = NULL, blkID = NULL)]

    # -- calculate large regions (for inblk ogs)
    h <- flag_lgRegionv(hits = h, synBuff = synBuff * 2)

    ############################################################################
    # 6. Make blocks (fine-scale within region synteny)
    h <- flag_blksInReg(
      gsParam = gsParam, hits = h, blkSize = blkSize, nGaps = nGaps)
    h[,isAnchor := !is.na(blkID)]

    # -- get hits in buffer from block anchors
    a <- pull_anchBuffv(hits = h, synBuff = synBuff)
    blkv <- a$blkID; names(blkv) <- with(a, paste(ofID1, ofID2))
    h[,`:=`(blkID = blkv[paste(ofID1, ofID2)],
            blkAnchor = isAnchor,
            blkBuffer = paste(ofID1, ofID2) %in% names(blkv),
            isAnchor = NULL, inBuffer = NULL)]
  }

  ##############################################################################
  # 7. Return data
  # -- merge with original hits
  colv <- c("ofID1", "ofID2", "regID" , "regBuffer", "regAnchor",
            "blkID", "blkBuffer", "blkAnchor", "lgRegID" )
  ho <- h[,colv, with = F]
  hitsOut <- merge(hits, ho, by = c("ofID1", "ofID2"), all.x = T)

  # -- setting NAs to FALSE where appropriate
  for(i in c("regBuffer", "regAnchor", "blkBuffer", "blkAnchor"))
    hitsOut[[i]][is.na(hitsOut[[i]])] <- FALSE
  return(hitsOut)
}

#' @title annotate_gff
#' @description
#' \code{annotate_gff} annotate_gff
#' @rdname synteny
#' @export
annotate_gff <- function(gsParam, genomeIDs){
  ##############################################################################
  # -- function to read in genespace-formatted gff files into memory
  read_gff <- function(gffFiles){
    genome <- NULL
    gff <- rbindlist(lapply(names(gffFiles), function(i){
      x <- fread(
        gffFiles[[i]],
        key = c("chr", "start","end","strand"))
      x[,genome := i]
      return(x)
    }))
    gffCols <- c("ord","chr","genome","start", "end", "strand","id")
    if(!all(gffCols %in% colnames(gff)))
      stop(paste(gffCols, collapse = ", "),
           " must all be column names in gff\n")
    return(gff)
  }

  ##############################################################################
  # -- function to add orthofinder ID to gff
  add_ofID2gff <- function(gff,
                           blastDir){
    id <- ofID <- genomeNum <- genome <- NULL
    specIDs <- read_orthofinderSpeciesIDs(blastDir)
    gv <- names(specIDs); names(gv) <- as.character(specIDs)
    seqIDs <- read_orthofinderSequenceIDs(blastDir)
    seqIDs[,genome :=  gv[as.character(genomeNum)]]
    idv <- seqIDs$ofID; names(idv) <- with(seqIDs, paste(genome, id))
    gff[,ofID := idv[paste(genome, id)]]
    return(gff)
  }

  ##############################################################################
  # -- function to add global orthofinder orthogroups to gff
  add_ogs2gff <- function(gsParam, gff){
    ogs <- parse_ogs(gsParam)
    gff <- merge(gff, ogs, by = c("genome","id"), all.x = T)
    gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
    setnames(gff, "ogID", "globOG")
    return(gff)
  }

  ##############################################################################
  # -- function to add peptide length to gff
  add_pepLen2gff <- function(gff,
                             gsParam){
    pepLen <- id <- ofID <-  NULL
    spl <- split(gff, by = "genome")
    naa <- rbindlist(lapply(spl, function(x){
      x[,pepLen := get_nAA(gsParam$paths$peptide[x$genome[1]], raw = T)[id]]
      return(x[,c("ofID","pepLen")])
    }))
    nao <- naa$pepLen; names(nao) <- naa$ofID
    gff[,pepLen := nao[ofID]]
    return(gff)
  }

  # -- read in the gff
  gff <- read_gff(gsParam$paths$gff)

  # -- add orthofinder ids
  gff <- add_ofID2gff(gff, gsParam$paths$blastDir)

  # -- subset to genomeIDs
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)

  # -- add peptide length
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)

  # -- add global orthogroups
  gff <- add_ogs2gff(gff = gff, gsParam = gsParam)

  return(gff)
}

#' @title pipe_synteny
#' @description
#' \code{pipe_synteny} pipe_synteny
#' @rdname synteny
#' @export
pipe_synteny <- function(gsParam,
                         gff,
                         nCores = NULL,
                         verbose = TRUE,
                         overwrite = FALSE,
                         useExistingSynHits = overwrite){

  ############################################################################
  # -- function to read in blast results
  parse_blast4synteny <- function(gsParam,
                                  genomeIDs,
                                  gff,
                                  ogColumn){

    score <- ofID1 <- ofID2 <- isOg <- og1 <- scrRank1 <- scrRank2 <- NULL
    chr1 <- chr2 <- ord1 <- ord2 <- NULL

    genome1 <- genomeIDs[1]
    genome2 <- genomeIDs[2]

    # -- get vectors from the gff
    if(!"isArrayRep" %in% colnames(gff))
      stop("must run add_arrayReps2gff\n")
    g <- subset(gff, genome %in% c(genome1, genome2))
    g[,og4syn := g[[ogColumn]]]
    gv <- g$genome; cv <- g$chr; ogv <- g$og4syn
    ov <- g$ord; sv <- g$start; ev <- g$end
    names(gv) <- names(cv) <- names(ogv) <- g$ofID
    names(ov) <- names(sv) <- names(ev) <- g$ofID

    # -- read in the blast hits
    ofsp <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
    if(genome1 == genome2){
      bl <- read_blast(
        path = gsParam$paths$blastDir,
        ofID1 = ofsp[genome1],
        ofID2 = ofsp[genome1])
    }else{
      bl <- rbind(read_blast(
        path = gsParam$paths$blastDir,
        ofID1 = ofsp[genome1],
        ofID2 = ofsp[genome2]),
        with(read_blast(
          path = gsParam$paths$blastDir,
          ofID1 = ofsp[genome2],
          ofID2 = ofsp[genome1]),
          data.table(
            ofID1 = ofID2, ofID2 = ofID1, score = score)))
    }

    # -- choose only the best scoring non-duplicated hits
    setorder(bl, -score)
    bl <- subset(bl, !duplicated(paste(ofID1, ofID2)))

    # -- add annotation information in
    bl[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
             start1 = sv[ofID1], start2 = sv[ofID2],
             end1 = ev[ofID1], end2 = ev[ofID2],
             chr1 = cv[ofID1], chr2 = cv[ofID2],
             ord1 = ov[ofID1], ord2 = ov[ofID2],
             scrRank1 = 1, scrRank2 = 1,
             isOg = ogv[ofID1] == ogv[ofID2])]

    bl <- subset(bl, !is.na(ord1) & !is.na(ord2))

    # -- get the number of unique hits by chr (for later culling)
    bl[,`:=`(nChr1 = uniqueN(ofID1[isOg]),
             nChr2 = uniqueN(ofID2[isOg])),
       by = c("gen1","gen2","chr1","chr2")]

    # -- optionally add the rank scores by geneID
    setorder(bl, ofID1, -score)
    bl[,scrRank1 := 1:.N, by = "ofID1"]
    setorder(bl, ofID2, -score)
    bl[,scrRank2 := 1:.N, by = "ofID2"]

    return(bl)
  }

  ############################################################################
  # -- function to split blocks by RLEs
  split_ovlBlks <- function(blk1, blk2, hits, dropInterleavesSmallerThan){
    blkID <- ofID1 <- ofID2 <- ord1 <- rl1 <- ord2 <- rl2 <- n <- NULL
    y <- subset(hits, blkID %in% c(blk1, blk2))
    y <- subset(y, !duplicated(ofID1))
    y <- subset(y, !duplicated(ofID2))
    setkey(y, ord1)
    y[,rl1 := add_rle(blkID)]
    setkey(y, ord2)
    y[,rl2 := add_rle(blkID)]
    y[,blkID := sprintf("%s_%s%s", blkID, rl1, rl2)]
    y[,n := .N, by = blkID]
    y <- subset(y, n >= dropInterleavesSmallerThan)
    y[,`:=`(rl1 = NULL, rl2 = NULL, n = NULL)]
    return(y)
  }

  ############################################################################
  # -- function to count and find overlapping blocks
  count_ovlHits <- function(hits, minPropDup){
    blkID <- i.blkID <- blk1 <- blk2 <- r1 <- r2 <- propDup <- n <- hasPriority <- NULL
    # -- calculate block coordinates
    hitstmp <- data.table(hits)
    bc <- calc_blkCoords(hitstmp)
    bc1 <- with(bc, data.table(
      blkID = blkID, chr = chr1, start = startOrd1, end = endOrd1,
      key = c("chr","start","end")))
    bc2 <- with(bc, data.table(
      blkID = blkID, chr = chr2, start = minOrd2, end = maxOrd2,
      key = c("chr","start","end")))

    # -- find non-self overlapping regions
    fo1 <- subset(foverlaps(bc1, bc1), blkID != i.blkID)
    fo2 <- subset(foverlaps(bc2, bc2), blkID != i.blkID)

    # -- add overlapping info into an output object
    out <- data.table(
      blk1 = c(fo1$blkID, fo2$blkID),
      blk2 = c(fo1$i.blkID, fo2$i.blkID))

    # -- choose the non-repeated overlaps
    u <- unique(c(out$blk1, out$blk2))
    un <- order(u); names(un) <- u
    out[,`:=`(r1 = un[blk1], r2 = un[blk2])]
    out[,u := ifelse(r1 > r2,  paste(blk1, blk2), paste(blk2, blk1))]
    out <- subset(out, !duplicated(u))
    out[,`:=`(r1 = NULL, r2 = NULL, u = NULL)]

    # -- count the number of duplicated hits between each pair of blocks
    if(nrow(out) > 0){
      tmp <- split(subset(hitstmp, blkID %in% unlist(out)), by = "blkID")
      out[,propDup := sapply(1:nrow(out), function(j){
        tb1 <- tmp[[out$blk1[j]]]
        tb2 <- tmp[[out$blk2[j]]]
        tot <- uniqueN(c(tb1$ofID1, tb1$ofID2, tb2$ofID1, tb2$ofID2))
        dup <- sum(tb1$ofID1 %in% tb2$ofID1, tb1$ofID2 %in% tb2$ofID2,
                   tb2$ofID1 %in% tb1$ofID1, tb2$ofID2 %in% tb1$ofID2)
        return(dup/tot)
      })]
      out <- subset(out, propDup <= minPropDup)
      if(nrow(out) > 0){
        # -- get the total number of hits in each block
        out[,`:=`(n1 = sapply(tmp[blk1], nrow),
                  n2 = sapply(tmp[blk2], nrow))]

        # -- prioritize by the pairs with the largest small block
        out[,n := apply(.SD, 1, min), .SDcols = c("n1", "n2")]

        setorder(out, -n)
        out[,n := NULL]
        v <- matrix(duplicated(as.character(unlist(t(out[,1:2])))), ncol = 2, byrow = T)
        out[,hasPriority := !apply(v, 1, any)]
        return(out)
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  }

  ##############################################################################
  # -- Function to split overlapping regions / blocks
  split_overlaps <- function(hits,
                             dropInterleavesSmallerThan,
                             minPropDup = 0.05,
                             maxIter = 20,
                             verbose){
    hasPriority <- blkID <- NULL
    anchs <- data.table(hits)
    nov <- 1
    iter <- 0
    while(nov > 0 && iter <= maxIter){
      iter <- iter + 1
      if(verbose)
        cat(sprintf("\titer %s: %s hits / %s blks\n",
                    iter, nrow(anchs), uniqueN(anchs$blkID)))
      ovlBlks <- count_ovlHits(anchs, minPropDup = minPropDup)
      if(!is.null(ovlBlks)){
        ovlBlks <- subset(ovlBlks, hasPriority)
        nov <- nrow(ovlBlks)
        if(is.null(nov))
          nov <- 0
        if(nov > 0){
          splAnch <- rbindlist(lapply(1:nrow(ovlBlks), function(j)
            split_ovlBlks(
              blk1 = ovlBlks$blk1[j],
              blk2 = ovlBlks$blk2[j],
              hits = anchs,
              dropInterleavesSmallerThan)))
          if(!is.null(splAnch)){
            anchs <- rbind(
              splAnch,
              subset(anchs, !blkID %in% unique(unlist(ovlBlks[,c(1:2)]))))
          }else{
            nov <- 0
          }
        }else{
          nov <- 0
        }
      }else{
        nov <- 0
      }
    }
    if(verbose)
      cat("\nfound no additional overlapping blocks\n")
    anchs[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2,
                            as.numeric(as.factor(blkID)))]
    return(anchs)
  }

  ##############################################################################
  # -- Function to write the syntenic hits
  write_synhits <- function(hits, gsParam){
    genome1 <- hits$gen1[1]
    genome2 <- hits$gen2[1]
    hitsFile <- file.path(gsParam$paths$results,
                          sprintf("%s_%s_synHits.txt.gz", genome1, genome2))
    fwrite(hits, file = hitsFile, quote = F, showProgress = F, sep = "\t")
  }

  ##############################################################################
  # -- Function to plot three plots
  plot_synHits <- function(hits, gsParam){
    tp <- data.table(hits)

    pm <- par()["mfrow"]
    par(mfrow = c(1,1))
    pmar <- par()["mar"]
    par(mar = c(3,2,2,2))

    pdf(file.path(gsParam$paths$results,
                  sprintf("%s_%s_dotplots.pdf", tp$gen1[1], tp$gen2[1])),
        height = 6, width = 6)
    par(mfrow = c(1,1))

    plot_hits(hits = tp, plotType = "allOG",  cols = "blue2", round2 = 10)
    plot_hits(hits = tp, plotType = "regBuffer")
    plot_hits(hits = tp, plotType = "blkAnchor")
    dev.off()

    par(mfrow = pm)
    par(mar = pmar)
  }

  ##############################################################################
  # -- Function to split up the synParam for chunk-wise parallelization
  split_synParam2chunks <- function(gsParam,
                                    nCores = NULL,
                                    overwrite = FALSE){
    if(is.null(nCores))
      nCores <- gsParam$params$nCores
    synp <- subset(data.table(gsParam$params$synteny), runBlast)

    # -- check if the synHits files are there
    if(!overwrite){
      synp[,synHitsFile := file.path(
        gsParam$paths$results,
        sprintf("%s_%s_synHits.txt.gz", genome1, genome2))]
      synp[,hasRes := file.exists(synHitsFile)]
      synp <- subset(synp, !hasRes)
      synp[,hasRes := NULL]
      chunks <- NULL
    }
    if(nrow(synp) > 0){
      synp[,wt := ifelse(
        genome1 == genome2 & (nSecondHits1 + nSecondHits2 + ploidy1 + ploidy2) == 1,
        nGenes1, nGenes1 + nGenes2)]
      setorder(synp, -wt)
      synp[,chunk := rep(1:nrow(synp), each = nCores)[1:.N]]
      chunks <- split(synp, by = "chunk")
    }
    return(chunks)
  }

  ##############################################################################
  ##############################################################################
  # 1. Split the synparamters into chunks based on nCores
  if(is.na(gsParam$paths$orthogroupsDir)){
    if(verbose)
      cat("Indexing location of orthofinder results ... ")
    gsParam <- find_orthofinderResults(gsParam)
    if(verbose)
      cat("Done!\n")
  }

  if(is.null(nCores))
    nCores <- gsParam$params$nCores
  splSynp <- split_synParam2chunks(
    gsParam, nCores = nCores, overwrite = overwrite)

  if(is.null(splSynp)){
    cat("All synHits results have been generated and !overwrite\n\tNot running initial synteny step\n")
  }else{
    # -- build synHits file by chunk
    if(verbose)
      cat(sprintf("Pulling synteny for %s unique pairwise combinations of genomes\n",
                  sum(sapply(splSynp, nrow))))
    if(verbose)
      cat(sprintf("\tRunning %s chunks of %s combinations each:\n",
                  length(splSynp), nCores))
    blks <- rbindlist(lapply(1:length(splSynp), function(i){
      if(verbose)
        cat(sprintf("\tChunk %s / %s (%s) ... ",
                    i,length(splSynp), format(Sys.time(), "%X")))
      chunkOut <- mclapply(1:nrow(splSynp[[i]]), mc.cores = nCores, function(j){
        x <- splSynp[[i]][j,]

        # -- check if we should run 2step
        run2step <- with(x, (genome1 != genome2 & (ploidy1 + ploidy2) > 2))

        # -- read in the hits
        synHitsFile <- with(x, file.path(
          gsParam$paths$results,
          sprintf("%s_%s_synHits.txt.gz", genome1, genome2)))
        if(useExistingSynHits && file.exists(synHitsFile)){
          allHits <- fread(synHitsFile, sep = "\t", na.strings = c("NA", ""))
        }
        allHits <- with(x, parse_blast4synteny(
          gsParam = gsParam, genomeIDs = c(genome1, genome2),
          gff = gff, ogColumn = "globOG"))

        # -- find the hits in synteny
        synHits <- with(x, flag_synteny(
          hits = allHits, gff = gff, nhits1 = nhits1, nhits2 = nhits2, nGaps = nGaps,
          synBuff = synBuff, blkSize = blkSize,
          onlyOgAnchors = onlyOgAnchors, twoStep = run2step))

        # -- check if we need to run again masking above
        runAgain <- with(
          x, nSecondHits1 > 0 | nSecondHits2 > 0 |
            (genome1 == genome2 & (ploidy1 + ploidy2) > 2))

        # -- if so, do it again, need to flag new blocks etc.
        if(runAgain){
          ns1 <- with(x, (ploidy1 - 1) + nSecondHits1)
          ns2 <- with(x, (ploidy2 - 1) + nSecondHits2)
          if(ns1 < 1) ns1 <- 1
          if(ns2 < 1) ns2 <- 2
          rrHits <- subset(synHits, !regBuffer)[,colnames(allHits), with = F]
          z2 <- with(x, flag_synteny(
            hits = rrHits, gff = gff, nhits1 = ns1, nhits2 = ns2,
            nGaps = nGapsSecond, synBuff = synBuffSecond, blkSize = blkSizeSecond,
            selfOnly = FALSE, onlyOgAnchors = onlyOgAnchorsSecond, twoStep = T))
          z2[,`:=`(regID = paste0(regID, "_second"),
                   blkID = paste0(blkID, "_second"),
                   lgRegID = paste0(lgRegID, "_second"))]
          synHits <- rbind(subset(synHits, !regBuffer), z2)
        }

        # -- make the plots
        plot_synHits(gsParam = gsParam, hits = synHits)

        # -- write the hits and return
        write_synhits(gsParam = gsParam, hits = synHits)
        return(synHits)
      })
      if(verbose)
        cat("Done!\n")
      # -- Print the summaries
      tmp <- rbindlist(lapply(chunkOut, function(x){
        if(verbose)
          cat(sprintf(
            "\t%s-%s: %s (tot), %s/%s (reg), %s/%s (blk)\n",
            pull_strWidth(x$gen1[i], 8), pull_strWidth(x$gen2[i], 8),  nrow(x),
            sum(x$regBuffer), uniqueN(x$regID[x$regAnchor]),
            sum(x$blkAnchor), uniqueN(x$blkID[x$blkAnchor])))
      }))
      return(bcout)
    }))
  }
  return(gsParam)
}

