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
#' @param minRbhScore integer specifying the minimum blast bit score to allow
#' for a RBH to be included
#' @param path2mcscanx character string file.path pointing to the install
#' directory for the MCScanX program. In particular, there must be an
#' executable in this directory for MCScanX_h.
#' @param overwrite logical, should existing directories be overwritten?
#' @param minGenes4of integer of length 1, specifying the minimum number of
#' unique genes on each sequence for orthofinder to be run
#' @param gff data.table containing the gff-like annotatin information
#' @param ogColumn character length 1, specifying the column in the gff
#' data.table to use as the orthogroup vector
#' @param nCores integer of length 1, specifying the number of parallel
#' processes to run
#' @param verbose logical, should updates be printed to the console
#' @param nhits1 integer of length 1, specifying the top n hits to retain for
#' each gene on genome 1
#' @param nhits2 integer of length 1, specifying the top n hits to retain for
#' each gene on genome 2
#' @param selfOnly logical, should only self hits be used? This is much faster.
#' @param maskTheseHits character vector paste of the pairwise blast orthofinder
#' gene IDs. Used to mask hits in the hits object
#' @param maxIter integer of length 1 specifying the maximum number of
#' iterations
#' @param blks data.table containing the block coordinates. See calc_blkCoords
#' @param radius numeric of length 1 specifying the search radius.
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
#' Synteny produces two measures of collinearity between pairs of genomes:
#' the collinear block, flagged as "blkID" and the large syntenic region,
#' flagged as "regID". The syntenic anchors for each of these types of synteny
#' are flagged in the logical "isAnchor" column (TRUE). Hits within the syntenic
#' buffer that are not anchors are flagged with the appropriate regID and blkID.
#' Both anchors and those in the syntenic buffer are flagges in the logical
#' column "inBuffer" (TRUE).
#' both blkIDs and regIDs and
#'
#' Parameterization of synteny is complex and varies by the type of
#' contrast desired. To simplify parameterization, we have build a convenience
#' function that infers various parameters based on basic
#' genespace parameters and ploidy. Set set_syntenyParams for more information
#' about the input synParam data.table.
#'
#' @return Synteny functions produce three types of outputs: (1) a 'hits'
#' data.table, (2) a 'blocks' data.table, and (3) the genespace parameter list
#' 'gsParam' in tutorials. Here, are is the metadata for the first two outputs.
#'
#' Hits output stored as results/$genome1_$genome2_synHits.txt.gz. Columns names
#' ending in '1' are the target gene and '2' are the query.
#' \enumerate{
#' \item ofID1/2: orthofinder unique ID for each gene.
#' \item score: the diamond blast-like bit score
#' \item gen1/2: genome IDs, as specified in init_genespace
#' \item start1/2: gff-derived gene start position in basepairs
#' \item end1/2: gff-derived gene end position in basepairs
#' \item chr1/2: gff-derived chromosome IDs
#' \item ord1/2: gene-rank order across the genome
#' \item scrRank1/2: within-gene score rank, higher score = lower rank
#' \item arrayID1/2: collinear array ID, defined as genes in the same orthogroup
#' within the synteny buffer within a chromosome.
#' \item arrayOrd1/2: the condensed order of genes, by arrayID
#' \item isRep1/2: is gene1/2 a representative of its array, defined as the most
#' centrally located gene within the array and the longest peptide (tiebreaker)
#' \item isOg: are both query and target in the same orthogroup?
#' \item nChr1/2: number of unique genes on chr1/2.
#' \item blkID: syntenic block ID
#' \item regID: large syntenic region ID (regions may be interleafed)
#' \item isAnchor: is this hit a syntenic block anchor?
#' \item inBuffer: is this hit within the synBuffer of an anchor?
#' }
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
                    overwrite = F,
                    minGenes4of = 40,
                    ...){

  ##############################################################################
  # 1. Set up the environment and parameter checking
  gp <- gsParam
  gsParam <- NULL

  # -- set data.table threads to 1
  setDTthreads(1)

  # -- genome IDs
  if(is.null(genomeIDs))
    genomeIDs <- gp$genomes$genomeIDs
  if(!any(genomeIDs %in% gp$genomes$genomeIDs))
      stop("specified genomeIDs dont look right\n")

  # -- shortcuts and output files
  verbose <- gp$params$verbose
  writeTo <- gp$paths$results
  nCores <- gp$params$nCores
  gffFile <- file.path(writeTo, "gffWithOgs.txt.gz")
  blksFile <- file.path(writeTo, "syntenicBlocks.txt.gz")

  # -- set the synteny parameters
  if(is.data.table(gp$params$synteny))
    if(!all(genomeIDs %in% gp$params$synteny$genome1))
      gp$params$synteny <- NULL
  if(!is.data.table(gp$params$synteny)){
    gp <- set_syntenyParams(gp, ...)
    if(identical(gp$params$synteny, gp$params$synteny)){
      cat("Synteny Parameters have not been set! Setting to defaults\n")
    }else{
      cat("Synteny Parameters have not been set! Using user-defined settings\n")
      gp$params$synteny <- gp$params$synteny
    }
  }

  # -- find the orthogroups
  if(is.na(gp$paths$orthogroupsDir)){
    if(verbose)
      cat("Indexing location of orthofinder results ... ")
    gp <- find_orthofinderResults(gp)
    if(verbose)
      cat("Done!\n")
  }

  ##############################################################################
  # 2. load and parse the gff
  # -- add in global syntenic orthogroup arrays, orthofinder IDs, etc
  if(verbose)
    cat("Parsing the gff files ... \n\tReading the gffs and adding orthofinder IDs ... ")
  gf <- annotate_gff(
    gsParam = gp,
    genomeIDs = genomeIDs)
  if(verbose)
    cat(sprintf(
      "Done!\n\tFound %s global OGs for %s genes\n",
      uniqueN(gf$globOG), nrow(gf)))

  # -- QC gff
  ng <- gf[,list(n = .N, nog = uniqueN(globOG)), by = c("genome", "chr")]
  np <- ng[,list(nchrPass = sum(nog >= 5), nogPass = sum(nog[nog >= 5]),
                 nchrFail = sum(nog < 5), nogFail = sum(nog[nog < 5]),
                 ngenePass = sum(n[nog >= 5]), ngeneFail = sum(n[nog < 5])),
           by = "genome"]
  np[,rat := nogPass/(nogPass + nogFail)]
  cat("\tQC-ing genome to ensure chromosomes/scaffolds are big enough...
      \t\tGenome: n. chrs PASS/FAIL, n. genes PASS/FAIL, n. OGs PASS/FAIL\n")
  for(i in 1:nrow(np))
    with(np[i,], cat(sprintf(
      "\t\t%s: %s/%s, %s/%s, %s/%s\n",
      genome, nchrPass, nchrFail, ngenePass, ngeneFail, nogPass, nogFail)))
  if(any(np$rat < 0.9)){
    warning("Some genomes have a large number of genes on scaffolds with < 5 genes
            these may be problematic for synteny\n")
  }else{
    cat("\tAll look good!\n")
  }

  ##############################################################################
  # -- add arrays. This takes the specifed column in the gff and parses array
  # and array representatives. We do this again below if inblk OG is used.
  # With verbose = T, prints the counts of reps and total members / genome.
  if(verbose)
    cat("\tDefining collinear orthogroup arrays ... \n")
  gf <- add_arrayReps2gff(
    gff = gf,
    synBuff = max(gp$params$synteny$synBuff)+1,
    ogColumn = "globOG",
    verbose = verbose)

  ##############################################################################
  # 3. Run the initial synteny builder
  # -- This is the full pipeline: split up synParams into chunks each with
  # nCores entries, read in blast results, parse to collinear hits, cull to
  # syntenic regions, form large regions, build syntenic blocks within regions.
  gp <- pipe_synteny(
    gsParam = gp,
    gff = gf,
    overwrite = overwrite)

  ##############################################################################
  # 4. Build syntenic orthogroups from syntenic hits
  if(verbose)
    cat("Defining synteny-constrained orthogroups ... \n")
  gf <- add_synOg2gff(
    gff = gf,
    gsParam = gp,
    genomeIDs = genomeIDs)
  if(verbose)
    cat(sprintf(
      "\tFound %s synteny-split OGs for %s genes\n",
      uniqueN(gf$synOG), nrow(gf)))

  ##############################################################################
  # 5. If runOrthofinderInBlk, do so here:
  if(!gp$params$orthofinderInBlk){
    synOG <- inblkOG <- NULL
    gf[,`:=`(inblkOG = synOG, og = synOG)]
  }else{
    gf <- blkwise_orthofinder(
      gsParam = gp,
      gff = gf,
      genomeIDs = genomeIDs,
      minGenes4of = minGenes4of)

    # -- make an new OG column with syntenic and inBlk orthogroups combined
    gf <- combine_inblkSynOG(
      genomeIDs = genomeIDs,
      gff = gf,
      gsParam = gp)
    isArrayRep <- arrayID <- NULL
    gf[,`:=`(isArrayRep = NULL, arrayID = NULL)]

    # -- re-call arrays with the new og column
    gf <- add_arrayReps2gff(
      gff = gf,
      synBuff = max(gp$params$synteny$synBuff)+1,
      ogColumn = "og",
      verbose = verbose)
  }

  # -- write the final gff
  fwrite(
    gf, file = gffFile, sep = "\t",
    quote = FALSE, showProgress = FALSE)
  if(verbose)
    cat(sprintf(
      "Found %s OGs across %s genes. gff3-like text file written to:\n\t%s\n",
      uniqueN(gf$og), nrow(gf), gffFile))

  ##############################################################################
  # 8. Calculate block coordinates
  if(verbose)
    cat("Calculating syntenic block breakpoints ... \n")
  runBlast <- isAnchor <- blkID <- NULL
  synHitsFiles <- with(subset(gp$params$synteny, runBlast), file.path(
    writeTo, sprintf("%s_%s_synHits.txt.gz", genome1, genome2)))

  blks <- rbindlist(mclapply(synHitsFiles, mc.cores = nCores, function(i)
    calc_blkCoords(subset(
      fread(i, sep = "\t", na.strings = c("NA", ""), showProgress = F),
      isAnchor & !is.na(blkID)), mirror = T)))

  fwrite(
    blks, file = blksFile, sep = "\t",
    quote = F, showProgress = F)
  if(verbose)
    cat(sprintf(
      "\tFound %s blocks. Text file written to:\n\t%s:\n",
      nrow(blks), blksFile))

  return(gp)
}

#' @title pipe_synteny
#' @description
#' \code{pipe_synteny} pipe_synteny
#' @rdname synteny
#' @importFrom grDevices pdf dev.off
#' @export
pipe_synteny <- function(gsParam,
                         gff,
                         ogColumn = "globOG",
                         nCores = NULL,
                         verbose = TRUE,
                         overwrite = FALSE){
  ##############################################################################
  # -- Function to write the syntenic hits
  write_synhits <- function(hits,
                            gsParam){
    setDTthreads(1)
    genome1 <- hits$gen1[1]
    genome2 <- hits$gen2[1]
    hitsFile <- file.path(gsParam$paths$results, sprintf(
      "%s_%s_synHits.txt.gz",  genome1, genome2))
    fwrite(
      hits, file = hitsFile, quote = F, showProgress = F, sep = "\t")
  }

  ##############################################################################
  # -- Function to plot three plots
  plot_synHits <- function(hits,
                           gsParam){
    setDTthreads(1)
    tp <- data.table(hits)

    pdf(file.path(gsParam$paths$results, sprintf(
      "%s_%s_dotplots.pdf",
      tp$gen1[1], tp$gen2[1])),
      height = 6, width = 6)

    plot_hits(hits = tp, plotType = "allOG",  cols = "blue2", round2 = 50)
    plot_hits(hits = tp, plotType = "regAnchor", useOrder = F)
    plot_hits(hits = tp, plotType = "blkAnchor")
    dev.off()
  }

  ##############################################################################
  # -- Function to split up the synParam for chunk-wise parallelization
  split_synParam2chunks <- function(gsParam,
                                    nCores = NULL,
                                    overwrite = FALSE){
    runBlast <- synHitsFile <- wt <- nGenes2 <- nGenes2 <- chunk <- og4syn <-
      u <- NULL
    setDTthreads(1)
    if(is.null(nCores))
      nCores <- gsParam$params$nCores
    synp <- subset(data.table(gsParam$params$synteny), runBlast)

    # -- check if the synHits files are there
    genome1 <- genome2 <- NULL
    if(!overwrite){
      synp[,synHitsFile := file.path(
        gsParam$paths$results,
        sprintf("%s_%s_synHits.txt.gz", genome1, genome2))]
      synp[,hasRes := file.exists(synHitsFile)]
      hasRes <- NULL
      synp <- subset(synp, !hasRes)
      synp[,hasRes := NULL]
      chunks <- NULL
    }

    if(nrow(synp) > 0){
      nGenes1 <- genome1 <- genome2 <-  nGenes1 <- ploidy2 <- ploidy1 <-
        nSecondHits2 <- nSecondHits1 <- wt <- chunk <- NULL
      synp[,wt := ifelse(
        genome1 == genome2 &
          (nSecondHits1 + nSecondHits2 + ploidy1 + ploidy2) == 1,
        nGenes1, nGenes1 + nGenes2)]
      setorder(synp, -wt)
      synp[,chunk := rep(1:nrow(synp), each = nCores)[1:.N]]
      chunks <- split(synp, by = "chunk")
    }
    return(chunks)
  }

  ##############################################################################
  # -- function to read in blast results
  parse_blast4synteny <- function(gsParam,
                                  genomeIDs,
                                  gff,
                                  ogColumn){
    setDTthreads(1)
    score <- ofID1 <- ofID2 <- isOg <- og1 <- scrRank1 <- scrRank2 <- chr1 <-
      chr2 <- ord1 <- ord2 <- og4syn <- NULL

    genome1 <- genomeIDs[1]
    genome2 <- genomeIDs[2]

    # -- get vectors from the gff
    if(!"isArrayRep" %in% colnames(gff))
      stop("must run add_arrayReps2gff\n")
    genome <- NULL
    g <- subset(gff, genome %in% c(genome1, genome2))
    g[,og4syn := g[[ogColumn]]]
    gv <- g$genome; cv <- g$chr; ogv <- g$og4syn; ov <- g$ord; av <- g$arrayID
    rv <- g$isArrayRep; sv <- g$start; ev <- g$end; arv <- g$arrayOrd
    names(rv) <- names(gv) <- names(cv) <- names(ogv) <- names(av) <- g$ofID
    names(arv) <- names(ov) <- names(sv) <- names(ev) <- g$ofID

    # -- read in the blast hits
    ofsp <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
    ofID1 <- ofID2 <- NULL
    if(genome1 == genome2){
      bl <- read_blast(
        path = gsParam$paths$blastDir,
        ofID1 = ofsp[genome1],
        ofID2 = ofsp[genome1])
    }else{
      bl <- rbind(read_blast(
        path = gsParam$paths$blastDir,
        ofID1 = ofsp[genome1], ofID2 = ofsp[genome2]),
        with(read_blast(path = gsParam$paths$blastDir,
                        ofID1 = ofsp[genome2],ofID2 = ofsp[genome1]),
             data.table(ofID1 = ofID2, ofID2 = ofID1, score = score)))
    }

    # -- choose only the best scoring non-duplicated hits
    ofID1 <- ofID2 <- score <- NULL
    setorder(bl, -score)
    bl <- subset(bl, !duplicated(paste(ofID1, ofID2)))

    # -- add annotation information in
    bl[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2], start1 = sv[ofID1],
             start2 = sv[ofID2], end1 = ev[ofID1], end2 = ev[ofID2],
             chr1 = cv[ofID1], chr2 = cv[ofID2], ord1 = ov[ofID1],
             ord2 = ov[ofID2], scrRank1 = 1, scrRank2 = 1,
             arrayID1 = av[ofID1], arrayID2 = av[ofID2],
             arrayOrd1 = arv[ofID1], arrayOrd2 = arv[ofID2],
             isRep1 = rv[ofID1], isRep2 = rv[ofID2],
             isOg = ogv[ofID1] == ogv[ofID2])]
    ord1 <- ord2 <- NULL
    bl <- subset(bl, !is.na(ord1) & !is.na(ord2))

    # -- get the number of unique hits by chr (for later culling)
    ofID1 <- ofID2 <- isOg <- NULL
    bl[,`:=`(nChr1 = uniqueN(ofID1[isOg]),
             nChr2 = uniqueN(ofID2[isOg])),
       by = c("gen1","gen2","chr1","chr2")]

    # -- optionally add the rank scores by geneID
    ofID1 <- ofID2 <- score <- NULL
    setorder(bl, ofID1, -score)
    bl[,scrRank1 := 1:.N, by = "ofID1"]
    setorder(bl, ofID2, -score)
    bl[,scrRank2 := 1:.N, by = "ofID2"]
    bl[,u := paste(ofID1, ofID2)]
    return(bl)
  }

  setDTthreads(1)
  ##############################################################################
  # 1. Get the workspace set up
  # -- name of output columns
  hitsNames <- c(
    "ofID1","ofID2","score","gen1","gen2","start1","start2","end1","end2",
    "chr1","chr2","ord1","ord2","scrRank1","scrRank2","arrayID1","arrayID2",
    "arrayOrd1","arrayOrd2","isRep1","isRep2","isOg","nChr1","nChr2",
    "blkID","regID","isAnchor","inBuffer")

  # -- load the orthogroup directory if needed
  if(is.na(gsParam$paths$orthogroupsDir)){
    if(verbose)
      cat("Indexing location of orthofinder results ... ")
    gsParam <- find_orthofinderResults(gsParam)
    if(verbose)
      cat("Done!\n")
  }

  # -- set the number of cores
  if(is.null(nCores) || is.na(as.numeric(nCores)))
    nCores <- gsParam$params$nCores
  nCores <- nCores[1]

  # -- split of the synteny parameters by chunks
  splSynp <- split_synParam2chunks(
    gsParam, nCores = nCores, overwrite = overwrite)
  u <- NULL
  if(is.null(splSynp)){
    cat("All synHits results have been generated and !overwrite\n\tNot running initial synteny step\n")
  }else{
    ############################################################################
    # -- build synHits file by chunk
    if(verbose)
      cat(sprintf("Pulling synteny for %s unique pairwise combinations of genomes\n",
                  sum(sapply(splSynp, nrow))))
    if(verbose)
      cat(sprintf("\tRunning %s chunks of %s combinations each:\n",
                  length(splSynp), nCores))
    blks <- lapply(1:length(splSynp), function(i){
      setDTthreads(1)
      if(verbose)
        cat(sprintf("\tChunk %s / %s (%s) ... ",
                    i,length(splSynp), format(Sys.time(), "%X")))
      chunkOut <- mclapply(1:nrow(splSynp[[i]]), mc.cores = nCores, function(j){
        setDTthreads(1)
        ########################################################################
        # 2. read in the hits
        x <- splSynp[[i]][j,]
        g1 <- x$genome1
        g2 <- x$genome2

        # -- read in the hits
        synHitsFile <- with(x, file.path(
          gsParam$paths$results,
          sprintf("%s_%s_synHits.txt.gz", g1, g2)))

        allHits <- parse_blast4synteny(
          gsParam = gsParam,
          genomeIDs = c(g1, g2),
          gff = gff,
          ogColumn = ogColumn)
        slf <- msk <- pri <- sec <- hme <- NULL
        if(g1 == g2){

          ######################################################################
          # 3. self only regions
          blkID <- regID <- inBuffer <- NULL
          slf <- flag_synteny(
            hits = allHits,
            gsParam = gsParam,
            blkSize = 1,
            synBuff = x$synBuff,
            nhits1 = 1,
            nhits2 = 1,
            nGaps = 1,
            selfOnly = TRUE,
            maskTheseHits = NULL,
            onlyOgAnchors = TRUE)
          slf <- subset(slf, inBuffer & !is.na(blkID) & !is.na(regID))
          slf[,`:=`(blkID = paste0("self_", blkID),
                    regID = paste0("self_", regID))]

          ######################################################################
          # 4. if ploidy > 1 or second hits needed mask self to larger region
          if(x$ploidy1 > 1 | x$nSecondHits > 0){
            inBuffer <- NULL
            msk <- flag_synteny(
              hits = allHits,
              gsParam = gsParam,
              blkSize = 1,
              synBuff = 500,
              nhits1 = 1,
              nhits2 = 1,
              nGaps = 1,
              selfOnly = TRUE,
              maskTheseHits = NULL,
              onlyOgAnchors = TRUE)
            msk <- subset(msk, inBuffer)
          }

          ######################################################################
          # 5. if polyploid, get homeologs
          if(x$ploidy1 > 1){
            onlyOgAnchors <- inBuffer <- inBuffer <- blkID <- regID <-
              nGaps <- nhits2 <- nhits1 <- synBuff <- blkSize <- NULL
            hme <- with(x, flag_synteny(
              hits = allHits,
              gsParam = gsParam,
              blkSize = blkSize,
              synBuff = synBuff,
              nhits1 = nhits1 - 1,
              nhits2 = nhits2 - 1,
              nGaps = nGaps,
              selfOnly = FALSE,
              maskTheseHits = msk,
              onlyOgAnchors = onlyOgAnchors))
            hme <- subset(hme, inBuffer)
            hme[,`:=`(
              blkID = ifelse(!is.na(blkID), paste0("homeo_", blkID), NA),
              regID = ifelse(!is.na(regID), paste0("homeo_", regID), NA))]
          }

          ######################################################################
          # 6. if second hits, get paralogs
          if(x$nSecondHits > 0){
            msk <- rbind(msk, hme)
            blkID <- onlyOgAnchorsSecond <-  nGapsSecond <- nSecondHits <-
              synBuffSecond <- blkSizeSecond <- inBuffer <- regID <- NULL
            sec <- with(x, flag_synteny(
              hits = allHits,
              gsParam = gsParam,
              blkSize = blkSizeSecond,
              synBuff = synBuffSecond,
              selfOnly = FALSE,
              nhits1 = nSecondHits,
              nhits2 = nSecondHits,
              maskTheseHits = msk,
              nGaps = nGapsSecond,
              onlyOgAnchors = onlyOgAnchors))
            sec <- subset(sec, inBuffer)
            sec[,`:=`(
              blkID = ifelse(!is.na(blkID), paste0("second_", blkID), NA),
              regID = ifelse(!is.na(regID), paste0("second_", regID), NA))]
          }

        }else{
          ######################################################################
          # 7. if intergenomics get synteny
          onlyOgAnchors <-  nGaps <- nhits2 <- nhits1 <- synBuff <- blkSize <-
            blkID <- regID <- inBuffer <- NULL
          pri <- with(x, flag_synteny(
            hits = allHits,
            gsParam = gsParam,
            blkSize = blkSize,
            synBuff = synBuff,
            nhits1 = nhits1,
            nhits2 = nhits2,
            nGaps = nGaps,
            selfOnly = FALSE,
            maskTheseHits = NULL,
            onlyOgAnchors = onlyOgAnchors))
          pri <- subset(pri, inBuffer)
          pri[,`:=`(
            blkID = ifelse(!is.na(blkID), paste0("primary_", blkID), NA),
            regID = ifelse(!is.na(regID), paste0("primary_", regID), NA))]

        ########################################################################
        # 8. if second hits, mask and rerun
          if(x$nSecondHits > 0){
            nGaps <- nhits2 <- nhits1 <- synBuff <- blkSize <-
              blkID <- regID <- inBuffer <- NULL
            sec <- with(x, flag_synteny(
              hits = allHits,
              gsParam = gsParam,
              blkSize = blkSizeSecond,
              synBuff = synBuffSecond,
              selfOnly = FALSE,
              nhits1 = nSecondHits,
              nhits2 = nSecondHits,
              maskTheseHits = pri,
              nGaps = nGapsSecond,
              onlyOgAnchors = onlyOgAnchorsSecond))
            sec <- subset(sec, inBuffer)
            sec[,`:=`(
              blkID = ifelse(!is.na(blkID), paste0("second_", blkID), NA),
              regID = ifelse(!is.na(regID), paste0("second_", regID), NA))]
          }
        }

        ########################################################################
        # 9. Process and plot

        # -- make sure there is synteny
        synHits <- rbind(slf, pri, hme, sec)
        if(nrow(synHits) < 1)
          with(synHits, stop(
            gen1[1], " vs. ", gen2[1], " does not have any detectable synteny ...
            either relax synteny parameters or drop one of these genomes\n"))

        # -- combine hits
        nHits <- subset(allHits, !u %in% synHits$u)
        nHits[,`:=`(blkID = NA, isAnchor = FALSE, regID = NA, inBuffer = FALSE)]
        synHits <- rbind(synHits, nHits)
        if(!all(hitsNames %in% colnames(synHits)))
          stop("hits object does not have the right column names\n")
        synHits <- synHits[,hitsNames, with = F]

        # -- make the plots
        plot_synHits(
          gsParam = gsParam,
          hits = synHits)

        # -- write the hits and return
        write_synhits(
          gsParam = gsParam,
          hits = synHits)

        return(synHits)
      })
      if(verbose)
        cat("Done!\n")

      # -- Print the summaries
      if(verbose)
        tmp <- lapply(chunkOut, function(x)
          cat(sprintf(
            "\t%s-%s: %s (tot), %s/%s (reg), %s/%s (blk)\n",
            pull_strWidth(x$gen1[i], 8),
            pull_strWidth(x$gen2[i], 8),
            nrow(x),
            sum(x$inBuffer),
            uniqueN(x$regID[x$inBuffer]),
            sum(x$isAnchor),
            uniqueN(x$blkID[x$isAnchor]))))
    })
  }
  #`1`
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
                              onlyOgAnchorsSecond = FALSE,
                              blkSize = 5,
                              blkSizeSecond = blkSize,
                              nGaps = 5,
                              nGapsSecond = nGaps*4,
                              nSecondHits = 0,
                              synBuff = 100,
                              synBuffSecond = synBuff,
                              selfRegionMask = synBuff * 5,
                              minRbhScore = 50,
                              arrayBuffer = sqrt(2*(synBuff^2))){
  setDTthreads(1)
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
  minRbhScore <- as.integer(minRbhScore)[1]
  if(is.null(minRbhScore) || is.na(minRbhScore) || length(minRbhScore) == 0)
    stop("cannot coerce minRbhScore to integer\n")

  # -- shorten some names
  genomeIDs <- gsParam$genomes$genomeIDs
  ploidy <- gsParam$genomes$ploidy[genomeIDs]
  blastDir <- gsParam$paths$orthofinder
  pepFiles <- gsParam$paths$peptide

  # -- Check the peptide files
  if(!all(file.exists(pepFiles)))
    stop(sprintf(
      "the following genome annotations have not been correctly parsed:\n\t%s\n\tRerun parse_annotations\n",
      paste(names(pepFiles)[which(!file.exists(pepFiles))], collapse = "\n\t")))

  # -- Calculate the number of genes for each genome
  nGenes <- sapply(pepFiles, function(x) length(readAAStringSet(x)))
  if(min(nGenes) == 0)
    stop(sprintf(
      "the following genome annotations have not been correctly parsed:\n\t%s\n\tRerun parse_annotations\n",
      paste(names(pepFiles)[which(nGenes == 0)], collapse = "\n\t")))

  # make the database of unique combinations of genomes
  cmb <- data.table(
    genome1 = names(ploidy),
    genome2 = names(ploidy))
  cmb <- cmb[,CJ(genome1, genome2)]
  setnames(cmb, c("genome1", "genome2"))
  # -- add ploidy, orthofinder IDs, n genes
  genome1 <- genome2 <- NULL
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
  genome1 <- genome2 <- nGenes1 <- nGenes2 <- NULL
  cmb[,tiebreak := as.numeric(factor(genome1, levels = genomeIDs)) -
        as.numeric(factor(genome2, levels = genomeIDs))]
  cmb[,`:=`(query = ifelse(nGenes1 > nGenes2, genome1,
                           ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                  genome1, genome2)),
            target = ifelse(nGenes1 > nGenes2, genome2,
                            ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                   genome2, genome1))),
      by = "u"]

  # -- flag which combinations should be run for blast etc.
  genome1 <- query <- target <- tiebreak <- runBlast <- NULL
  cmb[,`:=`(runBlast = genome1 == query,
            mirrorBlast = genome1 == target & genome1 != query,
            tiebreak = NULL)]

  # -- choose the number of secondary hits
  nhits2 <- nhits1 <- genome1 <- genome2 <- nSecondHits1 <- nSecondHits2 <-
    genome1 <- genome2 <- NULL
  cmb <<- cmb
  cmb[,nSecondHits1 := ifelse(
    nhits1 == 1 & nhits2 == 1, nSecondHits,
    ifelse(genome1 == genome2 & nhits2 > 1 & nSecondHits > 0, (nSecondHits * nhits2/2) + (nhits2/2),
           ifelse(nhits2 > 1 & nSecondHits > 0, nSecondHits * nhits2, nSecondHits)))]
  cmb[,nSecondHits2 := ifelse(
    nhits1 == 1 & nhits2 == 1, nSecondHits,
    ifelse(genome1 == genome2 & nhits1 > 1 & nSecondHits > 0, (nSecondHits * nhits1/2) + (nhits1/2),
           ifelse(nhits1 > 1 & nSecondHits > 0, nSecondHits * nhits1, nSecondHits)))]

  ##############################################################################
  # 3. Add all the parameters into the output data.table.
  ##############################################################################
  params <- c(
    "onlyOgAnchors","onlyOgAnchorsSecond","blkSize", "blkSizeSecond",
    "nGaps", "nGapsSecond", "synBuff", "synBuffSecond", "selfRegionMask",
    "nSecondHits", "minRbhScore")
  for(i in params)
    cmb[[i]] <- get(i)

  nGenes1 <- nGenes2 <- NULL
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
  setDTthreads(1)
  mcsID1 <- mcsID2 <- NULL
  ##############################################################################
  # parameter argument checking
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
    paste0("tmp_", paste(
      sample(c(letters,LETTERS), 20, replace = T),
      collapse = "")))
  if(dir.exists(tmpDir))
    unlink(tmpDir, recursive = T)
  dir.create(tmpDir)
  on.exit(expr = unlink(tmpDir, recursive = T))

  ##############################################################################
  # convert gene locations to 'gff' like mcscanx file
  ord1 <- ord2 <- chr1 <- chr2 <- ofID1 <- ofID2 <- NULL
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
    MCScanX_h,
    nGaps,
    blkSize,
    file.path(tmpDir,"mcs"))
  system(mcsCom)

  ##############################################################################
  # parse collinearity file
  idg <- strsplit(as.character(hits$ofID1[1]), "_")[[1]][1]
  suppressWarnings(collin <- fread(
    cmd = sprintf("cat %s | grep %s_ | grep :", colFile, idg),
    col.names = c("blkID","gn1","gn2"),
    select = 1:3,
    showProgress = FALSE,
    header = FALSE))

  if(nrow(collin) > 1){
    blkID <- gn1 <- gn2 <- NULL
    collin[,blkID := as.numeric(sapply(blkID, function(x)
      strsplit(x, "-")[[1]][1])) + 1]
    collin[,gn2 := gsub("xxxx$", "", gn2)]
    mcsb <- collin$blkID
    names(mcsb) <- with(collin, paste(gn1, gn2))
    return(mcsb)
  }
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
    setDTthreads(1)
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
    setDTthreads(1)
    id <- ofID <- genomeNum <- genome <- NULL
    specIDs <- read_orthofinderSpeciesIDs(blastDir)
    gv <- names(specIDs)
    names(gv) <- as.character(specIDs)

    genomeNum <- NULL
    seqIDs <- read_orthofinderSequenceIDs(blastDir)
    seqIDs[,genome :=  gv[as.character(genomeNum)]]
    idv <- seqIDs$ofID

    id <- genome <- NULL
    names(idv) <- with(seqIDs, paste(genome, id))
    gff[,ofID := idv[paste(genome, id)]]
    return(gff)
  }

  ##############################################################################
  # -- function to add global orthofinder orthogroups to gff
  add_ogs2gff <- function(gsParam, gff){
    setDTthreads(1)
    ogs <- parse_ogs(gsParam)
    gff <- merge(
      gff, ogs,
      by = c("genome","id"), all.x = T)
    gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
    setnames(gff, "ogID", "globOG")
    return(gff)
  }

  ##############################################################################
  # -- function to add peptide length to gff
  add_pepLen2gff <- function(gff,
                             gsParam){
    setDTthreads(1)
    id <- NULL
    spl <- split(gff, by = "genome")
    naa <- rbindlist(lapply(spl, function(x){
      x[,pepLen := get_nAA(gsParam$paths$peptide[x$genome[1]], raw = T)[id]]
      return(x[,c("ofID","pepLen")])
    }))
    nao <- naa$pepLen; names(nao) <- naa$ofID
    pepLen <- ofID <- NULL
    gff[,pepLen := nao[ofID]]
    return(gff)
  }

  # -- read in the gff
  gff <- read_gff(gsParam$paths$gff)

  # -- add orthofinder ids
  gff <- add_ofID2gff(gff, gsParam$paths$blastDir)

  # -- subset to genomeIDs
  genome <- ord <- NULL
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)

  # -- add peptide length
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)

  # -- add global orthogroups
  gff <- add_ogs2gff(gff = gff, gsParam = gsParam)

  return(gff)
}

#' @title flag_synteny
#' @description
#' \code{flag_synteny} flag_synteny
#' @rdname synteny
#' @export
flag_synteny <- function(hits,
                         gsParam,
                         blkSize,
                         synBuff,
                         nhits1,
                         nhits2,
                         nGaps,
                         selfOnly,
                         maskTheseHits = NULL,
                         onlyOgAnchors){
  setDTthreads(1)
  ##############################################################################
  # 1. subset the hits to those that could potentially be block anchors
  # -- "raw" will never be subset and gets columns added to it in situ
  raw <- data.table(hits)
  isAnchor <- blkID <- inBuffer <- arrayID1 <- arrayID2 <- regID <- NULL
  suppressWarnings(raw[,`:=`(blkID = NULL, isAnchor = NULL, inBuffer = NULL)])

  # -- "anchs" is the initial dataset that contains hits to be used for anchors
  # this never gets subset either
  # -- anchors must be array representatives
  isRep1 <- isRep2 <- NULL
  anchs <- subset(hits, isRep1 & isRep2)

  # -- optionally mask some genes
  # these genes are excluded from potentially being syntenic anchors
  ofID1 <- ofID2 <- u <- NULL
  if(is.data.table(maskTheseHits)){
    uMask <- with(maskTheseHits, unique(paste(ofID1, ofID2)))
    anchs <- subset(anchs, !u %in% uMask)
  }

  # -- re-rank the order based on the array rep order
  arrayOrd1 <- arrayOrd2 <- NULL
  anchs[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
              ord2 = frank(arrayOrd2, ties.method = "dense"))]

  # -- make a copy of anchs. This gets subset all the time.
  a <- data.table(anchs)

  ##############################################################################
  # 2. If self only, just return the vectors of anchors, buffer and blks based
  # on the chromosomes
  if(selfOnly & a$gen1[1] == a$gen2[1]){
    ofID1 <- ofID2 <- arrayOrd1 <- arrayOrd2 <- chr1 <- inBuffer <-
      isAnchor <- NULL
    raw[,`:=`(isAnchor = ofID1 == ofID2,
              inBuffer = abs(arrayOrd1 - arrayOrd2) <= synBuff |
                arrayID1 == arrayID2,
              blkID = as.numeric(as.factor(chr1)),
              regID = as.numeric(as.factor(chr1)))]
    raw$blkID[with(raw, !inBuffer | chr1 != chr2)] <- NA
    raw$regID[with(raw, !inBuffer | chr1 != chr2)] <- NA
  }else{

    ############################################################################
    # 3. Otherwise, do the full synteny finder ...
    # First off, parse and re-rank the "a" anchors data.table
    # -- if only OG anchors, subset to these
    isOg <- NULL
    if(onlyOgAnchors)
      a <- subset(a, isOg)

    # -- pull top n scoring hits
    ofID1 <- ofID2 <- score <- scrRank1 <- scrRank2 <- NULL
    setorder(a, ofID1, -score)
    a[,scrRank1 := 1:.N, by = "ofID1"]
    setorder(a, ofID2, -score)
    a[,scrRank2 := 1:.N, by = "ofID2"]
    if(nhits1 < 1) nhits1 <- 1
    if(nhits2 < 1) nhits2 <- 1
    a <- subset(a, scrRank1 <= nhits1 & scrRank2 <= nhits2)

    # -- compress order of culled array reps
    arrayOrd1 <- arrayOrd2<- NULL
    a[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
            ord2 = frank(arrayOrd2, ties.method = "dense"))]

    # -- mcscan pruning to anchors
    anchu <- run_mcscanx(
      nGaps = nGaps,
      gsParam = gsParam,
      blkSize = blkSize,
      hits = a,
      path2mcscanx = gsParam$paths$mcscanxCall)

    if(is.null(anchu)){
      raw[,`:=`(blkID = NA, inBuffer = FALSE, isAnchor = FALSE, regID = NA)]
    }else{
      ##########################################################################
      # 4. If there is synteny, pull hits close to anchors and re-do  anchors
      # this is so that, if ploidy is miss-specified, the blocks don't get
      # chopped up

      # -- keep only chrs with synteny
      chr1 <- chr2 <- NULL
      uchrs <- with(subset(a, u %in% names(anchu)), unique(paste(chr1, chr2)))
      a <- subset(anchs, paste(chr1, chr2) %in% uchrs)

      # -- optionally drop non-OG hits
      isOg <- NULL
      if(onlyOgAnchors)
        a <- subset(anchs, isOg)

      # -- pull hits proximate to the anchors
      u <- NULL
      anchs[,isAnchor := u %in% names(anchu)]
      a <- find_hitsInBuff(hits = anchs, synBuff = synBuff)

      # -- mcscan pruning to anchors
      anchu <- run_mcscanx(
        nGaps = nGaps,
        gsParam = gsParam,
        blkSize = blkSize,
        hits = a,
        path2mcscanx = gsParam$paths$mcscanxCall)

      if(is.null(anchu)){
        raw[,`:=`(blkID = NA, inBuffer = FALSE, isAnchor = FALSE, regID = NA)]
      }else{
        ########################################################################
        # 4. If there is synteny, conduct anchor clustering
        # -- round 1 - large radius
        u <- ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- chr1 <- chr2 <- blkID <- NULL
        a[,isAnchor := u %in% names(anchu)]
        a <- subset(a, isAnchor)
        a[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                ord2 = frank(arrayOrd2, ties.method = "dense"))]
        a[,blkID := dbscan(frNN(
          x = cbind(ord1, ord2),
          eps = synBuff),
          minPts = blkSize)$cluster,
          by = c("chr1", "chr2")]
        a <- subset(a, blkID > 0)
        a[,blkID := as.numeric(as.factor(paste(chr1, chr2, blkID)))]

        # -- round 2 - recalculate order within blocks and re-calculate anchors
        u <- ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- blkID <- NULL
        a[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                ord2 = frank(arrayOrd2, ties.method = "dense")),
          by = "blkID"]
        spl <- split(a, by = "blkID")
        anchu <- unlist(lapply(spl, function(inBlk)
          names(run_mcscanx(
            nGaps = nGaps,
            gsParam = gsParam,
            blkSize = blkSize,
            hits = inBlk,
            path2mcscanx = gsParam$paths$mcscanxCall))))
        a <- subset(a, u %in% anchu)

        # -- round 3 - recluster anchors with large radius in block
        chr1 <- chr2 <- u <- ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- blkID <- NULL
        a[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                ord2 = frank(arrayOrd2, ties.method = "dense"))]
        a[,blkID := dbscan(frNN(
          x = cbind(ord1, ord2),
          eps = synBuff),
          minPts = blkSize)$cluster,
          by = c("chr1", "chr2")]
        a <- subset(a, blkID > 0)
        a[,blkID := as.numeric(as.factor(paste(chr1, chr2, blkID)))]

        # -- round 4 - cluster hits within blocks using small block radius
        tmp <- ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- blkID <- NULL
        a[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                ord2 = frank(arrayOrd2, ties.method = "dense")),
          by = "blkID"]
        a[,tmp := dbscan(frNN(
          x = cbind(ord1, ord2),
          eps = blkSize),
          minPts = blkSize)$cluster,
          by = "blkID"]
        a <- subset(a, tmp > 0)
        a[,blkID := as.numeric(as.factor(paste(tmp, blkID)))]
        a[,tmp := NULL]

        # -- make sure that the unique gene number is right
        ofID1 <- ofID2 <- n1 <- n2 <- NULL
        a[,`:=`(n1 = uniqueN(ofID1), n2 = uniqueN(ofID2)), by = "blkID"]
        a <- subset(a, n1 >= blkSize & n2 >= blkSize)
        a[,`:=`(n1 = NULL, n2 = NULL)]

        # -- finalize anchor-block IDs by splitting overlapping non-duplicated
        # blocks
        ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- NULL
        a[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                ord2 = frank(arrayOrd2, ties.method = "dense"))]
        b <- split_blks(hits = a, blkSize = blkSize, maxIter = 5)
        b[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                ord2 = frank(arrayOrd2, ties.method = "dense"))]
        ########################################################################
        # 5. Add all block and buffer information back into the hits object
        # -- BLOCKS: add block anchor IDs to the raw data
        blkv <- b$blkID; names(blkv) <- b$u
        u <- NULL
        raw[,`:=`(blkID = blkv[u], isAnchor = u %in% names(blkv))]
        tmp <- data.table(raw)

        # -- calculate block coordinates
        ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- blkID <- NULL
        tmp[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                  ord2 = frank(arrayOrd2, ties.method = "dense"))]
        bc <- calc_blkCoords(subset(tmp, !is.na(blkID)))

        # -- pull hits within synteny buffer in the block coordinates
        u <- NULL
        tmp <- get_hitsInBlks(blks = bc, hits = tmp)
        tmp <- find_hitsInBuff(hits = tmp, synBuff = synBuff)
        blkv <- tmp$blkID; names(blkv) <- tmp$u
        raw[,blkID := blkv[u]]

        # -- REGIONS
        # -- cluster anchors
        u <- isAnchor <- ord1 <- ord2 <- arrayOrd1 <- arrayOrd2 <- blkID <- chr1 <- chr2 <- NULL
        tmp <- subset(raw, isAnchor)
        tmp[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
                  ord2 = frank(arrayOrd2, ties.method = "dense"))]
        tmp[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = synBuff),
                           minPts = blkSize)$cluster, by = c("chr1", "chr2")]
        tmp <- subset(tmp, blkID > 0)
        tmp[,blkID := as.numeric(as.factor(paste(chr1, chr2, blkID)))]
        blkv <- tmp$blkID; names(blkv) <- tmp$u
        raw[,regID := blkv[u]]

        # -- get buffer hits for large regions
        ord1 <- ord2 <- arrayOrd1 <- arrayOrd2<- NULL
        tmp <- data.table(raw)
        tmp[,`:=`(blkID = regID,
                  ord1 = frank(arrayOrd1, ties.method = "dense"),
                  ord2 = frank(arrayOrd2, ties.method = "dense"))]

        # -- pull hits bounded by regions
        blkID <- u <- NULL
        bc <- calc_blkCoords(subset(tmp, !is.na(blkID)))
        tmp <- get_hitsInBlks(blks = bc, hits = tmp)
        tmp <- find_hitsInBuff(hits = tmp, synBuff = synBuff)
        blkv <- tmp$blkID; names(blkv) <- tmp$u
        raw[,`:=`(regID = blkv[u], inBuffer = u %in% tmp$u)]
      }
    }
  }
  return(raw)
}

#' @title add_arrayReps2gff
#' @description
#' \code{add_arrayReps2gff} add_arrayReps2gff
#' @rdname synteny
#' @export
add_arrayReps2gff <- function(gff,
                              synBuff,
                              ogColumn,
                              verbose){
  pull_synArray <- function(gff, synBuff, ogColumn){
    rng <- NULL
    setDTthreads(1)
    # -- make global arrays from orthogroups
    gff[,arrayID := gff[[ogColumn]]]
    genome <- chr <- arrayID <- NULL
    gff[,arrayID := sprintf("%s_%s_%s", genome, chr, arrayID)]
    gff[,n := .N, by = "arrayID"]

    # -- combine 1x ogs with ogs in regions < synBuff
    n <- NULL
    g1 <- subset(gff, n == 1)
    g2 <- subset(gff, n > 1)

    if(nrow(g2) > 1){
      ord <- NULL
      g2[,rng := diff(range(ord)),  by = "arrayID"]
      g2[,clus := 1]
      g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
      g2 <- subset(g2, rng > synBuff)
      if(nrow(g2) > 1){
        # -- combine above with ogs without max gap < synBuff
        ord <- NULL
        g2[,rng := max(diff(ord[order(ord)])), by = "arrayID"]
        g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
        g2 <- subset(g2, rng > synBuff)

        # -- split ogs with gaps
        ord <- arrayID <- clus <- NULL
        g2[,clus := dbscan(frNN(cbind(ord, ord), eps = synBuff), minPts = 1)$cluster,
           by = "arrayID"]
      }
      g2[,arrayID := sprintf("%s_%s", arrayID, clus)]
      gff <- rbind(g1, g2[,colnames(g1), with = F])
    }
    gff[,n := NULL]
    return(gff)
  }

  pull_arrayRep <- function(gff){
    med <- medbp <- rnk <- NULL
    setDTthreads(1)
    if(!"arrayID" %in% colnames(gff))
      stop("Cannot find a column named arrayID ... can't add arrayReps\n")

    if(!all(c("start", "ord", "pepLen") %in% colnames(gff)))
      stop("cannot find start, ord and/or pepLen in data.table columns\n")

    # -- split single and multi-member arrays
    n <- genome <- ord <- NULL
    gff[,n := .N, by = "arrayID"]
    out <- subset(gff, n == 1)
    out[,isArrayRep := TRUE]
    tmp <- subset(gff, n > 1)
    if(nrow(tmp) > 1){
      setkey(tmp, genome, ord)

      # -- calulate the distance to the median for each gene
      ord <- start <- NULL
      tmp[,`:=`(med = as.numeric(median(ord)),
                medbp = as.numeric(median(start)))]
      tmp[,`:=`(dist2med = abs(med - ord),
                dist2bp = abs(medbp - start))]

      # -- rank and choose representatives
      dist2med <- dist2bp <- pepLen <- NULL
      setorder(tmp, dist2med, dist2bp, -pepLen)
      tmp[,rnk := 1:.N, by = "arrayID"]
      tmp[,`:=`(isArrayRep = rnk == 1, rnk = NULL, dist2med = NULL,
                dist2bp = NULL, med = NULL, medbp = NULL)]

      # -- combine and return
      out <- rbind(out, tmp)
      out[,n := NULL]
      genome <- ord <- NULL
      setkey(out, genome, ord)
    }
    return(out)
  }

  # -- get collinear array IDs
  setDTthreads(1)
  arrayOrd <- NULL
  ovf <- gff$ord; names(ovf) <- gff$ofID
  gff <- data.table(gff)
  gff <- pull_synArray(
    gff = gff,
    synBuff = synBuff,
    ogColumn = ogColumn)

  # -- add array representatives
  gff <- pull_arrayRep(gff)
  tmp <- data.table(gff)
  no <- 0; ni <- 1
  for(i in 1:10){
    if(ni > no){
      # -- get vector of gene order for array reps and convert order column
      ni <- uniqueN(tmp$arrayID)
      isArrayRep <- genome <- ord <- arrayID <- NULL
      tmps <- subset(tmp, isArrayRep)
      setkey(tmps, genome, ord)
      tmps[,ord := 1:.N, by = "genome"]
      ov <- tmps$ord; names(ov) <- tmps$arrayID
      tmp[,ord := ov[arrayID]]

      # -- re-calculate arrayIDs with new order
      tmp <- pull_synArray(
        gff = tmp,
        synBuff = synBuff,
        ogColumn = ogColumn)
      tmp <- pull_arrayRep(tmp)
      no <- uniqueN(tmp$arrayID)
    }
  }
  gff <- data.table(tmp)

  if(verbose){
    cat("\tFound the following counts of arrays / genome:\n")
    arrayID <- NULL
    spl <- subset(gff, arrayID %in% subset(gff, duplicated(arrayID))$arrayID)
    nu <- lapply(split(spl, by = "genome"), function(x)
      cat(sprintf("\t\t%s: %s genes in %s collinear arrays\n",
                  x$genome[1], nrow(x), uniqueN(x$arrayID))))
  }

  # -- adding in the array order here
  ofID <- genome <- ord <- arrayID <- NULL
  gff[,ord := ovf[ofID]]
  setkey(gff, genome, ord)
  gff[,arrayOrd := as.numeric(factor(arrayID, levels = unique(arrayID))),
      by = "genome"]

  # -- lower-memory coding of the arrays and OGs
  arrayID <- globOG <- NULL
  gff[,arrayID := as.integer(factor(arrayID, levels = unique(arrayID)))]
  gff[,globOG := as.integer(factor(globOG, levels = unique(globOG)))]
  return(gff)
}

#' @title add synteny-constrained orthogroups to gff
#' @description
#' \code{add_synOg2gff} add synteny-constrained orthogroups to gff
#' @rdname synteny
#' @importFrom parallel mclapply
#' @export
add_synOg2gff <- function(gff,
                          gsParam,
                          genomeIDs){

  setDTthreads(1)
  # -- find hits files
  eg <- CJ(genomeIDs, genomeIDs)
  fs <- file.path(gsParam$paths$results,
                  sprintf("%s_%s_synHits.txt.gz",eg[[1]], eg[[2]]))
  fs <- fs[file.exists(fs)]

  # -- read all hits
  nCores <- gsParam$params$nCores
  setDTthreads(1)
  ofID1 <- ofID2 <- isOg <- inBuffer <- NULL
  hts <- rbindlist(lapply(fs, function(i){
    x <- fread(
      i, select = c("ofID1","ofID2","isOg","inBuffer"),
      na.strings = c("","NA"), showProgress = F)
    x <- subset(x, inBuffer & isOg & ofID1 != ofID2)
    return(x[,c("ofID1", "ofID2")])
  }))

  # -- summarize by arrays
  ofID1 <- ofID2 <- arrayID <- NULL
  gff[,arrv := as.character(as.numeric(as.factor(arrayID)))]
  av <- gff$arrv; names(av) <- gff$ofID
  hts[,`:=`(a1 = as.character(av[ofID1]),
            a2 = as.character(av[ofID2]))]

  # -- convert to syntenic orthogroups
  a1 <- a2 <- NULL
  ic <- with(hts, clus_igraph(
    id1 = c(a1, a2), id2 = c(a2, a1)))

  # -- add in genes with missing synOgs
  arrv <- NULL
  tmp <- subset(gff, !arrv %in% as.character(names(ic)))
  nc <- 1:uniqueN(tmp$arrv) + uniqueN(ic); names(nc) <- unique(tmp$arrv)
  ic <- c(ic, nc)
  gff[,`:=`(synOG = ic[as.character(arrv)], arrv = NULL)]
  return(gff)
}

#' @title combine inblk and syn OGs
#' @description
#' \code{combine_inblkSynOG} combine inblk and syn OGs
#' @rdname synteny
#' @export
combine_inblkSynOG <- function(genomeIDs,
                               gff,
                               gsParam){
  setDTthreads(1)
  if(gsParam$params$verbose)
    cat("Combining synteny-constrained and inblock orthogroups ...\n")

  if(gsParam$params$verbose)
    cat(sprintf("\tsyn OGs: %s, inblk OGs: %s",
                uniqueN(gff$synOG, na.rm = T),
                uniqueN(gff$inblkOG, na.rm = T)))
  noInBlk <- all(is.na(gff$inblkOG))

  synOG <- NULL
  if(noInBlk)
    gff[,inblkOG := synOG]

  ofID <- ofID1 <- ofID2 <- clus <- og <- combOG <- inblkOG <- synOG <- NULL
  inblk <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "inblkOG"]
  syn <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "synOG"]
  u <- with(inblk, paste(ofID1, ofID2))

  syn <- subset(syn, !paste(ofID1, ofID2) %in% u & ofID1 != ofID2)
  tmp <- rbind(syn[,c("ofID1", "ofID2")], inblk[,c("ofID1", "ofID2")])
  tmp[,clus := clus_igraph(ofID1, ofID2)]

  ov <- with(tmp, c(clus, clus))
  names(ov) <- with(tmp, c(ofID1, ofID2))

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
    gff[,inblkOG := NA]
  return(gff)
}

#' @title Split overlapping non-duplicated blocks
#' @description
#' \code{split_blks} function to deal with overlapping / interleafed blocks
#' @rdname synteny
#' @importFrom dbscan dbscan frNN
#' @export
split_blks <- function(hits, blkSize, verbose = F, maxIter = 5){
  setDTthreads(1)
  # -- function to count the overlapping hits
  count_ovlHits <- function(hits, maxDup){
    blkID <- i.blkID <- blk1 <- blk2 <- r1 <- r2 <- nOvlp <- totHits <- n1 <-
      n2 <- NULL
    bc <- calc_blkCoords(hits)
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
    blkOvl <- subset(out, !duplicated(out))
    blkOvl[,u := apply(.SD, 1, function(x) paste(x[order(x)], collapse = "_")),
           .SDcols = c("blk1", "blk2")]
    blkOvl <- subset(blkOvl, !duplicated(u))
    blkOvl[,u := NULL]
    tmp <- subset(hits, blkID %in% unique(c(blkOvl$blk1, blkOvl$blk2)))
    spl <- split(tmp, by = "blkID")
    blkOvl[,nOvlp := apply(blkOvl, 1, function(z){
      x <- spl[[as.character(z[[1]])]]
      y <- spl[[as.character(z[[2]])]]
      return(sum(x$ofID1 %in% y$ofID1, x$ofID2 %in% y$ofID2,
                 x$ofID1 %in% y$ofID1, x$ofID2 %in% y$ofID2))
    })]
    nv <- with(bc, nHits1 + nHits2); names(nv) <- bc$blkID
    blkOvl[,`:=`(n1 = nv[as.character(blk1)],
                 n2 = nv[as.character(blk2)])]
    blkOvl[,totHits := n1 + n2]

    setorder(blkOvl, -totHits)
    return(subset(blkOvl, nOvlp <= maxDup))
  }

  spl_ovlblks <- function(blk1, blk2, hits, blkSize){
    blkID <- isAnchor <- arrayOrd1 <- arrayOrd2 <- tmpb <- ord1 <- ord2 <- NULL
    tmp <- subset(hits, blkID %in% c(blk1, blk2) & isAnchor)
    tmp[,`:=`(ord1 = frank(arrayOrd1, ties.method = "dense"),
              ord2 = frank(arrayOrd2, ties.method = "dense"))]
    tmp[,tmpb := dbscan(frNN(cbind(ord1, ord2), eps = blkSize),
                        minPts = blkSize)$cluster]
    tmp <- subset(tmp, tmpb > 0)
    tmp[,`:=`(blkID = paste(tmpb, blkID), tmpb = NULL)]
    return(tmp)
  }

  isAnchor <- blk1 <- blk2 <- blkID <- NULL
  tmp <- subset(hits, isAnchor)
  runit <- T
  niter <- 1
  excl <- character()
  while(niter < maxIter & runit){
    if(verbose)
      cat(sprintf("%s: %s", niter, uniqueN(tmp$blkID)))
    niter <- niter + 1
    novl <- count_ovlHits(hits = tmp, maxDup = blkSize)
    novl <- subset(novl, !paste(blk1, blk2) %in% excl)
    novl <- subset(novl, !blk1 %in% novl$blk1[duplicated(novl$blk1)])
    novl <- subset(novl, !blk2 %in% novl$blk2[duplicated(novl$blk2)])
    novl <- subset(novl, !blk1 %in% blk2)
    novl <- subset(novl, !blk2 %in% blk1)
    noovl <- subset(tmp, !blkID %in% c(novl$blk1, novl$blk2))
    ovl <- subset(tmp, blkID %in% c(novl$blk1, novl$blk2))
    if(nrow(ovl) < blkSize){
      runit <- FALSE
    }else{
      ovls <- rbindlist(lapply(1:nrow(novl), function(k)
        spl_ovlblks(
          blk1 = novl$blk1[k], blk2 = novl$blk2[k],
          hits = ovl, blkSize = ceiling(blkSize/2))))
      tmp <- rbind(noovl, ovls)
      excl <- unique(c(excl, with(novl, unique(
        c(paste(blk1, blk2), paste(blk2, blk1))))))
    }
    if(verbose)
      cat(sprintf(" --> %s\n", uniqueN(tmp$blkID)))
  }
  tmp[,blkID := as.numeric(as.factor(blkID))]
  return(tmp)
}


#' @title pull hits within a buffer
#' @description
#' \code{find_hitsInBuff} pull hits within a buffer to any flagged with "isAnchor"
#' @rdname synteny
#' @importFrom dbscan dbscan frNN
#' @export
find_hitsInBuff <- function(hits, synBuff){
  setDTthreads(1)
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

#' @title pull hits within block coordinates
#' @description
#' \code{get_hitsInBlks} pull hits within block coordinates
#' @rdname synteny
#' @export
get_hitsInBlks <- function(blks, hits){
  setDTthreads(1)
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

