#' @title Find synteny between two genomes
#' @author JT Lovell
#' @description
#' \code{extract_synteny} Typically only used internally as part of the
#' syntenic_orthogroups pipeline.
#'
#' @name extract_synteny
#'
#' @param synParam file.path to the directory storing the input orthofinder
#' blast files and orthogroups.tsv. The orthogroups file can be in its
#' original subdirectory. Genesppace will only use the most recently modified
#' occurance of orthogroups.tsv in all subdirectories of blastDir.
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run) and ploidy
#' (integer vector of genome ploidies named by genomeIDs).
#' @param gsAnnot list of file locations for genespace annotation input files.
#' Alternatively, a named list containing element 'gff' which is a named vector
#' of the file paths to the gff-like genespace input file.
#'
#' @param blastDir file.path to the directory storing the input orthofinder
#' blast files and orthogroups.tsv. The orthogroups file can be in its
#' original subdirectory. Genesppace will only use the most recently modified
#' occurance of orthogroups.tsv in all subdirectories of blastDir.
#' @param path2mcscanx character string file.path pointing to the install
#' directory for the MCScanX program. In particular, there must be an
#' executable in this directory for MCScanX_h.
#' @param blastDir file.path to the directory storing the input orthofinder
#' blast files and orthogroups.tsv. The orthogroups file can be in its
#' original subdirectory. Genesppace will only use the most recently modified
#' occurance of orthogroups.tsv in all subdirectories of blastDir.
#' @param hits data.table containing annotated blast-format pairwise hits
#' @param allHits data.table with the full set of hits to search through
#' @param nHits1 integer of length 1 specifying the number of hits to retain for
#' each gene in genome1
#' @param nHits2 integer of length 1 specifying the number of hits to retain for
#' each gene in genome2
#' @param onlyOg logical, should only orthogroups be considered?
#' @param nGaps integer of length 1 specifying the -m 'gaps' MCScanX paramerter
#' @param blkSize integer of length 1 specifying the minimum size for a syntenic
#' block and the -s 'size' MCScanX parameter
#' @param synBuff integer of length 1 specifying the maximum euclidean distance
#' from an 'anchor' so that it can be considered syntenic
#' @param radius numeric of length 1 specifying the eps dbscan parameter; the
#' search radius within which to count clustered density-based xy points.
#' @param dropSmallBlks logical, should hits in blocks that are too small be
#' dropped from the dataset?
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#' @param rerank logical, should hits be reranked?
#' @param type character string added to the 'type' column of the hits object
#' @param gffFiles named vector of file paths pointing to the parsed gff-like
#' annotation files
#' @param genomeID1 character string specifying the first genome to consider
#' @param genomeID2 character string specifying the second genome to consider
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
#' @return a 'hits' data.table for a pairwise combination of genomes.
#' @seealso set_syntenyParams pairwise_syntenyUtils
#' @examples
#' \dontrun{
#' # coming soon
#' }

#' @title Extract synteny pipeline
#' @description
#' \code{extract_synteny} Main engine for synteny extraction in GENESPACE.
#' @rdname extract_synteny
#' @importFrom parallel mclapply
#' @importFrom grDevices pdf dev.off
#' @import data.table
#' @export
extract_synteny <- function(synParam,
                            gsParam,
                            gsAnnot,
                            nCores = NULL,
                            overwrite = FALSE){

  genomeIDs <- gsParam$genomeIDs[!gsParam$genomeIDs %in% gsParam$outgroup]
  runBlast <- blkID <- NULL
  on.exit(expr = setDTthreads(getDTthreads()))
  verbose <- gsParam$verbose
  if(is.null(nCores))
    nCores <- gsParam$nCores

  sp <- subset(synParam, runBlast)
  sp <- subset(sp, !genome1 %in% gsParam$outgroup & !genome2 %in% gsParam$outgroup)
  sf <- file.path(gsParam$synteny, with(sp, sprintf("%s_vs_%s_synHits.txt.gz", genome1, genome2)))
  if(all(file.exists(sf)) & !overwrite){
    warning("all synHits.txt.gz exist in /synteny\n\tSet overwrite = TRUE to run anyways.")
    sp[,synHitFiles := sf]
    return(sp)
  }else{
    if(verbose)
      cat("Pruning pairwise hits to syntenic regions:\n\tGenome Comp.\tnPrim. Hits / Secondary\n")

    # load full gff
    gff <- add_ofID2gff(
      read_gff(gsAnnot$gff),
      blastDir = gsParam$blast)

    # get vector of orthogroups
    ogv <- pull_orthogroups(
      path = gsParam$blast,
      genomeIDs = genomeIDs)$ogv

    setDTthreads(nCores)
    pdf(file.path(gsParam$results, "syntenicRegionHeatmap.pdf"), height = 8, width = 8)
    synHitFiles <- unlist(lapply(1:nrow(sp), function(i){
      g1 <- sp$genome1[i]; g2 <- sp$genome2[i]

      if(gsParam$verbose){
        s1 <- substr(g1, 1, 5)
        s2 <- substr(g2, 1, 5)
        if(nchar(s1) < 5)
          s1 <- sprintf("%s%s", s1, paste(rep(" ", 5-nchar(s1)), collapse = ""))
        if(nchar(s2) < 5)
          s2 <- sprintf("%s%s", s2, paste(rep(" ", 5-nchar(s2)), collapse = ""))
        cat(sprintf("\t%s vs %s:  ", s1, s2))
      }

      synHits <- parse_hits2synteny(
        blastDir = gsParam$blast,
        gffFiles = gsAnnot$gff,
        path2mcscanx = gsParam$path2mcscanx,
        synParam = synParam,
        genomeID1 = g1,
        genomeID2 = g2,
        gff = gff,
        nCores = nCores,
        ogv = ogv)

      tp <- data.table(synHits); tp[,blkID := NULL]
      pltdat <- plot_hits(
        hits = tp,
        gsParam = gsParam,
        gsAnnot = gsAnnot,
        reorderChrs = FALSE,
        minGenes = 50,
        nHeatmapBins = 2000)

      fo <- file.path(gsParam$synteny, sprintf("%s_vs_%s_synHits.txt.gz", g1, g2))
      fwrite(synHits, file = fo, sep = "\t")
      if(verbose)
        with(synHits,
             cat(sprintf("%s / %s\n",
                         sum(type == "primary"), sum(type != "primary"))))
      return(fo)
    }))
    dev.off()
    sp[,synHitFiles := synHitFiles]
    return(sp)
  }
}

#' @title Pull out all genes near self hits
#' @description
#' \code{pull_selfRegion} For an intra-genomic hits data.table only, a faster method to
#' extract just the region around the self hits.
#' @rdname extract_synteny
#' @import data.table
#' @export
pull_selfRegion <- function(hits, synBuff){
  genome2 <- blkID <- ord2 <- ofID1 <- ofID2 <- ord1 <- chr1 <- chr2 <- NULL

  buff <- sqrt(synBuff^2 * 2) + 1
  x <- subset(
    hits,
    chr1 == chr2 &
      abs(ord1 - ord2) <= buff)
  x[,blkID := sprintf("blk_%s_%s_%s", chr1, chr2, 1)]
  x[,`:=`(type = "primary", isAnchor = ofID1 == ofID2)]
  return(x)
}

#' @title Find hits within a buffer around anchors
#' @description
#' \code{pairwise_syntenyUtils} Use fixed-radius nearest-neighbor searches for
#' all hits within a buffer around anchors. Hits data.table must have a column
#' called 'isAnchor' specifying whether the hit is an anchor or one to be
#' assigned as inside or outside of the buffer.
#' @rdname extract_synteny
#' @importFrom dbscan frNN
#' @importFrom parallel mclapply
#' @import data.table
#' @export
split_hitsInBuffer <- function(hits,
                               synBuff,
                               nCores){
  inBuffer <- isAnchor <- ord2 <- ofID1 <- ofID2 <- ord1 <- chr1 <- chr2 <- NULL

  f_hib <- function(rb, x, y, ia){
    nn <- frNN(x = data.frame(x, y), eps = rb)
    return(1:length(x) %in% unique(unlist(nn$id[ia])))
  }

  anchorChrs <- with(subset(hits, isAnchor), unique(paste(chr1, chr2)))
  chrHits <- subset(hits, paste(chr1, chr2) %in% anchorChrs)
  spl <- split(chrHits, by = c("chr1", "chr2"))
  setDTthreads(1)
  chrHits <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){
    x[,inBuffer := f_hib(
      x = ord1,
      y = ord2,
      ia = isAnchor,
      rb = synBuff),
      by = c("chr1", "chr2")]
    return(x)
  }))
  setDTthreads(nCores)

  inBuff <- subset(chrHits, inBuffer)
  buffHits <- with(inBuff, paste(ofID1, ofID2))
  outBuff <- subset(hits, !paste(ofID1, ofID2) %in% buffHits)
  return(list(inBuffer = inBuff, outBuffer = outBuff))
}

#' @title Extract top-scoring hits
#' @description
#' \code{cull_hitsByScore} Subset a hits data.table based on the rank order
#' of the bit score ('score' column) for each unique gene ID ('ofID1/2' columns)
#' @rdname extract_synteny
#' @import data.table
#' @export
cull_hitsByScore <- function(hits, nHits1, nHits2, onlyOg){
  score <- scoreRank1 <- scoreRank2 <- ofID1 <- ofID2 <- ord1 <- og2 <- og1 <- NULL

  if(onlyOg){
    hits <- subset(hits, og1 == og2)
  }
  setorder(hits, ofID1, -score)
  hits[,scoreRank1 := 1:.N, by = "ofID1"]
  setorder(hits, ofID2, -score)
  hits[,scoreRank2 := 1:.N, by = "ofID2"]
  hits <- subset(hits, scoreRank1 <= nHits1 & scoreRank2 <= nHits2)
  return(hits)
}

#' @title Form dbscan clusters
#' @description
#' \code{cluster_dbscan} Build density-based 2d clusters based on the ord1 and
#' ord2 columns in a hits data.table.
#' @rdname extract_synteny
#' @importFrom dbscan frNN dbscan
#' @importFrom parallel mclapply
#' @import data.table
#' @export
cluster_dbscan <- function(hits,
                           nCores,
                           radius,
                           blkSize,
                           dropSmallBlks){
  n2 <- n1 <- blkID <- ord2 <- ofID1 <- ofID2 <- ord1 <- chr1 <- chr2 <- NULL

  x <- split(hits, by = c("chr1","chr2"))
  setDTthreads(1)
  x <- rbindlist(mclapply(x, mc.cores = nCores, function(y){
    y[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = radius),
                       minPts = blkSize)$cluster]
    return(y)
  }))
  setDTthreads(nCores)
  x <- subset(x, blkID > 0)
  x[,blkID := sprintf("blk_%s_%s_%s", chr1, chr2, blkID)]
  if(dropSmallBlks){
    x[, `:=`(n1 = length(unique(ofID1)), n2 = length(unique(ofID2))),
      by = "blkID"]
    x <- subset(x, n1 >= blkSize & n2 >= blkSize)
    x[,`:=`(n1 = NULL, n2 = NULL)]
  }
  return(x)
}

#' @title Pull buffered syntenic hits
#' @description
#' \code{pull_blkBuff} Given a hits object and synteny parameters, extract
#' synteny regions by culling to syntenic 'anchors' then pull all hits within
#' a buffer around the anchors.
#' @rdname extract_synteny
#' @import data.table
#' @export
pull_blkBuff <- function(hits,
                         allHits,
                         blkSize,
                         nGaps,
                         path2mcscanx,
                         synBuff,
                         nCores,
                         rerank){
  isAnchor <- u <- ord2 <- nAnc <- ord1 <- NULL

  if(rerank){
    hits[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
               ord2 = frank(ord2, ties.method = "dense"))]
  }
  mcs <- run_mcscanx(
    hits = hits,
    blkSize = blkSize,
    nGaps = nGaps,
    path2mcscanx = path2mcscanx)
  allHits[,isAnchor := u %in% names(mcs[!is.na(mcs)])]
  allHits[,nAnc := sum(isAnchor), by = c("chr1","chr2")]
  allHits <- subset(allHits, nAnc >= blkSize)
  inBuff <- split_hitsInBuffer(
    hits = allHits,
    nCores = nCores,
    synBuff = synBuff)$inBuffer
  out <- cluster_dbscan(
    hits = inBuff,
    radius = synBuff,
    nCores = nCores,
    blkSize = blkSize,
    dropSmallBlks = T)
  out[, `:=` (nAnc = NULL, inBuffer = NULL)]
  return(out)
}

#' @title Extract synteny homeologous regions
#' @description
#' \code{pull_selfHomeologSynteny} Given ploidy for a genome and a self-hit
#' region masked hits data table, pull out syntenic regions.
#' @rdname extract_synteny
#' @import data.table
#' @export
pull_buffSynteny <- function(nHits1,
                             nHits2,
                             hits,
                             onlyOg,
                             blkSize,
                             nGaps,
                             nCores,
                             path2mcscanx,
                             synBuff,
                             type){

  cullHits <- cull_hitsByScore(
    hits = hits,
    nHits1 = nHits1,
    nHits2 = nHits2,
    onlyOg = onlyOg)
  synHits <-  pull_blkBuff(
    allHits = hits,
    hits = cullHits,
    blkSize = blkSize,
    nGaps = nGaps,
    nCores = nCores,
    path2mcscanx = path2mcscanx,
    synBuff = synBuff,
    rerank = TRUE)
  synHits[,type := type]
  return(synHits)
}

#' @title Extract synteny homeologous regions
#' @description
#' \code{pull_selfHomeologSynteny} Given ploidy for a genome and a self-hit
#' region masked hits data table, pull out syntenic regions.
#' @rdname extract_synteny
#' @import data.table
#' @importFrom parallel mclapply
#' @export
cull_hitsWithinAnchorBounds <- function(hits, nCores){
  isAnchor <- ord2 <- ord1 <- NULL

  spl <- split(hits, by = "blkID")
  setDTthreads(1)
  hits <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){
    rng1 <- with(subset(x, isAnchor), range(ord1))
    rng2 <- with(subset(x, isAnchor), range(ord2))
    x <- subset(x, ord1 >= rng1[1] & ord2 >= rng2[1] & ord1 <= rng1[2] & ord2 <= rng2[2])
    return(x)
  }))
  setDTthreads(nCores)
  return(hits)
}

#' @title Run MCScanX from R
#' @description
#' \code{run_mcscanx} Calculate MCScanX 'collinearity' from an annotated
#' blast-formatted data.table
#' @rdname extract_synteny
#' @import data.table
#' @export
run_mcscanx <- function(hits,
                        blkSize,
                        nGaps,
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
    getwd(),
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

#' @title Engine for pairwise synteny
#' @description
#' \code{parse_hits2synteny} Find synteny between two genomes and pull all hits
#' within specified buffer
#' @rdname extract_synteny
#' @import data.table
#' @export
parse_hits2synteny <- function(blastDir,
                               gffFiles,
                               path2mcscanx,
                               genomeID1,
                               genomeID2,
                               synParam,
                               gff,
                               ogv,
                               nCores){
  setDTthreads(nCores)
  blkID <- type <- ofID1 <- ofID2 <- u <- nu1 <- nu2 <- NULL

  ##############################################################################
  # 0. Check parameters
  sp <- subset(synParam, genome1 == genomeID1 & genome2 == genomeID2)[1,]
  genome1 <- genomeID1
  genome2 <- genomeID2
  ploidy1  <- sp$ploidy1
  ploidy2 <- sp$ploidy2
  onlyOg <- sp$onlyOg
  onlyOgHomeo <- sp$onlyOgHomeo
  onlyOgSecond <- sp$onlyOgSecond
  synBuff <- sp$synBuff
  synBuffHomeo <- sp$synBuffHomeo
  synBuffSecond <- sp$synBuffSecond
  blkSize <- sp$blkSize
  blkSizeHomeo <- sp$blkSizeHomeo
  blkSizeSecond <- sp$blkSizeSecond
  nGaps <- sp$nGaps
  nGapsHomeo <- sp$nGapsHomeo
  nGapsSecond <- sp$nGapsSecond
  nSecondHits1 <- sp$nSecondHits1
  nSecondHits2 <- sp$nSecondHits2

  ploidy <- max(c(ploidy1, ploidy2))
  n2hits <- max(c(nSecondHits1, nSecondHits2))
  if(n2hits > 0){
    if(nSecondHits1 == 0)
      nSecondHits1 <- n2hits
    if(nSecondHits2 == 0)
      nSecondHits2 <- n2hits
  }

  # load blast
  a <- read_blast(
    path = blastDir,
    ofID1 = read_orthofinderSpeciesIDs(blastDir)[genome1],
    ofID2 = read_orthofinderSpeciesIDs(blastDir)[genome2],
    onlyIDScore = TRUE)


  # add in missing orthogroups
  ov <- gff$ord;  gv <- gff$genome; cv <- gff$chr
  names(cv) <- names(gv) <- names(ov) <- gff$ofID
  uid <- unique(c(a$ofID1, a$ofID2))
  idNoOg <- uid[!uid %in% names(ogv)]
  ogz <- paste0("noOg_", 1:length(idNoOg)); names(ogz) <- idNoOg
  ogv <- c(ogv, ogz)

  # add genome positions to blast
  a[,`:=`(chr1 = cv[ofID1], chr2 = cv[ofID2],
          ord1 = ov[ofID1], ord2 = ov[ofID2],
          og1 = ogv[ofID1], og2 = ogv[ofID2],
          u = paste(ofID1, ofID2))]

  # drop chromosome combinations smaller than the block size
  a[,`:=`(nu1 = length(unique(ofID1)), nu2 = length(unique(ofID2))),
    by = c("chr1","chr2")]
  a <- subset(a, nu1 >= blkSize & nu2 >= blkSize)
  a[,`:=`(nu1 = NULL, nu2 = NULL)]

  ##############################################################################
  # 2. Extract synteny
  # Pull self synteny in intragenomic hits
  homeoOut <- NULL; secondOut <- NULL
  if(genome1 == genome2)
    primaryOut <- pull_selfRegion(
      hits = a,
      synBuff = synBuff)

  # Pull primary synteny in intergenomic hits
  if(genome1 != genome2)
    primaryOut <- pull_buffSynteny(
      nHits1 = ploidy2,
      nHits2 = ploidy1,
      onlyOg = onlyOg,
      hits = a,
      nCores = nCores,
      blkSize = blkSize,
      nGaps = nGaps,
      path2mcscanx = path2mcscanx,
      synBuff = synBuff,
      type = "primary")

  # Pull homeologous synteny in intragenomic polyploid hits
  if(ploidy > 1 & genome1 == genome2)
    homeoOut <- pull_buffSynteny(
      nHits1 = ploidy - 1,
      nHits2 = ploidy - 1,
      nCores = nCores,
      hits = subset(a, !u %in% primaryOut$u),
      onlyOg = onlyOgHomeo,
      blkSize = blkSizeHomeo,
      nGaps = nGapsHomeo,
      path2mcscanx = path2mcscanx,
      synBuff = synBuffHomeo,
      type = "homeolog")

  # pull 'secondary' hits after masking primary/homeologs
  if(n2hits > 0)
    secondOut <- pull_buffSynteny(
      nHits1 = nSecondHits1,
      nHits2 = nSecondHits2,
      onlyOg = onlyOgSecond,
      nCores = nCores,
      hits = subset(a, !u %in% unique(c(primaryOut$u, homeoOut$u))),
      blkSize = blkSizeSecond,
      nGaps = nGapsSecond,
      path2mcscanx = path2mcscanx,
      synBuff = synBuffSecond,
      type = "secondary")

  ############################################################################
  # 3. Combine results
  hits <- rbind(primaryOut, homeoOut, secondOut, fill = T)
  hits[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2],
             chr1 = cv[ofID1], chr2 = cv[ofID2],
             blkID = paste(blkID, type),
             genome1 = gv[ofID1], genome2 = gv[ofID2])]
  out <- cull_hitsWithinAnchorBounds(
    hits = hits, nCores = nCores)

  return(out)
}


