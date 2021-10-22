#' @title annotate_gff
#' @description
#' \code{annotate_gff} annotate_gff
#'
#' @name annotate_gff
#'
#' @param gsParam a list containing all parameters for a GENESPACE run. See
#' init_genespace
#' @param genomeIDs an optional vector of genomeIDs to consider. If not
#' specified (default) taken from gsParam$genomeIDs$genomeIDs
#' @param genomeID single value fromo genomeIDs
#' @param synBuff synteny buffer, see set_synParam
#' @param overwrite logical, should existing results be overwritten?
#' @param minGenes4of integer specifying the minimum number of genes allowed
#' for an orthofinder run
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run

#' @details ...
#'
#' @return ...
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#' @title add number of anchors to gff
#' @description
#' \code{annotate_gff} add number of anchors to gff
#' @rdname annotate_gff
#' @importFrom parallel mclapply
#' @export
annotate_gff <- function(gsParam,
                         genomeIDs = NULL,
                         overwrite = FALSE){

  genome <- ord  <- arrayID <- isArrayRep <- med <- pepLen <- rnk <- ofID <- NULL
  medbp <- start <- dist2med <- dist2bp <- og <- n <- chr <- globOG <- NULL

  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(file.exists(gffFile) & !overwrite){
    tmp <- fread(gffFile, na.strings = c("-", "NA", ""), showProgress = F)
    if(all(c("pepLen", "globOG", "arrayID","isArrayRep") %in% colnames(tmp)))
      stop("annotated gff file exists and !overwrite, so not running ...\n")
  }
  if(!is.data.table(gsParam$params$synteny))
    stop("Must run set_syntenyParams first!\n")

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs

  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  synBuff <- max(gsParam$params$synteny$synBuff)
  if(verbose)
    cat("Loading annotations ...\n")
  # -- get paths to the orthofinder run
  if(is.na(gsParam$paths$orthogroupsDir)){
    if(verbose)
      cat("\tIndexing location of orthofinder results ... ")
    gsParam <- find_orthofinderResults(gsParam)
    if(verbose)
      cat("Done!\n")
  }

  ##############################################################################
  # 1. Load the gff and add global metadata
  ##############################################################################
  # -- read in the gff
  if(verbose)
    cat("\tReading the gffs ... ")
  gff <- read_gff(gsParam$paths$gff)
  gff <- add_ofID2gff(gff, gsParam$paths$blastDir)
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)

  # -- add peptide length
  if(verbose)
    cat("Done!\n\tPulling gene lengths ... ")
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)

  # -- add global orthogroupsto the gff
  if(verbose)
    cat("Done!\n\tParsing global orthogroups ... ")
  ogs <- parse_ogs(gsParam)
  gff <- merge(gff, ogs, by = c("genome","id"), all.x = T)
  gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
  setnames(gff, "ogID", "globOG")

  ##############################################################################
  # 2. Define collinear arrays
  ##############################################################################
  # -- build arrays from method in gsParam
  if(verbose)
    cat("Done!\nDefining collinear orthogroup arrays ... \n")
  gff <- add_arrays2gff(gsParam = gsParam, gff = gff)
  gff[,arrayID := paste(arrayID, og)]
  gff[,n := .N, by = "arrayID"]
  gff$arrayID[gff$n == 1] <- NA
  gff[,`:=`(n = NULL, isArrayRep = TRUE, og = NULL)]

  # -- choose the array reps
  if(verbose)
    cat("\tChoosing array representative genes ... ")
  gffa <- subset(gff, !is.na(arrayID))
  setkey(gffa, genome, ord)
  gffn <- subset(gff, is.na(arrayID))
  gffa[,arrayID := sprintf(
    "%s_%s_%s",
    genome, chr, as.numeric(factor(arrayID, levels = unique(arrayID))))]
  gffa[,med := as.numeric(median(ord)), by = "arrayID"]
  gffa[,medbp := as.numeric(median(start)), by = "arrayID"]
  gffa[,`:=`(dist2med = abs(med - ord),
             dist2bp = abs(medbp - start))]
  setorder(gffa, arrayID, dist2med, dist2bp, -pepLen)
  gffa[,rnk := 1:.N, by = "arrayID"]
  gffa[,`:=`(isArrayRep = rnk == 1, rnk = NULL, dist2med = NULL, dist2bp = NULL,
             med = NULL, medbp = NULL)]
  gff <- rbind(gffa, gffn)
  setkey(gff, genome, ord)

  if(verbose)
    cat(sprintf("Done!\nWriting gff to file: %s", gffFile))
  gff[,`:=`(synOG = NA, inBlkOG = NA, combOG = NA, og = globOG, refCoord = NA)]
  fwrite(gff, file = gffFile, sep = "\t", quote = F, showProgress = F)
  return(gsParam)
}

#' @title add array representative to a gff object
#' @description
#' \code{add_arrayRep2gff} choose most central gene by orthogroup
#' @rdname annotate_gff
#' @importFrom parallel mclapply
#' @export
add_arrayRep2gff <- function(gff,
                             gsParam){
  synArr <- chr <- NULL
  # -- count peptides
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)
  genomeIDs <- unique(gff$genome)
  nCores <- gsParam$params$nCores

  # -- count number of orthologs
  di <- dir.exists(gsParam$paths$orthologuesDir)
  dl <- length(list.files(gsParam$paths$orthologuesDir)) > 1

  nGenome <- nGenes <- gen2 <- id1 <- gen1 <- genome <- id <- og <- NULL
  if(di && dl){
    ogcnt <- rbindlist(mclapply(genomeIDs, mc.cores = nCores, function(i){
      ogs <- parse_orthologues(gsParam = gsParam, refGenome = i)
      ogn <- ogs[,list(nGenome = uniqueN(gen2),
                       nGenes = .N), by = c("gen1","id1")]
      return(ogn)
    }))
    setorder(ogcnt, -nGenome, -nGenes)
    ogcnt <- subset(ogcnt, !duplicated(paste(gen1, id1)))
    nog <- ogcnt$nGenome; ng <- ogcnt$nGenes
    names(ng) <- names(nog) <- with(ogcnt, paste(gen1, id1))
    gff[,`:=`(nGenomeOrthologs = nog[paste(genome, id)],
              nTotalOrthologs = ng[paste(genome, id)])]
    gff$nGenomeOrthologs[is.na(gff$nGenomeOrthologs)] <- 0
    gff$nTotalOrthologs[is.na(gff$nTotalOrthologs)] <- 0
  }else{
    gff[,`:=`(nGenomeOrthologs = 0,
              nTotalOrthologs = 0)]
  }

  # -- split into single and multiple member arrays
  gffi <- data.table(gff)
  d2h <- gsParam$params$maxDistBtwPgHits
  gff[,synArr := as.integer(as.factor(paste(genome, chr, og)))]
  gff[,nog := .N, by = "synArr"]
  g1 <- subset(gff, nog == 1)
  g1[,nog := NULL]
  g2 <- subset(gff, nog > 1)

  # -- calculate the maximum distance between genes in an array
  maxJump <- ord <- NULL
  g2[,maxJump := max(diff(ord[order(ord)])), by = "synArr"]
  g2r <- subset(g2, maxJump > d2h)
  g2 <- subset(g2, maxJump <= d2h)

  # -- cluster genes in arrays with big jumps
  clus <- NULL
  if(nrow(g2r) > 1){
    g2[,maxJump := NULL]
    g2r[,clus := dbscan(frNN(cbind(ord, ord), eps = d2h), minPts = 0)$cluster,
        by = "synArr"]
    g2r[,synArr := paste(synArr, clus)]
    g2 <- rbind(g2,  g2r[,colnames(g2),with = F])
  }

  # -- calculate distance to the median
  dist2median <- nGenomeOrthologs <- pepLen <- nTotalOrthologs <- ofID <- NULL
  g2[,dist2median := abs(as.numeric(median(ord, na.rm = T)) - ord),
     by = c("synArr","genome","chr")]

  # -- order and rank genes, choosing representatives for each array
  setorder(g2, genome, chr, synArr, -nGenomeOrthologs, -nTotalOrthologs,
           dist2median, -pepLen, ord)
  arep <- rbind(g1, g2[,colnames(g1), with = F])
  sar <- as.numeric(as.factor(arep$synArr)); names(sar) <- arep$ofID
  arep <- subset(arep, !duplicated(synArr))
  gffi[,`:=`(synArray = sar[ofID],
             isArrayRep = ofID %in% arep$ofID)]
  return(gffi)
}


#' @title add_arrays2gff
#' @description
#' \code{add_arrays2gff} add_arrays2gff
#' @rdname annotate_gff
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
add_arrays2gff <- function(gsParam,
                           gff){

  ##############################################################################
  ##############################################################################
  arrayID <-  genome <-  chr <-  globOG <-  n <-  rng <- ord  <- clus  <- og  <-  collinearOG <-
  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  synBuff <- max(gsParam$params$synteny$synBuff)

  # -- make global arrays from orthogroups
  gff[,arrayID := sprintf("%s_%s_%s", genome, chr, globOG)]
  gff[,n := .N, by = "arrayID"]

  # -- combine 1x ogs with ogs in regions < synBuff
  g1 <- subset(gff, n == 1)
  g2 <- subset(gff, n > 1)
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

  # -- NA out arrays with just one member
  gff[,n := .N, by = "arrayID"]
  gff[,og := arrayID]
  gff$arrayID[gff$n == 1] <- NA
  gff[,n := NULL]
  setkey(gff, genome, ord)

  # -- print updates and number of global orthogroups
  if(verbose){
    cat("\tUsing collinear orthogroups for array identity:\n")
    nu <- lapply(split(subset(gff, !is.na(arrayID)), by = "genome"), function(x)
      cat(sprintf("\t%s%s: %s genes in %s collinear arrays\n",
                  app, x$genome[1], nrow(x), uniqueN(x$arrayID))))
  }

  return(gff)
}
