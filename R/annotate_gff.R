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
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param blastDir file.path to the location of the blast results.

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
                         minGenes4of = 40,
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
  fwrite(gff, file = gffFile, sep = "\t")
  return(gsParam)
}

#' @title add peptide length to a gff object
#' @description
#' \code{add_pepLen2gff} read the peptide lengths and merge with the gff
#' @rdname annotate_gff
#' @export
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

#' @title add orthofinder ID to a gff object
#' @description
#' \code{add_ofID2gff} read the orthofinder species and gene IDs and merge
#' these with the gff-like data.table
#' @rdname annotate_gff
#' @export
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

#' @title Call arrays from within-chromosome orthogroups
#' @description
#' \code{recall_arrays} More precise than global runs, since global orthogroups
#' can expand as a function of paralogy among genomes.
#' @rdname annotate_gff
#' @import data.table
#' @export
recall_arrays <- function(gsParam, minGenes4of = 50){
  genome1 <- genome2 <- genome <- synArray <- recallArray <- ofID <- NULL
  if(!check_orthofinderInstall("orthofinder"))
    stop("orthofinder not found. Open R from an environment with orthofinder in the path\n")
  if(gsParam$params$verbose)
    cat("Re-calling syntenic arrays from within-chromosome orthofinder runs ...\n\t(genome: n1x genes/n2+x / narrays for global .. rerun)\n")

  # -- read in the gff
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, showProgress = F, na.strings = c("", "NA"))

  # -- get the synteny buffer info
  tmp <- subset(gsParam$params$synteny, genome1 %in% gsParam$genomes$genomeIDs & genome1 == genome2)
  sb <- tmp$synBuff; names(sb) <- tmp$genome1; tmp <- NULL

  # -- for each genome, ...
  arrayDt <- rbindlist(lapply(gsParam$genomes$genomeIDs, function(i){

    # -- pull global array info and print (if verbose)
    if(gsParam$params$verbose){
      tmpg <- subset(gff, genome == i)
      hasArr <- tmpg$synArray[duplicated(tmpg$synArray)]
      noar <- subset(tmpg, !synArray %in% hasArr)
      hasArr <- subset(tmpg, synArray %in% hasArr)
      cat(sprintf(
        "\t%s: %s / %s / %s .. ",
        i, sum(!duplicated(tmpg$synArray)), nrow(hasArr),
        uniqueN(hasArr$synArray)))
    }

    arrout <- clus_arrays(
      gff = gff,
      gsParam = gsParam,
      synBuff = sb[i],
      minGenes4of = 50,
      genomeID = i)

    # -- parse and merge
    tmp <- merge(gff, arrout[,c("ofID", "array")], by = "ofID", all.y = T)
    tmp[,recallArray := as.numeric(as.factor(paste(genome, array)))]
    tmp$recallArray[is.na(tmp$recallArray)] <- paste("n", tmp$synArray[is.na(tmp$recallArray)])
    tmp[,`:=`(synArray = recallArray, recallArray = NULL)]

    # -- pull global array info and print (if verbose)
    if(gsParam$params$verbose){
      tmpg <- data.table(tmp)
      hasArr <- tmpg$synArray[duplicated(tmpg$synArray)]
      noar <- subset(tmpg, !synArray %in% hasArr)
      hasArr <- subset(tmpg, synArray %in% hasArr)
      cat(sprintf(
        "%s / %s / %s\n",
        sum(!duplicated(tmpg$synArray)), nrow(hasArr),
        uniqueN(hasArr$synArray)))
    }

    return(tmp)
  }))
  if(gsParam$params$verbose)
    cat("\tChoosing representatives ... ")

  ra <- arrayDt$synArray; names(ra) <- arrayDt$ofID
  gff[,tmp := ra[ofID]]
  gff$tmp[is.na(gff$tmp)] <- gff$synArray[is.na(gff$tmp)]
  gff[,`:=`(synArray = tmp, tmp = NULL)]
  gff2 <- add_arrayRep2gff(gff = gff, gsParam = gsParam)
  if(gsParam$params$verbose)
    cat("Done!\n")
  return(gff2)
}

#' @title Cluster arrays within a single genome
#' @description
#' \code{clus_arrays} Routine piped by recall_arrays
#' @rdname annotate_gff
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom Biostrings readAAStringSet
#' @export
clus_arrays <- function(gff,
                        gsParam,
                        synBuff,
                        minGenes4of = 40,
                        genomeID){

  ofID1 <- ofID2 <- chr1 <- chr2 <- n <- ord1 <- ord2 <- ofID <- maxJump <- NULL
  ord <- clus <- chr <- og <- NULL
  # -- read the gff and get vectors of chrs and orders
  chrv <- gff$chr; ordv <- gff$ord; names(chrv) <- names(ordv) <- gff$ofID

  # -- read the peptide file
  if(is.na(gsParam$paths$blastDir))
    gsParam <- find_orthofinderResults(gsParam)
  ofsp <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  pepSS <- readAAStringSet(file.path(gsParam$paths$blastDir,
                                     sprintf("Species%s.fa", ofsp[genomeID])))

  # -- read the blast file
  bl <- read_blast(
    path = gsParam$paths$blastDir,
    ofID1 = ofsp[genomeID],
    ofID2 = ofsp[genomeID],
    onlyIDScore = F)
  setnames(bl, 1:2, c("ofID1", "ofID2"))
  bl[,`:=`(chr1 = chrv[ofID1], ord1 = ordv[ofID1],
           chr2 = chrv[ofID2], ord2 = ordv[ofID2])]

  # -- subset the self hits to those in proximity and split by chr
  bl <- subset(bl, chr1 == chr2)
  bl[,n := uniqueN(c(ofID1, ofID2)), by = "chr1"]
  bl <- subset(bl, n >= minGenes4of)
  bl <- subset(bl,  abs(ord1 - ord2) <= synBuff*2)[,1:13]
  spl <- lapply(split(bl, by = "chr1"), function(x) x[,1:12])
  spl <- spl[order(-sapply(spl, nrow))]

  selfOgs <- rbindlist(mclapply(names(spl), mc.cores = gsParam$params$nCores, mc.preschedule = F, function(i){
    tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4arrayTMPdir", i))
    if(dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    x <- data.table(spl[[i]])
    og <- run_ofFromObj(
      blast00 = x, blast01 = x, blast10 = x, blast11 = x,
      pep0 = pepSS, pep1 = pepSS, writeDir = tmpDir)
    if(dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    return(og)
  }))

  selfOgs <- subset(selfOgs, !duplicated(ofID))
  selfOgs[,`:=`(chr = chrv[ofID], ord = ordv[ofID], maxJump = 1, clus = 1)]
  selfOgs[,n := uniqueN(ofID), by = c("chr","og")]
  so1 <- subset(selfOgs, n == 1)
  so2 <- subset(selfOgs, n > 1)
  so2[,maxJump := max(diff(ord)), by = c("chr", "og")]
  selfOgs <- rbind(so1, subset(so2, maxJump <= synBuff * 2))

  g2r <- subset(selfOgs, maxJump > synBuff * 2)
  # -- cluster and add back in those with broad hits
  if(nrow(g2r) > 1){
    g2r[,clus := dbscan(frNN(cbind(ord, ord), eps = synBuff * 2), minPts = 0)$cluster,
        by = c("synOg", "chr")]
    selfOgs <- rbind(selfOgs, g2r)
  }
  selfOgs[,array := paste(chr, clus, as.numeric(as.factor(og)), sep = "_"),
          by = "chr"]
  return(selfOgs)
}

#' @title add_arrays2gff
#' @description
#' \code{add_arrays2gff} add_arrays2gff
#' @rdname annotate_gff
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
add_arrays2gff <- function(gsParam,
                           gff,
                           minGenes4of = 40){

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
    app <- ifelse(gsParam$params$recallArrays, "###\t", "")
    if(gsParam$params$recallArrays){
      cat("\tCollinear orthogroups for array identity (ignoring these):\n")
    }else{
      cat("\tUsing collinear orthogroups for array identity:\n")
    }
    nu <- lapply(split(subset(gff, !is.na(arrayID)), by = "genome"), function(x)
      cat(sprintf("\t%s%s: %s genes in %s collinear arrays\n",
                  app, x$genome[1], nrow(x), uniqueN(x$arrayID))))
  }

  # -- re-call arrays within chrs if necessary
  if(gsParam$params$recallArrays){
    if(verbose)
      cat("\tRe-calling arrays within chromosomes:\n")
    gffArrays <- pull_ogsByChr(
      gsParam = gsParam,
      gff = gff,
      minGenes4of = minGenes4of,
      synBuff = synBuff,
      nCores = nCores)
    gff <- merge(gff, gffArrays, all.x = T, by = "ofID")
    gff$collinearOG[is.na(gff$collinearOG)] <- paste0("NOG",1:sum(is.na(gff$collinearOG)))
    gff[,arrayID := sprintf("%s_%s_%s", genome, chr, collinearOG)]
    gff[,n := .N, by = "arrayID"]
    gff$arrayID[gff$n == 1] <- NA
    gff[,n := .N, by = c("genome", "chr")]

    # -- use global orthogroups for chrs that are too small
    if(verbose)
      cat(sprintf(
        "\tUsing global orthogroups for %s genes on chrs with < %s loci\n",
        sum(gff$n < 40), minGenes4of))
    gff$collinearOG[gff$n < 40] <- gff$globOG[gff$n < 40]

    # -- rename and return
    gff[,arrayID := sprintf("%s_%s_%s", genome, chr, collinearOG)]
    gff[,`:=`(n = NULL, collinearOG = NULL)]
  }
  return(gff)
}

#' @title pull_ogsByChr
#' @description
#' \code{pull_ogsByChr} pull_ogsByChr
#' @rdname annotate_gff
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
pull_ogsByChr <- function(gsParam,
                          gff,
                          synBuff,
                          minGenes4of,
                          nCores){

  V1 <- V2 <- u1 <- u2 <- o1 <- o2 <- nGenes <- NULL
  verbose <- gsParam$params$verbose
  ofsp <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  cv <- with(gff, paste(genome, chr)); ov <- gff$ord
  names(ov) <- names(cv) <- gff$ofID

  # -- load the blast files and parse to collinear genes
  bl <- mclapply(gsParam$genomes$genomeIDs, mc.cores = nCores, function(i){
    tmp <- read_blast(
      path = gsParam$paths$blastDir,
      ofID1 = ofsp[i],
      ofID2 = ofsp[i],
      onlyIDScore = F)
    tmp[,`:=`(genome = i, u1 = cv[V1], u2 = cv[V2],
              o1 = ov[V1], o2 = ov[V2])]
    tmp <- subset(subset(
      tmp, u1 == u2),
      abs(o1 - o2) <= synBuff)
    tmp <- split(tmp, by = c("u1","u2"))
    tmp <- lapply(tmp, function(x) x[,1:12])
    tmp <- tmp[sapply(tmp, nrow) >= minGenes4of]
    return(tmp)
  })
  names(bl) <- gsParam$genomes$genomeIDs

  # -- load the peptides
  pepspl <- sapply(gsParam$genomes$genomeIDs, USE.NAMES = T, simplify = F, function(i)
    readAAStringSet(file.path(gsParam$paths$blastDir,
                              sprintf("Species%s.fa", ofsp[i]))))

  # -- run orthofinder for each chr
  ogs <- rbindlist(lapply(names(bl), function(i){
    if(verbose)
      cat(sprintf("\t\t%s: ", i))
    bo <- names(bl[[i]])[order(-sapply(bl[[i]], nrow))]
    out <- rbindlist(mclapply(bo, mc.cores = nCores, mc.preschedule = F, function(j){
      x <- bl[[i]][[j]]
      tmpDir <- file.path(gsParam$params$wd, sprintf(
        "%s_%s_og4inBlkTMPdir", i, x[[1]][1], x[[2]][1]))
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      tmp <- run_ofFromObj(
        blast00 = x,
        blast01 = x,
        blast10 = x,
        blast11 = x,
        pep0 = pepspl[[i]],
        pep1 = pepspl[[i]],
        writeDir = tmpDir)
      unlink(tmpDir, recursive = T)
      tmp[,`:=`(genome = i, chr = j)]
      return(subset(tmp, !duplicated(tmp)))
    }))
    out[,nGenes := .N, by = c("og","genome", "chr")]
    tmp <- subset(out, nGenes > 1)
    if(verbose)
      cat(sprintf("%s genes in %s collinear arrays\n",
                  nrow(tmp), uniqueN(paste(tmp$chr, tmp$og))))
    setnames(out, "og", "collinearOG")
    return(out[,c("ofID","collinearOG")])
  }))
  return(ogs)
}
