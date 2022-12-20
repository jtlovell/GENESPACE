#' @title Build genespace pangenome
#'
#' @description
#' \code{pangenome} Convert orthogroup and synteny information into a
#' pangenome database. Predict locations of orthogroups that are missing a
#' node in the reference.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param refGenome character string matching one of the genomeIDs in gsParam
#' @param genomeIDs character vector, specifying which genomes to use. Defaults
#' to all genomeIDs specification in gsParam.
#' @param propAssignThresh numeric of length 1, the minimum proportion of genes
#' in a pangenome entry needed to anchor the physical position to a new location
#' @param verbose logical, should updates be printed to the console?
#'
#' @details The pangenome annotation is a projection of syntenic orthogroups
#' on the physical coordinate system of a reference genome. The pangenome
#' function runs the following pipeline.
#' \enumerate{
#' \item within-block orthogroups and synteny-constrained global orthogroups
#' are merged.
#' \item the reference position is projected for all genes within syntenic
#' block bounds.
#' \item physical projected positions for each syntenic orthgroup is clustered,
#' permitting multiple reference locations (e.g. in a polyploid)
#' \item orthogroups with placements on the reference are populated and un-
#' placed orthgroups are added within reference position information
#' \item array members are added to the orthogroups
#' \item non-syntenic orthologs are added and flagged
#' }
#'
#' @return a data.table with lists of gene ids for each genome. Each row
#' corresponds to a unique combination of orthogroup and physical position
#' on the reference. Where multiple positions are inferred on a single
#' chromosome, the positions are broken out by the column 'clus'.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#
#' @import data.table
#' @importFrom stats complete.cases median
#' @export
pangenome <- function(gsParam,
                      genomeIDs = NULL,
                      refGenome = NULL,
                      propAssignThresh = .25,
                      maxNonSynPlaces = 10,
                      verbose = T){

  ##############################################################################
  # -- ad hoc function to cluster interpolated positions into a pangenome
  clus_interp2pg <- function(bed, interp, minScaleProp, synBuff){

    og <- interpChr <- interpOrd <- isAnchor <- jumpLeft <- clus <- nClus <-
      nOgMem <- ofID <- nClusMem <- propInClus <- scaleProp <- pgOrd <-
      pgID <- NULL

    bd <- bed[,c("ofID", "og")]
    bd <- subset(bd, !duplicated(bd))
    intp <- interp[,c("ofID", "interpGenome",
                      "interpChr", "interpOrd", "isAnchor")]
    intp <- subset(intp, !duplicated(intp))

    tmp <- merge(bd, intp, by = "ofID", allow.cartesian = TRUE)

    # -- cluster nearby hits
    setkey(tmp, og, interpChr, interpOrd, isAnchor)
    tmp[,jumpLeft := c(synBuff + 1, diff(interpOrd)), by = c("interpChr", "og")]
    tmp[,clus := as.integer(jumpLeft > synBuff), by = c("interpChr")]
    tmp[,clus := cumsum(clus), by = c("interpChr", "og")]
    tmp[,clus := paste(interpChr, clus)]
    tmp[,nClus := uniqueN(clus), by = "og"]

    # -- calculate the proportions in each placement, dropping those with < 25%
    tmp[,nOgMem := uniqueN(ofID), by = "og"]
    tmp[,nClusMem := uniqueN(ofID), by = c("og", "clus")]
    tmp[,propInClus := nClusMem/nOgMem]
    tmp[,scaleProp := scale_between(
      x = propInClus, min = 0, max = 1, scale1toMean = FALSE),
      by = "og"]
    mNopos <- subset(tmp, scaleProp < minScaleProp)
    tmp <- subset(tmp, scaleProp >= minScaleProp)

    # --- get median positions for each orthogroup
    pg <- tmp[,list(pgOrd = median(interpOrd)),
              by = c("og", "interpChr", "clus")]
    setkey(pg, pgOrd)
    pg[,pgID := sprintf("pg_%s", 1:.N)]
    setnames(pg, "interpChr", "pgChr")
    pg[,clus := NULL]
    setcolorder(pg, c("pgID", "pgChr", "pgOrd", "og"))

    return(pg)
  }

  ##############################################################################
  # -- ad hoc function to pull non-syntenic orthologs
  parse_orthoList <- function(pangenomeDir, bed, repOfIDs){
    gids <- unique(bed$genome)
    fs <- file.path(pangenomeDir, sprintf("%s_compiledOrthologs.txt", gids))
    out <- rbindlist(lapply(fs, function(x){
      d <- fread(x, fill = T, sep = "\t", showProgress = FALSE)
      d <- subset(d, ofID %in% repOfIDs)
      d <- subset(d, !duplicated(ofID))
      d <- d[,list(orthos = strsplit(orthos, "|", fixed = T)[[1]]), by = "ofID"]

      ogv <- bed$og; names(ogv) <- bed$ofID
      d[,sameOg := ogv[ofID] == ogv[orthos]]
      d <- subset(d, !sameOg)[,c("ofID", "orthos")]
      d <- subset(d, !duplicated(d))
      setnames(d, "orthos", "nsOrthos")
      return(d)
    }))
    out <- subset(out, !duplicated(out))
    tmp <- bed[,c("genome", "ofID", "id")]
    setnames(tmp, "ofID", "nsOrthos")
    out <- merge(out, tmp, by = "nsOrthos")
    out[,`:=`(isNSOrtho = TRUE, isArrayRep = FALSE, isRep = FALSE)]
    return(out)
  }
  ##############################################################################
  ##############################################################################
  # 1. Read in the required datasets

  genome1 <- genome2 <- hasRef <- genome <- hasRefChr <- interpGenome <-
    interpChr <- isArrayRep <- noAnchor <- og <- pgID <- pgOrd <- ord <-
    pgChr <- chr <- diffChr <- notArrayRep <- ordDiff <- isRep <- repGenome <-
    nGenome <- isNSOrtho <- flag <- gen1 <- id1 <- pgID <- pgID2 <- id <- id2 <-
    gen2 <- ofID <- medord <- id <- repGene <- NULL

  ##############################################################################
  # -- 1.1 get parameters in order
  if(is.null(refGenome))
    refGenome <- gsParam$genomeIDs[1]
  if(is.null(genomeIDs))
    genomeIDs <- c(refGenome, gsParam$genomeIDs[gsParam$genomeIDs != refGenome])

  synBuff <- gsParam$params$synBuff * sqrt(2)
  nCores <- gsParam$params$nCores
  # find bed file here
  bedFile <- file.path(gsParam$paths$results, "combBed.txt")

  # write pangenome here
  pgFile <- file.path(
    gsParam$paths$pangenome,
    sprintf("%s_refPangenomeAnnot.txt", refGenome))

  # # -- interpolated position files
  # spFiles <- file.path(gsParam$paths$pangenome, sprintf(
  #   "%sintegratedSynPos.txt", genomeIDs))

  ##############################################################################
  # -- 1.2 read in the bed file
  if(verbose)
    cat(sprintf("Building pangenome against %s coordinates\n", refGenome))
  bed <- read_combBed(bedFile)

  # get list of genes with synteny information
  bed[,hasRef := refGenome %in% genome & uniqueN(genome) > 1, by = "og"]
  bed[,hasRefChr := any(hasRef), by = c("genome", "chr")]
  u <- with(subset(bed, hasRefChr & genome == refGenome), unique(chr))

  ##############################################################################
  # -- 1.3 read in the interpolated position data
  interp <- rbindlist(mclapply(genomeIDs, mc.cores = nCores, function(i)
    subset(read_intSynPos(file.path(gsParam$paths$pangenome, sprintf(
      "%s_vs_%s.integratedSynPos.txt", i, refGenome))), interpChr %in% u)))

  if(verbose)
    cat(sprintf("\t%s of %s genes have interpolated positions against %s\n",
                sum(bed$ofID %in% interp$ofID), nrow(bed), refGenome))

  ##############################################################################
  # -- 1.4 compile orthologs
  ocFiles <- file.path(gsParam$paths$pangenome,
                       sprintf("%s_compiledOrthologs.txt", genomeIDs))
  names(ocFiles) <- genomeIDs
  if(!all(file.exists(ocFiles))){
    if(verbose)
      cat("\tcompiling raw orthologs ... \n")
    oFileList <- mclapply(genomeIDs, mc.cores = nCores,  function(i)
      list.files(gsParam$paths$results, pattern = sprintf("%s__v", i), full.names = T))
    names(oFileList) <- genomeIDs

    idv <- bed$ofID; names(idv) <- with(bed, paste(genome, id))
    labs <- align_charLeft(genomeIDs)
    names(labs) <- genomeIDs
    for(i in genomeIDs){
      if(verbose)
        cat(sprintf("\t%s: ", labs[i]))

      x <- oFileList[[i]]
      ogs <- rbindlist(mclapply(x, mc.cores = nCores, parse_orthologues))
      ogs[,`:=`(ofID1 = idv[paste(gen1, id1)],
                ofID2 = idv[paste(gen2, id2)])]
      ogo <- ogs[,list(orthos = list(unique(ofID2))), by = "ofID1"]
      setnames(ogo, "ofID1", "ofID")
      fwrite(ogo, file = ocFiles[i], showProgress = F, sep = "\t")

      if(verbose)
        cat(sprintf("%s genes with %s orthologs\n", nrow(ogo), nrow(ogs)))
    }
  }
  ##############################################################################
  ##############################################################################
  # 2. build pangenome skeleton

  ##############################################################################
  # -- 2.1 cluster interpolated positions into anchored skeleton
  pgAnch <- clus_interp2pg(
    bed = subset(bed, isArrayRep & !noAnchor),
    interp = interp, minScaleProp = propAssignThresh,
    synBuff = synBuff)

  ##############################################################################
  # -- 2.2 merge with the bed file by orthogroup
  bd <- bed[,c("og", "ofID")]
  bd <- subset(bd, !duplicated(bd))
  pga <- subset(pgAnch, !duplicated(pgAnch))
  pgl <- merge(pga, bd, by = "og", allow.cartesian = TRUE)

  ##############################################################################
  # -- 2.3 cluster non-anchor genes and combine
  pgRem <- clus_interp2pg(
    bed = subset(bed, !og %in% unique(pgl$og)),
    interp = interp,  minScaleProp = propAssignThresh,
    synBuff = synBuff)

  bd <- bed[,c("og", "ofID")]
  bd <- subset(bd, !duplicated(bd))
  pga <- subset(pgRem, !duplicated(pgRem))
  pglr <- merge(pga, bed[,c("og", "ofID")], by = "og", allow.cartesian = TRUE)
  pglr[,pgID := sprintf("%s_remain", pgID)]
  pgl <- rbind(pgl, pglr)
  setkey(pgl, pgOrd)

  ##############################################################################
  # -- 2.4 Add remaining orthogroups without positions

  # pull out the missing data from the bed
  bedMiss <- subset(bed, !og %in% pgl$og)

  # pull refgenome missing data
  bedRef <- subset(bedMiss, genome == refGenome)
  bedRef <- subset(bed, og %in% bedRef$og)

  # make skeleton with missing, non-ref etc data
  pgr <- with(bedRef, data.table(
    og = og, pgChr = chr, pgOrd = ord, ofID = ofID))
  pgr[,pgID := sprintf("%s_missRef", og)]
  pgm <- with(subset(bedMiss, !og %in% bedRef$og), data.table(
    og = og, pgID = og, pgChr = NA, pgOrd = NA, ofID = ofID))
  pgm[,pgID := sprintf("%s_missNon", og)]

  # combine skeleton
  pgl <- rbind(pgl, pgr, pgm)
  setorder(pgl, pgOrd, na.last = T)
  pgl[,pgID := sprintf(
    "pg_%s", as.numeric(factor(pgID, levels = unique(pgID))))]

  ##############################################################################
  # -- 2.5 Add in genome and array rep info
  # combine skeleton and bed
  pgla <- subset(pgl, !duplicated(pgl))
  bd <-  bed[,c("ofID", "genome", "isArrayRep", "ord", "chr", "id")]
  bd <- subset(bd, !duplicated(bd))
  pglo <- merge(pgla, bd, all.x = TRUE, allow.cartesian = TRUE)
  genomeOrd <- genomeIDs

  # determine the representative gene by several nested factors
  pglo[,`:=`(genome = factor(genome, levels = genomeOrd),
             ordDiff = abs(ord - pgOrd),
             diffChr = pgChr == chr,
             notArrayRep = !isArrayRep,
             pgGenome = refGenome)]
  setorder(pglo, pgID, genome, diffChr, notArrayRep, ordDiff, na.last = T)
  pglo[,isRep := c(TRUE, rep(FALSE, .N-1)), by = "pgID"]

  # add in other necessary columns
  pglo[,repGenome := genome[isRep], by = "pgID"]
  pglo[,nGenome := uniqueN(genome), by = "pgID"]
  pglo[,isNSOrtho := FALSE]

  # reformat
  pglo <- pglo[,c("pgID", "pgGenome", "pgChr", "pgOrd", "genome", "og",
                  "isRep","ofID", "id", "isNSOrtho", "isArrayRep")]
  setorder(pglo, pgOrd, na.last = T)

  ##############################################################################
  ##############################################################################
  # 3 flag and reformat to return
  # non-syntenic orthologs are added if:
  #     1. they are orthologous to a isRep gene
  #     2. they are not in the same orthogroup as isRep gene
  #     3. they are not themselves a rep with the isRep gene in that orthogroup

  ##############################################################################
  # -- 3.1 read and parse ortholog tables to build non-syntenic orthologs
  nsOrthos <- parse_orthoList(
    pangenomeDir = gsParam$paths$pangenome,
    bed = bed,
    repOfIDs = pglo$ofID[pglo$isRep])

  ##############################################################################
  # -- 3.2 get the data in format to merge with pangenome
  pgtmp <- subset(pglo, isRep & ofID %in% nsOrthos$ofID)
  pgtmp <- pgtmp[,c("pgID", "pgGenome", "pgChr", "pgOrd", "og", "ofID")]
  pgtmp <- subset(pgtmp, !duplicated(pgtmp))
  nso <- merge(pgtmp, nsOrthos, by = "ofID")
  nso[,ofID := NULL]
  nso <- merge(bed[,c("genome", "id", "ofID")], nso, by = c("genome", "id"))
  nsOrthos <- nso[,colnames(pglo), with = F]

  ##############################################################################
  # -- 3.3 drop and NS orthos that hit too many places
  if(is.finite(maxNonSynPlaces)){
    npl <- gsParam$ploidy * maxNonSynPlaces
    nsOrthos[,nPlaces := uniqueN(pgID), by = "ofID"]
    nsOrthos[,dropThis := nPlaces > npl[genome]]
    if(verbose)
      cat(sprintf(
        "\tdropped %s non-syntenic orthologs that hit > `%s * ploidy` places\n",
        sum(nsOrthos$dropThis), maxNonSynPlaces))
    nsOrthos <- subset(nsOrthos, !dropThis)
    nsOrthos[,`:=`(dropThis = NULL, nPlaces = NULL)]
  }

  ##############################################################################
  # -- 3.4 combine
  pgout <- rbind(pglo, nsOrthos)
  setorder(pgout, pgOrd, na.last = T)

  ##############################################################################
  # -- 3.5 re-format and return
  pgw <- data.table(pgout)
  pgw[,flag := ifelse(isNSOrtho, "*", ifelse(!isArrayRep, "+", ""))]
  pgw[,id := sprintf("%s%s", id, flag)]
  pgw[,repGene := id[isRep][1], by = "pgID"]
  pgw <- dcast(pgw, pgID + pgGenome + pgChr + pgOrd + og + repGene ~ genome,
               value.var = "id", fun.aggregate = list)
  setorder(pgw, pgOrd, pgChr, pgGenome, pgID, na.last = TRUE)
  write_pangenome(pgout, filepath = pgFile)

  if(verbose)
    with(pgout, cat(sprintf(
      "Built pangenome annotation with ...\n\t%s positions\n\t%s syntenic orthogroup anchors\n\t%s tandem array members\n\t%s non-syntenic orthologs\n",
      uniqueN(pgID), sum(isArrayRep & !isNSOrtho),
      sum(!isArrayRep & !isNSOrtho), sum(isNSOrtho))))
  return(pgw)
}
