#' @title Build genespace pangenome
#'
#' @description
#' \code{pangenome} Convert orthogroup and synteny information into a
#' pangenome database. Predict locations of orthogroups that are missing a
#' node in the reference.
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param refGenome character string matching one of the genomeIDs in gsParam
#' @param genomeIDs character vector, specifying which genomes to use. Defaults
#' to all genomeIDs specification in gsParam.
#' @param propAssignThresh numeric of length 1, the minimum proportion of genes
#' in a pangenome entry needed to anchor the physical position to a new location
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
                      verbose = T){

  clus_interp2pg <- function(bed, interp, minScaleProp, synBuff){
    bd <- bed[,c("ofID", "og")]
    bd <- subset(bd, !duplicated(bd))
    intp <- interp[,c("ofID", "interpGenome", "interpChr", "interpOrd", "isAnchor")]
    intp <- subset(intp, !duplicated(intp))

    tmp <- merge(bd, intp, by = "ofID", allow.cartesian = TRUE)

    # -- cluster nearby hits
    setkey(tmp, og, interpChr, interpOrd, isAnchor)
    tmp[,jumpLeft := c(synBuff + 1, diff(interpOrd)), by = c("interpChr", "og")]
    tmp[,clus := as.integer(jumpLeft > synBuff), by = c("interpChr")]
    tmp[,clus := cumsum(clus), by = c("interpChr", "og")]
    tmp[,clus := paste(interpChr, clus)]
    tmp[,nClus := uniqueN(clus), by = "og"]

    # -- alculate the proportions in each placement, dropping those with < 25%
    tmp[,nOgMem := uniqueN(ofID), by = "og"]
    tmp[,nClusMem := uniqueN(ofID), by = c("og", "clus")]
    tmp[,propInClus := nClusMem/nOgMem]
    tmp[,scaleProp := scale_between(x = propInClus, min = 0, max = 1, scale1toMean = FALSE), by = "og"]
    mNopos <- subset(tmp, scaleProp < minScaleProp)
    tmp <- subset(tmp, scaleProp >= minScaleProp)

    # --- get median positions for each orthogroup
    pg <- tmp[,list(pgOrd = median(interpOrd)), by = c("og", "interpChr", "clus")]
    setkey(pg, pgOrd)
    pg[,pgID := sprintf("pg_%s", 1:.N)]
    setnames(pg, "interpChr", "pgChr")
    pg[,clus := NULL]
    setcolorder(pg, c("pgID", "pgChr", "pgOrd", "og"))

    return(pg)
  }

  ##############################################################################
  ##############################################################################
  # 1. Read in the required datasets

  ##############################################################################
  # -- 1.1 get parameters in order
  if(is.null(refGenome))
    refGenome <- gsParam$genomeIDs[1]
  if(is.null(genomeIDs))
    genomeIDs <- c(refGenome, gsParam$genomeIDs[gsParam$genomeIDs != refGenome])

  synBuff <- gsParam$params$synBuff * sqrt(2)

  # find bed file here
  bedFile <- file.path(gsParam$paths$results, "combBed.txt")

  # write pangenome here
  pgFile <- file.path(
    gsParam$paths$pangenome,
    sprintf("%s_refPangenomeAnnot.txt", refGenome))

  # ortholog files
  orthoFiles <- data.table(CJ(
    genome1 = genomeIDs,
    genome2 = genomeIDs))
  orthoFiles[,file := file.path(gsParam$paths$results,
                                sprintf("%s__v__%s.tsv", genome1, genome2))]
  orthoFiles <- subset(orthoFiles, file.exists(file))

  # -- interpolated position files
  spFiles <- file.path(gsParam$paths$pangenome, sprintf(
    "%sintegratedSynPos.txt", genomeIDs))

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
  interp <- rbindlist(lapply(spFiles, function(x)
    subset(read_intSynPos(x),
           interpGenome == refGenome & interpChr %in% u)))
  if(verbose)
    cat(sprintf("\t%s of %s genes have interpolated positions against %s\n",
                sum(bed$ofID %in% interp$ofID), nrow(bed), refGenome))

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
  pglo <- pglo[,c("pgID", "pgGenome", "pgChr", "pgOrd", "genome", "og", "isRep","ofID", "id", "isNSOrtho", "isArrayRep")]
  setorder(pglo, pgOrd, na.last = T)

  ##############################################################################
  ##############################################################################
  # 3 deal with non-syntenic orthologs

  ##############################################################################
  # -- 3.1 pull all pangenome entries that should be considered

  # drop entries with a single location
  pglo[,flag := uniqueN(genome) == 1 & genome[isRep] != refGenome, by = "pgID"]
  pg2x <- subset(pglo, !flag)
  pg1x <- subset(pglo, flag)
  pg2x[,flag := NULL]
  pg1x[,flag := NULL]
  pglo[,flag := NULL]

  # pull all genes that should be considered as potential anchors
  anchGenes <- with(subset(pglo, isRep), unique(paste(genome, id)))

  ##############################################################################
  # -- 3.2 make data to complete the merge and add data to the orthologs
  pgi1 <- subset(pglo, isRep)
  pgi2 <- with(pglo, data.table(
    id2 = id, gen2 = genome, pgID2 = pgID))
  pgi2 <- subset(pgi2, !duplicated(pgi2))
  ov <- bed$ofID; names(ov) <- with(bed, paste(genome, id))

  ##############################################################################
  # -- 3.3 query each ortholog file
  nsOrthos <- rbindlist(lapply(orthoFiles$file, function(i){
    tmp <- parse_orthologues(i)
    tmp <-  subset(tmp, paste(gen1, id1) %in% anchGenes)
    setnames(tmp, c("gen1", "id1"), c("genome", "id"))

    pgi1 <- subset(pgi1, !duplicated(pgi1))
    pgi2 <- subset(pgi2, !duplicated(pgi2))
    tmp <- subset(tmp, !duplicated(tmp))

    tmp <- merge(
      pgi1,
      merge(pgi2, tmp, by = c("gen2", "id2"), allow.cartesian = TRUE),
      by = c("genome", "id"), allow.cartesian = TRUE)
    tmp <- subset(tmp, pgID != pgID2)
    tmp <- subset(tmp, !duplicated(paste(id, id2)))
    tmp[,`:=`(ofID = ov[paste(gen2, id2)], id = id2, genome = gen2,
              isNSOrtho = TRUE)]
    tmp <- tmp[,colnames(pglo), with = F]
    return(tmp)
  }))
  nsOrthos[,isRep := FALSE]

  ##############################################################################
  # -- 3.4 put real interpolated positions back in
  nsi <- interp[,c("ofID", "interpChr", "interpOrd")]
  setnames(nsi, c("ofID", "pgChr", "pgOrd"))
  nsOrthos[,`:=`(pgChr = NULL, pgOrd = NULL)]
  nsOrthos <- subset(nsOrthos, !duplicated(nsOrthos))
  nsi <- subset(nsi, !duplicated(nsi))
  nsOrthos <- merge(nsOrthos, nsi, by = "ofID",
                    all.x = T, allow.cartesian = TRUE)

  ##############################################################################
  # -- 3.5 drop 1x position orthogroups that have non-syntenic orthos to 2x+ ogs
  ns2x <- subset(nsOrthos, pgID %in% pg2x$pgID)
  pg2x <- rbind(pg2x, ns2x)
  ns1x <- subset(nsOrthos, pgID %in% pg1x$pgID[!pg1x$ofID %in% pg2x$ofID])
  pg1x[,flag := ofID %in% pg2x$ofID]
  pg1x[,flag := all(flag), by = "pgID"]
  pg1x <- subset(pg1x, !flag)
  pg1x[,flag := NULL]
  pg1x <- rbind(pg1x, ns1x)

  # combine
  pgout <- rbind(pg1x, pg2x)
  setorder(pgout, pgOrd, na.last = T)

  ##############################################################################
  # -- 3.6 check theres not a problem and some NSOrthos are actually syntenic
  cv <- bed$chr; names(cv) <- bed$ofID
  pgout[,chr := cv[ofID]]
  pgout[, medord := pgOrd[isRep], by = "pgID"]
  pgout[,`:=`(ordDiff = abs(medord - pgOrd),
             diffChr = pgChr == chr,
             notArrayRep = !isArrayRep)]
  setorder(pgout, pgID, ofID, diffChr, notArrayRep, ordDiff, isNSOrtho, na.last = T)
  pgout <- subset(pgout, !duplicated(paste(pgID, ofID)))

  # set those genes with interpolated positions close enough to not NS
  pgout$isNSOrtho[with(pgout, isNSOrtho & ordDiff < synBuff & !diffChr)] <- FALSE
  pgout <- pgout[,colnames(pg1x), with = F]
  pgout[,`:=`(pgChr = pgChr[isRep], pgOrd = pgOrd[isRep]), by = "pgID"]

  # rename
  setorder(pgout, pgOrd, pgID, -isRep, isNSOrtho, -isArrayRep, na.last = T)
  pgout[,pgID := as.numeric(factor(pgID, levels = unique(pgID)))]
  pgout[,pgID := sprintf("pg_%s", gsub(" ", "0", align_charRight(pgID)))]
  pgout$isArrayRep[pgout$isRep] <- TRUE

  ##############################################################################
  ##############################################################################
  # 4 flag and reformat to return
  pgw <- data.table(pgout)
  pgw[,flag := ifelse(isNSOrtho, "*", ifelse(!isArrayRep, "+", ""))]
  pgw[,id := sprintf("%s%s", id, flag)]
  pgw[,repGene := id[isRep][1], by = "pgID"]
  pgw <- dcast(pgw, pgID + pgGenome + pgChr + pgOrd + og + repGene ~ genome,
               value.var = "id", fun.aggregate = list)
  fwrite(pgout, file = pgFile, sep = "\t")

  if(verbose)
    with(pgout, cat(sprintf(
      "Built pangenome annotation with ...\n\t%s positions\n\t%s syntenic orthogroup anchors\n\t%s tandem array members\n\t%s non-syntenic orthologs\n",
      uniqueN(pgID), sum(isArrayRep & !isNSOrtho), sum(!isArrayRep & !isNSOrtho), sum(isNSOrtho))))
  return(pgw)
}
