#' @title generate syntenic pan-gene sets
#' @description
#' \code{syntenic_pangenes} Integrates interpolated gene order position with
#' orthogroups to generate pan-gene sets across multiple genomes.
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param refGenome character specifying the reference genome
#' @param maxPlacePerChr integer specifying the maximum number of placements
#' allowed on a single reference chromosome. This keeps things simpler in very
#' messy synteny maps.
#' @param minPropInterp2keep numeric specifying the minimum proportion of genes
#' in an orthogroup with positions clustered into a specific location.
#'
#' @details formerly called 'pangenome' this function decomposes orthogroups
#' and syntenic interpolated positions into a long-formatted pan-gene set. See
#' query_pangenes to access this output.
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom dbscan dbscan frNN
#' @importFrom stats median
#' @export
syntenic_pangenes <- function(gsParam,
                              refGenome,
                              maxPlacePerChr = 2,
                              minPropInterp2keep = .75){

  cluster_anchors <- function(gsParam,
                              bed,
                              refGenome,
                              maxPlacePerChr,
                              minPropInterp2keep){

    noAnchor <- isArrayRep <- ofID2 <- genome <- ofID <- ord <- nInChr <-
      clus <- ntot <- ninclus <- prop <- n <- rnk <- bad <- propOfMax <-
      og <- chr <- interpOrd <- pgID <- flag <- globHOG <- sameChr <-
      interpChr <- wt <- pgRepID <- nInArray <- glevs <- NULL

    bedAll <- data.table(bed)
    bed <- subset(bed, !noAnchor & isArrayRep)
    ogv <- bed$og; names(ogv) <- bed$ofID
    narrv <- bed$nInArray; names(narrv) <- bed$ofID

    # -- 1 convert to a scaffold
    scaff <- data.table(bed)
    scaff <- scaff[,c("genome","ofID", "chr", "ord", "og")]

    # -- 2 read interpolated positions
    ipFile <- file.path(gsParam$paths$pangenes,
                        sprintf("%s_integratedSynPos.txt", refGenome))
    synPos <- subset(read_intSynPos(ipFile), ofID2 %in% scaff$ofID)

    # -- 3 pull synPos for anchors
    synPos <- subset(synPos, ofID2 %in% scaff$ofID)
    synPos <- with(subset(synPos, !is.na(ofID2)), data.table(
      genome = genome2, ofID = ofID2, chr = chr1, ord = ord1, og = ogv[ofID2]))
    refa <- rbind(subset(scaff, genome == refGenome), synPos)
    refa <- subset(refa, !duplicated(paste(ofID, ord)))
    refa[,nInArray := narrv[ofID]]
    # -- 4 add ogs to interpolated positions
    refa[,nInChr := sum(!is.na(ord)), by =  c("chr", "og")]
    n1 <- subset(refa, nInChr <= 1)
    n2 <- subset(refa, nInChr > 1)
    n2[,diff := max(abs(diff(ord[!is.na(ord)]))), by = c("chr", "og")]

    # -- 5 drop syntenic placements > maxPlacePerChr/chr
    n3 <- subset(n2, diff > gsParam$params$synBuff * 2)
    n2[,`:=`(clus = 1, nInChr = NULL, diff = NULL)]
    n1[,`:=`(clus = 1, nInChr = NULL)]

    if(nrow(n3) > 0){
      n3[,clus := dbscan(frNN(
        cbind(ord, ord),
        eps = gsParam$params$synBuff*2), minPts = 0)$cluster,
        by = c("chr", "og")]

      n3[,ntot := uniqueN(ofID), by = c("og")]
      n3[,ninclus := uniqueN(ofID), by = c("chr", "og", "clus")]
      n3[,prop := ninclus / ntot]
      n3 <- subset(n3, prop >= minPropInterp2keep)
      n3[,n := uniqueN(clus), by = c("chr", "og")]
      n3[,rnk := frank(prop), by = c("chr", "og")]
      n3[,bad := n > maxPlacePerChr & rnk == 1]
      n3 <- subset(n3, !bad)

      n3[,`:=`(ntot = NULL, ninclus = NULL, prop = NULL,
               bad = NULL, rnk = NULL, n = NULL, nInChr = NULL, diff = NULL)]
      n3 <- rbind(n2, n3, n1)
    }else{
      n3 <- rbind(n2, n1)
    }
    n3 <- subset(n3, !duplicated(n3))

    # -- 6 drop placements < minPropInterp2keep
    n3[,ntot := uniqueN(ofID), by = c("og")]
    n3[,ninclus := uniqueN(ofID), by = c("chr", "og", "clus")]
    n3[,prop := ninclus / ntot]
    n3[,propOfMax := prop/max(prop), by = "og"]
    n3 <- subset(n3, prop >= minPropInterp2keep | propOfMax == 1)
    n3[,`:=`(ntot = NULL, ninclus = NULL, prop = NULL, propOfMax = NULL,
             genome = factor(genome, levels = glevs))]
    setorder(n3, genome, -nInArray)
    n3[,interpOrd := ord[1], by = c("og", "clus", "chr")]

    # -- 7 reformat

    n4 <- subset(n3, !duplicated(paste(og, clus, chr, interpOrd)))

    bedo <- subset(bed, og %in% unique(n3$og))
    bedo <- merge(
      bedo[,c("genome", "ofID", "og")],
      n4[,c("og", "clus", "interpOrd", "chr")],
      by = "og", allow.cartesian = T, all.x = T)

    bedo[,pgID := paste(chr, og, clus)]
    bedo <- bedo[,c("pgID", "chr", "interpOrd", "og", "genome", "ofID")]
    setkey(bedo, interpOrd, og, genome, pgID, ofID)
    bedo[,pgID := as.integer(factor(pgID, levels = unique(pgID)))]

    # -- 8 add in orthogroups that are missing syntenic positions
    u <- unique(bedAll$og[bedAll$isArrayRep])
    noPosOg <- subset(bedAll, og %in% u[!u %in% bedo$og])
    noPosOg[,`:=`(pgID = max(bedo$pgID) + as.numeric(as.factor(og)),
                  chr = NA, interpOrd = NA)]
    noPosOg <- noPosOg[,colnames(bedo), with = F]

    bedo <- rbind(bedo, noPosOg)
    setkey(bedo, pgID, genome)
    bedo[,flag := "PASS"]

    # -- 9 add bottom drawer
    bd <- subset(bedAll, !og %in% bedo$og)
    bd[,`:=`(pgID = max(bedo$pgID) + as.numeric(as.factor(globHOG)),
                  chr = NA, interpOrd = NA)]
    bd[,flag := "badOG"]
    bd <- bd[,colnames(bedo), with = F]

    bedo <- rbind(bedo, bd)
    setkey(bedo, pgID, genome)

    # -- 10 add pangenes rep genes
    setnames(bedo, "chr", "interpChr")
    bedo <- merge(bedo, bedAll[,c("ofID", "chr")],
                  by = "ofID", all.x = T, allow.cartesian = T)
    gids <- c(refGenome, gsParam$genomeIDs[gsParam$genomeIDs != refGenome])
    bedo[,genome := factor(as.character(genome), levels = gids)]
    bedo[,sameChr := chr == interpChr]
    setorder(bedo, pgID, genome, -sameChr, -flag, ofID, na.last = T)
    bedo[,wt := 1:.N, by = "pgID"]
    bedo[,pgRepID := ofID[wt == 1], by = "pgID"]
    bedo[,`:=`(chr = NULL, sameChr = NULL,  wt = NULL)]
    return(bedo)
  }

  isArrayRep <- flag <- gen1 <- id1 <- gen2 <- id2 <- ofID1 <- ofID2 <-
    pgID <- genome <- ofID <- nInArray <- NULL

  # 1. read in the combined bed file
  bed <- read_combBed(
    filepath = file.path(gsParam$paths$results, "combBed.txt"))
  bed[,nInArray := .N, by = "arrayID"]

  # 2. Cluster the anchors into pangenes entries
  pgScaff <- cluster_anchors(
    bed = bed,
    gsParam = gsParam,
    refGenome = refGenome,
    maxPlacePerChr = maxPlacePerChr,
    minPropInterp2keep = minPropInterp2keep)


  # 3. Add in the bottom drawer array members
  arrmem <- subset(bed, !isArrayRep)[,c("ofID", "og", "genome")]
  arrmem[,flag := "array"]
  tmp <- pgScaff[,c("pgID", "interpChr", "interpOrd", "og", "pgRepID")]
  tmp <- subset(tmp, !duplicated(tmp))
  pgArr <- merge(tmp, arrmem, by = "og", allow.cartesian = T)
  pgArr <- pgArr[,colnames(pgScaff), with = F]
  pgOut <- rbind(pgScaff, pgArr)

  # 4. Add non-syntenic orthologs
  dict <- bed$ofID; names(dict) <- with(bed, paste(genome, id))
  u <- unique(pgOut$pgRepID)
  hasu <- with(pgOut, paste(pgRepID, ofID))
  orthof <- unlist(gsParam$synteny$blast[,c("queryOrthologs", "targetOrthologs")])
  orthof <- orthof[!is.na(orthof)]
  ogs <- rbindlist(mclapply(orthof, mc.cores = gsParam$params$nCores, function(i){
    x <- parse_orthologues(i)
    x[,`:=`(ofID1 = dict[paste(gen1, id1)],
            ofID2 = dict[paste(gen2, id2)],
            gen1 = NULL, gen2 = NULL, id1 = NULL, id2 = NULL, orthID = NULL)]
    x <- subset(x, ofID1 %in% u)
    x <- subset(x, !paste(ofID1, ofID2) %in% hasu)
    return(x)
  }))
  setnames(ogs, c("pgRepID", "ofID"))
  ogs <- merge(bed[,c("genome", "og", "ofID")], ogs,
               by = "ofID", allow.cartesian = T)
  nsogs <- merge(
    pgOut[,c("pgID", "interpChr", "interpOrd", "pgRepID")],
    ogs, by = "pgRepID", allow.cartesian = T)
  nsogs[,flag := "NSOrtho"]

  pgOut <- rbind(pgOut, nsogs)
  setcolorder(pgOut, c(
    "pgID", "interpChr", "interpOrd", "pgRepID", "genome", "og", "ofID", "flag"))

  # 5. add in positions of genes and drop duplicates
  pgOut <- merge(pgOut, bed[,c("ofID", "id", "chr", "start", "end", "ord")],
                 by = "ofID", all.x = T, allow.cartesian = T)

  pgOut[,flag := factor(flag, levels = c("PASS", "array", "badOG", "NSOrtho"))]
  setkey(pgOut, pgID, flag, genome, ofID)
  pgOut <- subset(pgOut, !duplicated(paste(ofID, pgID)))

  outf <- file.path(gsParam$paths$pangenes,
                    sprintf("%s_pangenes.txt.gz", refGenome))
  write_pangenes(x = pgOut, filepath = outf)

  return(pgOut)
}
