#' @title Combine and annotate bed files
#' @description
#' \code{annotate_bed} Combines all orthogroup and exclustion data for each
#' gene
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @export
annotate_bed <- function(gsParam){

  ##############################################################################
  # 0. Add hoc functions
  # -- 0.1 add orthofinder ID to the bed file
  add_ofID2bed <- function(bed, gsParam){
    genomeNum <- genome <- id <- ofID <- NULL
    seqid <- read_orthofinderSequenceIDs(
      file.path(gsParam$paths$results, "SequenceIDs.txt"))
    speid <- read_orthofinderSpeciesIDs(
      file.path(gsParam$paths$results, "SpeciesIDs.txt"))
    seqid[,genome := names(speid)[match(genomeNum, speid)]]
    di <- seqid$ofID; names(di) <- with(seqid, paste(genome, id))
    bed[,ofID := di[paste(genome, id)]]
    return(bed)
  }

  # -- 0.2 add peptide lengths to the bed file
  add_pepLen2bed <- function(bed, gsParam){
    id <- pepLen <- NULL
    genomeIDs <- unique(bed$genome)
    pepFiles <- file.path(gsParam$paths$peptide, sprintf("%s.fa", genomeIDs))
    pepFiles <- pepFiles[file.exists(pepFiles)]
    if(length(pepFiles) > 1){
      naa <- character()
      for(i in pepFiles){
        gid <- gsub(".fa$", "", basename(i))
        tmp <- get_nAA(i)
        names(tmp) <- paste(gid, names(tmp))
        naa <- c(naa, tmp)
      }
      bed[,pepLen := naa[paste(genome, id)]]
    }else{
      bed[,pepLen := NA]
    }
    return(bed)
  }

  # -- 0.3 add orthogroup IDs to the bed file
  add_og2bed <- function(bed, gsParam){
    genome <- id <- NULL
    ogPath <- file.path(gsParam$paths$results, "Orthogroups.tsv")
    ogs <- parse_ogs(filepath = ogPath, genomeIDs = genomeIDs)
    di <- ogs$ogID; names(di) <- with(ogs, paste(genome, id))
    bed[,globOG := di[paste(genome, id)]]
    wh <- which(is.na(bed$globOG))
    lab <- gsub(" ", "0", align_charRight(c("00000", 1:length(wh))))[-1]
    bed$globOG[wh] <- sprintf("NoOG%s", lab)
    return(bed)
  }

  # -- 0.4 add HOG IDs to the bed file
  add_hog2bed <- function(bed, gsParam){
    genome <- id <- NULL
    ogPath <- file.path(gsParam$paths$results, "N0.tsv")
    ogs <- parse_hogs(filepath = ogPath, genomeIDs = genomeIDs)
    di <- ogs$hogID; names(di) <- with(ogs, paste(genome, id))
    bed[,globHOG := di[paste(genome, id)]]
    wh <- which(is.na(bed$globHOG))
    lab <- gsub(" ", "0", align_charRight(c("00000", 1:length(wh))))[-1]
    bed$globHOG[wh] <- sprintf("N0.NoHOG%s", lab)
    return(bed)
  }

  # -- 0.5 add the number of places the orthogroup hits to the bed file
  add_nOGplaces2bed <- function(bed, arrayJump){
    genome <- chr <- og <- ord <- jumpLeft <- NULL
    setkey(bed, genome, chr, og, ord)
    bed[,jumpLeft := c(arrayJump + 1, diff(ord)), by = c("genome", "chr", "og")]
    bed[,nOGPlaces := sum(jumpLeft > arrayJump), by = c("genome","og")]
    bed[,`:=`(jumpLeft = NULL)]
    setkey(bed, genome, ord)
    return(bed)
  }

  # -- 0.6 add array ID to the bed file
  add_array2bed <- function(bed, synBuff, maxIter = 10, reorder = TRUE){

    # -- we can make arrays for everything except that we want to exclude the
    # arrays that are huge and problematic
    genome <- chr <- id <- arrayID <- nOGPlaces <- tord <- NULL
    tmp <- data.table(bed)
    tmp[,tord := as.numeric(ord)]

    # -- set up the iteration
    cnt <- 1
    diffn <- 1
    tmp[,arrayID := sprintf(
      "tmp%s", as.integer(as.factor(paste(genome, chr, id))))]
    while(cnt <= maxIter && diffn > 0){

      # -- for each iteration, calculate clusters by the size of jump between
      # genes larger than the synBuffer
      cnt <- cnt + 1
      initn <- uniqueN(tmp$arrayID)
      genome <- chr <- og <- ord <- jumpLeft <- clus <- n <- NULL
      setkey(tmp, genome, chr, og, tord)
      if(reorder)
        tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
      tmp[,jumpLeft := c(synBuff + 1, diff(tord)),
          by = c("genome", "chr", "og")]
      tmp[,clus := as.integer(jumpLeft > synBuff), by = c("genome", "chr")]
      tmp[,clus := cumsum(clus), by = c("genome", "chr", "og")]
      tmp[,`:=`(
        arrayID = sprintf("tmp%s",
                          as.integer(as.factor(paste(genome, chr, og, clus)))),
        jumpLeft = NULL, clus = NULL)]

      # -- get the new order of genes based on array ID
      tmp[,n := .N, by = "arrayID"]
      tmp1 <- subset(tmp, n == 1)
      tmp2 <- subset(tmp, n > 1)
      tmp2[,tord := as.numeric(tord)]
      tmp2[,tord := mean(as.numeric(tord)), by = "arrayID"]
      tmp <- rbind(tmp1, tmp2)
      tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
      newn <- uniqueN(tmp$arrayID)
      diffn <- initn - newn
    }

    # -- relabel
    ogID <- NULL
    lab <- gsub(" ", "0",
                align_charRight(
                  as.numeric(factor(tmp$arrayID,
                                    levels = unique(tmp$arrayID)))))
    arrv <- sprintf("Arr%s", lab); names(arrv) <- tmp$ofID
    bed[,arrayID := arrv[ofID]]

    tmp <- subset(bed, is.na(arrayID))
    wh <- with(tmp,
               as.numeric(as.factor(paste(genome, chr, og))))
    lab <- gsub(" ", "0", align_charRight(wh))
    wh <- which(is.na(bed$arrayID))
    bed$arrayID[wh] <- sprintf("NoArr%s", lab)
    return(bed)
  }

  # -- 0.7 calculate the array representatives
  add_arrayReps2bed <- function(bed){
    n <- ord <- medOrd <- medDiff <- pepLen <- arrayID <- isRep <-
      isArrayRep <- ofID <- NULL
    bed[,n := .N, by = "arrayID"]
    tmp <- subset(bed, n > 1)
    tmp[,medOrd := as.numeric(median(ord)), by = "arrayID"]
    tmp[,medDiff := as.numeric(abs(medOrd - ord))]
    setorder(tmp, arrayID, medDiff, -pepLen)
    tmp$noAnchor[duplicated(tmp$arrayID)] <- TRUE
    tmp[,isRep := !duplicated(arrayID)]
    bed[,isArrayRep := ofID %in% tmp$ofID[tmp$isRep] | n == 1]
    tmp <- subset(bed, isArrayRep)
    tmp[,ord := frank(ord, ties.method = "dense"), by = "genome"]
    di <- tmp$ord; names(di) <- tmp$ofID
    return(bed)
  }

  # -- 0.8 get summary stats for og/hog matching with ploidy
  check_ploidyOgMatch <- function(ploidy, bed){
    genome <- nGenome <- nGenes <- ofID <- n <- nOGs <- smPloid <- NULL
    cat("\t##############\n\tChecking which type of orthogroup best matches ploidy\n")
    tmp <- data.table(genome = names(ploidy), ploidy = ploidy)
    bed[,nGenome := uniqueN(genome), by = "globOG"]
    cat(sprintf("\t%s OGs (of %s) with %s genes (of %s) hit all genomes\n",
                uniqueN(bed$globOG[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$globOG),
                uniqueN(bed$ofID[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$ofID)))
    cnts <- subset(bed, nGenome == nrow(tmp))
    cnts <- cnts[,list(nGenes = uniqueN(ofID)), by = c("globOG", "genome")]
    cnts <- cnts[,list(n = .N), by = c("genome", "nGenes")]
    ogCnt <- merge(tmp, cnts, by = "genome", all = T)
    ogCnt[,smPloid := ploidy == nGenes]
    ogCnt <- ogCnt[,list(nOGs = sum(n), type = "OG"), by = c("genome", "smPloid")]

    bed[,nGenome := uniqueN(genome), by = "globHOG"]
    cat(sprintf("\t%s HOGs (of %s) with %s genes (of %s) hit all genomes\n",
                uniqueN(bed$globHOG[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$globHOG),
                uniqueN(bed$ofID[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$ofID)))
    cnts <- subset(bed, nGenome == nrow(tmp))
    cnts <- cnts[,list(nGenes = uniqueN(ofID)), by = c("globHOG", "genome")]
    cnts <- cnts[,list(n = .N), by = c("genome", "nGenes")]
    hogCnt <- merge(tmp, cnts, by = "genome", all = T)
    hogCnt[,smPloid := ploidy == nGenes]
    hogCnt <- hogCnt[,list(nOGs = sum(n), type = "HOG"), by = c("genome", "smPloid")]
    cnts <- rbind(ogCnt, hogCnt)
    tp <- dcast(cnts, genome ~ type + smPloid, value.var = "nOGs")
    tp <- rbind(tp,
                matrix(c(genome = "total",
                         colSums(tp[,-1,with = F])),nrow =  1), use.names = F)

    pOG <- pHOG <- genome <- HOG <- OG <- HOG_TRUE <- HOG_FALSE <- OG_TRUE <-
      OG_FALSE <- NULL
    tp[,`:=`(pHOG = as.numeric(HOG_TRUE)/(as.numeric(HOG_FALSE) + as.numeric(HOG_TRUE)),
             pOG = as.numeric(OG_TRUE)/(as.numeric(OG_FALSE) + as.numeric(OG_TRUE)))]
    tp[,`:=`(pHOG = round(pHOG*100, 2), pOG = round(pOG*100,2))]

    md <- tp[,list(genome = align_charLeft(c("#", sprintf("...%s",genome))),
                   HOG = align_charLeft(c("match, not, % (HOGs)",
                                          sprintf("%s, %s, %s", HOG_TRUE, HOG_FALSE, pHOG))),
                   OG = align_charLeft(c("match, not, % (OGs)",
                                         sprintf("%s, %s, %s", OG_TRUE, OG_FALSE, pOG))))]
    nu <- apply(md, 1, function(x){
      x <- unlist(x)
      cat(sprintf("\t%s: %s || %s\n", x[1], x[2], x[3]))
    })
    tmp <- subset(tp, genome == "total")
    if(tmp$pHOG > tmp$pOG){
      useHOGs <- TRUE
      cat("\t**NOTE** HOGs better match ploidy than OGs. Setting useHOGs = TRUE")
    }else{
      useHOGs <- FALSE
      cat("\t**NOTE** OGs better match ploidy than HOGs. Setting useHOGs = FALSE")
    }
    return(useHOGs)
  }

  # -- 0.9 add small chromosome flag to bed file
  add_smallChr2bed <- function(bed, blkSize){
    smallChr <- nGenes <- nChrs <- NULL
    cat(sprintf(
      "\n\t##############\n\tFlagging chrs. w/ < %s unique orthogroups\n",
      blkSize * 2))
    bed[,smallChr := uniqueN(og) < (blkSize * 2), by = c("genome", "chr")]
    tab <- bed[,list(nGenes = sum(smallChr),
                     nChrs = uniqueN(chr[smallChr]),
                     flg = ifelse(sum(smallChr)/.N > 0.05, "***", "")),
               by = "genome"]
    tab[,`:=`(lab = align_charLeft(genome),
              nGenes = align_charRight(nGenes),
              nChrs = align_charRight(nChrs))]
    with(tab, cat(sprintf(
      "\t...%s: %s genes on %s small chrs. %s\n", lab, nGenes, nChrs, flg)))
    anyBad <- any(tab$flg == "***")
    if(anyBad)
      cat(strwrap(
        "NOTE! Genomes flagged *** have > 5% of genes on small chrs. These are
      likely not great assemblies and should be examined carefully",
        indent = 8, exdent = 16), sep = "\n")
    return(bed)
  }

  # -- 0.10 add overdispersed og flag to bedfile
  add_dispOG2bed <- function(bed, arrayJump){
    dispOG <- nGenes <- nOgs <- NULL
    cat("\t##############\n\tFlagging over-dispered OGs\n")
    bed <- add_nOGplaces2bed(bed = bed, arrayJump = arrayJump)
    bed[,dispOG := nOGPlaces > maxOgPlaces[genome]]

    tab <- bed[,list(nGenes = sum(dispOG),
                     nOgs = uniqueN(og[dispOG]),
                     flg = ifelse(sum(dispOG)/.N > 0.05, "***", "")),
               by = "genome"]
    tab[,`:=`(lab = align_charLeft(genome),
              nGenes = align_charRight(nGenes),
              nOgs = align_charRight(nOgs))]
    with(tab, cat(sprintf(
      "\t...%s: %s genes in %s OGs hit > %s unique places %s\n",
      lab, nGenes, nOgs, maxOgPlaces[genome], flg)))
    anyBad <- any(tab$flg == "***")
    if(anyBad)
      cat(strwrap(
        "NOTE! Genomes flagged *** have > 5% of genes in over-dispersed
      orthogroups. These are likely not great annotations and should be examined
      carefully",
        indent = 8, exdent = 16), sep = "\n")
    bed[,nOGPlaces := NULL]
    return(bed)
  }

  ##############################################################################
  # 1. get the parameters in order

  # -- shorten the names of the parameters
  resDir      <- gsParam$paths$results
  genomeIDs   <- gsParam$genomeIDs
  useHOGs     <- gsParam$params$useHOGs
  nGaps       <- gsParam$params$nGaps
  synBuff     <- gsParam$params$synBuff
  blkSize     <- gsParam$params$blkSize
  maxOgPlaces <- gsParam$ploidy * 8
  arrayJump   <- gsParam$params$synBuff / 2
  nCores      <- gsParam$params$nCores
  ploidy      <- gsParam$ploidy

  # -- bed files
  bedFiles <- file.path(gsParam$paths$bed, sprintf("%s.bed", genomeIDs))
  if(!all(file.exists(bedFiles)))
    stop(sprintf("could not find all bed files in %s\n", bedDir))

  # -- output file
  combBedFile <- file.path(gsParam$paths$results, "combBed.txt")

  # -- null out env variables
  chrn <- chr <- ord <- start <- bed <- noAnchor <- arrayID <- ofID <- pepLen <-
    isArrayRep <- globOG <- globHOG <- synOG <- inblkOG <- og <- globHOG <-
    genome <- nOGPlaces <- tmp <- n <- maxPlaces <- bedDir <- smallChr <-
    dispOG <- nGenesInArray <- nArr <- nGenes <- nOGs <- NULL

  ##############################################################################
  # 2. combine bed files
  # -- 2.1 read in the bed files
  bed <- rbindlist(mclapply(bedFiles, mc.cores = nCores, function(i){
    tmp <- fread(
      i, col.names = c("chr", "start", "end", "id"),
      colClasses = c("character", "integer", "integer", "character"),
      showProgress = FALSE, verbose = FALSE)
    tmp[,genome := gsub(".bed$", "", basename(i))]
    return(tmp)
  }))

  # -- 2.3 add in gene orders
  bed[,chrn := as.numeric(gsub('\\D+','', chr))]
  bed[,n := .N, by = c("genome", "chr")]
  setorder(bed, genome, chrn, chr, -n, start)
  bed[,ord := 1:.N, by = "genome"]

  # -- 2.4 add in all other columns
  bed[,`:=`(
    chrn = NULL, n = NULL, ofID = NA, pepLen = NA, arrayID = NA,
    isArrayRep = FALSE,  globOG = NA, globHOG = NA, synOG = NA,
    inblkOG = NA, noAnchor = FALSE, smallChr = FALSE, dispOG = FALSE)]

  ##############################################################################
  # 3. Add in necesary data columns if they do not exist
  # -- 3.1 orthofinder IDs
  bed <- add_ofID2bed(bed = bed, gsParam = gsParam)

  # -- 3.2 peptide lengths
  bed <- add_pepLen2bed(bed = bed, gsParam = gsParam)

  # -- 3.3 orthogroup IDs
  bed <- add_og2bed(bed = bed, gsParam = gsParam)

  # -- 3.4 hier orthogroup IDs
  bed <- add_hog2bed(bed = bed, gsParam = gsParam)

  ##############################################################################
  # 4. Arrays and multi-placed orthogroups
  if(is.null(useHOGs) || is.na(useHOGs))
    useHOGs <- check_ploidyOgMatch(bed = bed, ploidy = ploidy)
  if(useHOGs){
    bed[,og := globHOG]
  }else{
    bed[,og := globOG]
  }

  # -- 4.1 Flag problematic chromosomes with < blkSize * 2 Ogs
  bed <- add_smallChr2bed(bed = bed, blkSize = blkSize)

  # -- 4.2 n. orthogroup placements
  bed <- add_dispOG2bed(bed = bed, arrayJump = arrayJump)
  bed[,noAnchor := smallChr | dispOG]

  # -- 4.3 find the arrays
  beda <- add_array2bed(
    bed = subset(bed, !noAnchor),
    synBuff = synBuff,
    maxIter = 10,
    reorder = T)

  bedi <- add_array2bed(
    bed = subset(bed, noAnchor),
    synBuff = synBuff,
    maxIter = 10,
    reorder = F)

  bedi[,arrayID := sprintf("%s_noAnch", arrayID)]
  bed <- rbind(beda, bedi)
  setkey(bed, genome, ord)

  # -- 4.4 choose the array representative genes
  bed <- add_arrayReps2bed(bed)

  cat("\t##############\n\tAnnotation summaries (after exclusions):\n")
  tab <- subset(bed, !noAnchor)[,list(
    nGenes = .N,
    nOGs = uniqueN(og),
    nArr = uniqueN(arrayID[!isArrayRep]),
    nGenesInArray = sum(!isArrayRep) + uniqueN(arrayID[!isArrayRep])),
    by = "genome"]
  tab[,`:=`(lab = align_charLeft(genome),
            nGenesInArray = align_charRight(nGenesInArray),
            nArr = align_charRight(nArr),
            nGenes = align_charRight(nGenes),
            nOGs = align_charRight(nOGs))]
  with(tab, cat(sprintf(
    "\t...%s: %s genes in %s OGs || %s genes in %s arrays\n",
    lab, nGenes, nOGs, nGenesInArray, nArr)))

  bed[,noAnchor := noAnchor | !isArrayRep]

  bed[,`:=`(
    nGenome = NULL, n = NULL)]

  write_combBed(
    x = bed, filepath = combBedFile)
  return(bed)
}
