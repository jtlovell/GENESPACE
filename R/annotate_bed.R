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
    ogs <- parse_hogs(filepath = ogPath)
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

  # -- 0.9 add small chromosome flag to bed file
  add_smallChr2bed <- function(bed, blkSize){
    smallChr <- nGenes <- nChrs <- NULL
    cat(sprintf(
      "\t##############\n\tFlagging chrs. w/ < %s unique orthogroups\n",
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
    bed[,dispOG := nOGPlaces > max(maxOgPlaces)]

    tab <- bed[,list(nGenes = sum(dispOG),
                     nOgs = uniqueN(og[dispOG]),
                     flg = ifelse(sum(dispOG)/.N > 0.05, "***", "")),
               by = "genome"]
    tab[,`:=`(lab = align_charLeft(genome),
              nGenes = align_charRight(nGenes),
              nOgs = align_charRight(nOgs))]
    with(tab, cat(sprintf(
      "\t...%s: %s genes in %s OGs hit > %s unique places %s\n",
      lab, nGenes, nOgs, max(maxOgPlaces), flg)))
    anyBad <- any(tab$flg == "***")
    if(anyBad)
      cat(strwrap(
        "NOTE! Genomes flagged *** have > 5% of genes in over-dispersed
      orthogroups. These are likely not great annotations, or the synteny run
      contains un-specified WGDs. Regardless, these should be examined
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
    isArrayRep <- globOG <- globHOG <-  og <- globHOG <-
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
    isArrayRep = FALSE, globOG = NA, globHOG = NA,
     noAnchor = FALSE, smallChr = FALSE, dispOG = FALSE)]

  ##############################################################################
  # 3. Add in necesary data columns if they do not exist
  # -- 3.1 orthofinder IDs
  bed <- add_ofID2bed(bed = bed, gsParam = gsParam)

  # -- 3.2 peptide lengths
  bed <- add_pepLen2bed(bed = bed, gsParam = gsParam)

  # -- 3.3 orthogroup IDs
  bed <- add_og2bed(bed = bed, gsParam = gsParam)

  # -- 3.4 hier orthogroup IDs
  if(useHOGs){
    bed <- add_hog2bed(bed = bed, gsParam = gsParam)
  }else{
    bed[,globHOG := NA]
  }

  if(useHOGs){
    bed[,og := globHOG]
  }else{
    bed[,og := globOG]
  }

  ##############################################################################
  # 4. Arrays and multi-placed orthogroups


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

  bed[,`:=`(n = NULL)]

  write_combBed(
    x = bed, filepath = combBedFile)
  return(bed)
}
