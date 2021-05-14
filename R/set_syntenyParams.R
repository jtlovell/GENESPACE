#' @title Calculate paramters for pairwise synteny search
#' @description
#' \code{set_syntenyParams} Calculate paramters for pairwise synteny search.
#' Other synteny functions require this as input.
#'
#' @param gsAnnot list of length 2, containing the genespace annotation paths
#' 'gff' and 'peptide' -- file path character vector with the locations of
#' the peptide and gff-like annotation files. Each element is named by the
#' associated genomeID in gsParam. Can be made manually as a list with the
#' 'peptide' entry as a named vector of paths to the parsed peptide fastas.
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run) and ploidy
#' (integer vector of genome ploidies named by genomeIDs).
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
#' @param onlyOg Logical of length 1 specifying if blast anchors
#' should only be formed from hits in the same orthogroups
#' @param nGapsHomeo see nGaps, but passed to homeolog  hits after masking
#' primary hits.
#' @param nGapsSecond see nGaps, but passed to secondary hits after masking
#' primary hits.
#' @param blkSizeHomeo see blkSize, but passed to the homeolog scan if
#' nSecondaryHits > 0.
#' @param blkSizeSecond see blkSize, but passed to the secondary scan if
#' nSecondaryHits > 0.
#' @param synBuffHomeo see syntenyBuffer. Applied only to synteny
#' construction of homeolog hits.
#' @param synBuffSecond see syntenyBuffer. Applied only to synteny
#' construction of secondary hits.
#' @param onlyOgHomeo see onlyOg. Applied only to
#' synteny construction of homeolog hits.
#' @param onlyOgSecond  see onlyOg. Applied only to
#' synteny construction of secondary hits.
#'
#' @details ...
#'
#' @return ...
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @importFrom Biostrings readAAStringSet
#' @import data.table
#' @export
set_syntenyParams <- function(
  gsParam,
  gsAnnot,
  nGaps = 5,
  blkSize = 5,
  nSecondHits = 0,
  synBuff = 100,
  onlyOg = TRUE,
  nGapsHomeo = nGaps,
  nGapsSecond = nGaps,
  blkSizeHomeo = blkSize,
  blkSizeSecond = blkSize,
  synBuffHomeo = synBuff,
  synBuffSecond = synBuff,
  onlyOgHomeo = TRUE,
  onlyOgSecond = FALSE)
{

  target <- genome1 <- genome2 <- ord2 <- ofID1 <- ofID2 <- ord1 <- chr1 <- chr2 <- NULL
  nSecondHits1 <- nSecondHits2 <- nGenes1 <- nGenes2 <- tiebreak <- query <- NULL
  ploidy1 <- ploidy2 <- nhits2 <- nhits1 <- mirrorBlast <- runBlast <- u <- NULL

  genomeIDs <- gsParam$genomeIDs[!gsParam$genomeIDs %in% gsParam$outgroup]
  ploidy <- gsParam$ploidy[genomeIDs]
  blastDir <- gsParam$blast
  pepFiles <- gsAnnot$peptide[genomeIDs]

  ##############################################################################
  # -- 0. Argument checking
  onlyOgSecond <- check_logicalArg(onlyOgSecond)
  onlyOgHomeo <- check_logicalArg(onlyOgHomeo)
  onlyOg <- check_logicalArg(onlyOg)

  ##############################################################################
  # -- 1. make the database of unique combinations of genomes
  cmb <- data.table(genome1 = names(ploidy), genome2 = names(ploidy))
  cmb <- cmb[,CJ(genome1, genome2)]
  cmb[,u := paste(genome1, genome2)]
  cmb[,ploidy1 := ploidy[genome1]]
  cmb[,ploidy2 := ploidy[genome2]]

  speciesIDs <- read_orthofinderSpeciesIDs(blastDir)
  cmb[,ofID1 := speciesIDs[genome1]]
  cmb[,ofID2 := speciesIDs[genome2]]

  ##############################################################################
  # -- 2. Choose query and target genomes, based on the total number of genes
  nGenes <- sapply(pepFiles, function(x) length(readAAStringSet(x)))
  cmb[,nGenes1 := nGenes[genome1]]
  cmb[,nGenes2 := nGenes[genome2]]
  cmb[,tiebreak := as.numeric(factor(genome1, levels = genomeIDs)) - as.numeric(factor(genome2, levels = genomeIDs))]
  cmb[,query := ifelse(nGenes1 > nGenes2,
                       genome1,
                       ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                              genome1,
                              genome2)),
      by = "u"]
  cmb[,target := ifelse(nGenes1 > nGenes2,
                        genome2,
                        ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                               genome2,
                               genome1)),
      by = "u"]

  cmb[,runBlast := genome1 == query]
  cmb[,mirrorBlast := genome1 == target & genome1 != query]
  cmb[,tiebreak := NULL]

  ##############################################################################
  # -- 3. Choose the number of hits, based on ploidy
  cmb[,nhits1 := ploidy2]
  cmb[,nhits2 := ploidy1]
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
  # -- 4. Check, then add all the parameters into the output data.table.
  nGaps <- as.integer(nGaps)[1]
  if(is.null(nGaps) || is.na(nGaps) || length(nGaps) == 0)
    stop("cannot coerce nGaps to integer\n")
  blkSize <- as.integer(blkSize)[1]
  if(is.null(blkSize) || is.na(blkSize) || length(blkSize) == 0)
    stop("cannot coerce blkSize to integer\n")
  nSecondHits <- as.integer(nSecondHits)[1]
  if(is.null(nSecondHits) || is.na(nSecondHits) || length(nSecondHits) == 0)
    stop("cannot coerce nSecondHits to integer\n")
  synBuff <- as.integer(synBuff)[1]
  if(is.null(synBuff) || is.na(synBuff) || length(synBuff) == 0)
    stop("cannot coerce synBuff to integer\n")
  onlyOg <- as.logical(onlyOg)[1]
  if(is.null(onlyOg) || is.na(onlyOg) || length(onlyOg) == 0)
    stop("cannot coerce onlyOg to logical\n")

  nGapsHomeo <- as.integer(nGapsHomeo)[1]
  if(is.null(nGapsHomeo) || length(nGapsHomeo) == 0)
    stop("cannot coerce nGapsHomeo to integer\n")
  nGapsSecond <- as.integer(nGapsSecond)[1]
  if(is.null(nGapsSecond) || length(nGapsSecond) == 0)
    stop("cannot coerce nGapsSecond to integer\n")

  blkSizeHomeo <- as.integer(blkSizeHomeo)[1]
  if(is.null(blkSizeHomeo) || length(blkSizeHomeo) == 0)
    stop("cannot coerce blkSizeHomeo to integer\n")
  blkSizeSecond <- as.integer(blkSizeSecond)[1]
  if(is.null(blkSizeSecond) || length(blkSizeSecond) == 0)
    stop("cannot coerce blkSizeSecond to integer\n")

  synBuffHomeo <- as.integer(synBuffHomeo)[1]
  if(is.null(synBuffHomeo) || is.na(synBuffHomeo) || length(synBuffHomeo) == 0)
    stop("cannot coerce synBuffHomeo to integer\n")
  synBuffSecond <- as.integer(synBuffSecond)[1]
  if(is.null(synBuffSecond) || is.na(synBuffSecond) || length(synBuffSecond) == 0)
    stop("cannot coerce synBuffSecond to integer\n")

  onlyOgHomeo <- as.logical(onlyOgHomeo)[1]
  if(is.null(onlyOgHomeo) || length(onlyOgHomeo) == 0)
    stop("cannot coerce onlyOgHomeo to logical\n")
  onlyOgSecond <- as.logical(onlyOgSecond)[1]
  if(is.null(onlyOgSecond) || length(onlyOgSecond) == 0)
    stop("cannot coerce onlyOgSecond to logical\n")

  params <- c("nGaps", "blkSize", "nSecondHits", "synBuff", "onlyOg",
              "nGapsHomeo", "nGapsSecond" ,
              "blkSizeHomeo", "blkSizeSecond",
              "synBuffHomeo", "synBuffSecond",
              "onlyOgHomeo", "onlyOgSecond")
  for(i in params)
    cmb[[i]] <- get(i)

  cmb <- cmb[with(cmb, order(-(nGenes1 + nGenes2))),]

  return(cmb)
}
