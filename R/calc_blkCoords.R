#' @title Build custom syntenic blocks based on regions or reference genomes
#'
#' @description
#' \code{recall_synBlks} Faster re-calling of syntenic block breakpoints to let
#' users alter the performance of the riparian plotter. NOTE: this function will
#' only re-call syntenic blocks that are the same size or smaller than the
#' parameters used in synteny. If you want larger syntenic blocks (esp. with
#' larger blkRadius), you will need to re-run the synteny function of GENESPACE.
#'
#' @name calc_blkCoords
#'
#' @param genomeIDs character vector of length > 1, matching length
#' of speciesIDs, versions and ploidy. Specifies the name to assign
#' to each genome. This vector must be unique and can be any string
#' that begins with a letter (a-z, A-Z) and is alphanumeric. '.' and '_' are
#' allowed as long as they are not the first character.
#' @param synHitsDir file path to the directory containing the syntenic hits
#' @param nCores integer specifying the number of parallel processes to run
#' @param niter integer specifying the number of iterations to re-call syntenic
#' blocks. Must be > 0. 1 should be fine, but 2 may resolve some issues in some
#' cases. It is not clear that higher numbers will result in better calls.
#' @param verbose logical, should updates be printed to the console?
#' @param nGaps integer of length 1, specifying the -m param to mcscanx
#' for the primary MCScanX run. This acts on the results from the initial
#' MCScanX run.
#' @param blkSize integer of length 1, specifying the -s param to mcscanx
#' @param blkRadius integer of length 1, specifying the search radius in 2d
#' clustering to assign hits to the same block. This is a sensitive parameter
#' as smaller values will result in more blocks, gaps and SV. Typically using
#' 2x or greater blkSize is fine.
#'
#' @details There are six behavoiors of this function, that are determined by
#' the parameters given.
#'
#' (1) The default performance (only providing a gsParam object) simply
#' re-calculates block coordinates. If the blkSize, blkRadius and nGaps
#' parameters are all NULL, this will simply return the coordinates from hits
#' in the syntenicBlast directory.
#'
#' (2) If given a reference genome and a gsParam object, the syntenic blocks
#' are re-calculated, phased against the reference genome.
#'
#' (3) If given only a gsParam object and nGaps/blkSize/blkRadius parameters
#' that do not match the gsParam object, will re-call syntenic blocks WITHIN
#' existing blocks.
#'
#' (4) If given a reference genome and (3), re-calls and phases syntenic blocks
#'
#' (5) If given a region bed data.table/data.frame, will extract each region,
#' phase against the reference for each region.
#'
#' (6) If given nGaps/blkSize/blkRadius parameters that do not match the gsParam
#' object and (5), re-calculates block coordinates for each region.

#' @return A data.table containing block coordinates
#' @examples
#' \dontrun{
#' # coming soon
#' }


#' @title calculate syntenic block coordinates
#' @description
#' \code{calc_blkCoords} from a hits object, determine block coordinates,
#' orientation and membership
#' @rdname calc_blkCoords
#' @import data.table
#' @importFrom parallel mclapply
#' @export
# recall_blks <- function(gsParam = NULL,
#                         regBed = NULL,
#                         path2combBed = NULL,
#                         synHitsDir = NULL,
#                         pangenomeDir = NULL,
#                         genomeIDs = NULL,
#                         refGenome = NULL,
#                         blkRadius = NULL,
#                         blkSize = NULL,
#                         nGaps = NULL,
#                         nCores = NULL,
#                         niter = 2,
#                         onlyAnchors = TRUE,
#                         maxInterpOrdDiff = 100,
#                         verbose = TRUE){
#
#   ##############################################################################
#   ##############################################################################
#   ##############################################################################
#   # 1. check general parameters
#
#   ##############################################################################
#   # -- 1.1 regBed
#   if(!is.null(regBed)){
#     if(!is.data.frame(regBed)){
#       stop("region bed object provided but is not a data.table or data.frame\n")
#     }else{
#       bedCols <- c("genome", "chr", "start", "end", "id", "color")
#       if(!all(bedCols %in% colnames(regBed))){
#         stop(sprintf(
#           "region bed object provided but does not contain %s columns\n",
#           paste(bedCols, collapse = ", ")))
#       }
#     }
#   }
#
#   ##############################################################################
#   # -- 1.2 refGenome, gsParam, pangenomeDir and synHitsDir
#   if(!is.null(refGenome) && !is.null(regBed))
#     stop("can only recall blocks for a reference genome OR a set of intervals, not both\n")
#   if(is.null(gsParam))
#     if(is.null(synHitsDir) || !dir.exists(synHitsDir))
#       stop("either the genespace parameters (gsParam) or a valid synHitsDir must be specified\n")
#   if(!is.null(gsParam))
#     synHitsDir <- gsParam$paths$syntenicHits
#   if(!dir.exists(synHitsDir))
#     stop("there is a problem with your gsParam object - can't find the path to syntenic hits\n")
#   if(is.null(pangenomeDir) && !is.null(gsParam))
#     pangenomeDir <- gsParam$paths$pangenome
#
#   ##############################################################################
#   # -- 1.3 blkSize
#   if(is.null(blkSize)){
#     blkSize <- gsParam$params$blkSize
#   }else{
#     blkSize <- check_integer(blkSize, onlySingleValue = TRUE)
#     if(!is.integer(blkSize))
#       stop("blkSize parameter is malformed\n")
#   }
#
#   ##############################################################################
#   # -- 1.4 blkRadius
#   if(is.null(blkRadius)){
#     blkRadius <- gsParam$params$blkRadius
#   }else{
#     blkRadius <- check_integer(blkRadius, onlySingleValue = TRUE)
#     if(!is.integer(blkRadius))
#       stop("blkRadius parameter is malformed\n")
#   }
#
#   ##############################################################################
#   # -- 1.5 nGaps
#   if(is.null(nGaps)){
#     nGaps <- gsParam$params$nGaps
#   }else{
#     nGaps <- check_integer(nGaps, onlySingleValue = TRUE)
#     if(!is.integer(nGaps))
#       stop("nGaps parameter is malformed\n")
#   }
#
#   ##############################################################################
#   # -- 1.6 nCores
#   if(is.null(nCores)){
#     nCores <- gsParam$params$nCores
#   }else{
#     nCores <- check_integer(nCores, onlySingleValue = TRUE)
#     if(!is.integer(nCores))
#       stop("nCores parameter is malformed\n")
#   }
#
#   ##############################################################################
#   # -- 1.7 niter
#   niter <- check_integer(niter, onlySingleValue = TRUE)
#   if(!is.integer(niter))
#     stop("niter parameter is malformed\n")
#
#   ##############################################################################
#   # -- 1.8 verbose
#   if(!is.logical(verbose))
#     verbose <- TRUE
#
#   ##############################################################################
#   # -- 1.9 combBed file
#   if(is.null(path2combBed))
#     if(is.null(gsParam))
#       stop("either gsParam or path2combBed needs to be specified")
#   if(is.null(path2combBed))
#     path2combBed <- file.path(gsParam$paths$results, "combBed.txt")
#   if(!file.exists(path2combBed))
#     stop("could not find a valid path to the combBed.txt file.")
#
#   if(!is.null(path2combBed)){
#     if(!file.exists(path2combBed)){
#       bed <- NULL
#       path2combBed <- NULL
#     }else{
#       bed <- read_combBed(path2combBed)
#     }
#   }else{
#     bed <- NULL
#     path2combBed <- NULL
#   }
#
#   ##############################################################################
#   ##############################################################################
#   ##############################################################################
#   # 2. Determine function behavior
#
#   ##############################################################################
#   # -- 2.1 get list of syntenic hits files
#   fs <- list.files(
#     path = synHitsDir, pattern = "synBlast.txt.gz", full.names = T)
#
#   ##############################################################################
#   # -- 2.2 get list of unique genomeIDs available
#   names(fs) <- align_charLeft(
#     gsub(".synBlast.txt.gz", "", gsub("_vs_", " v. ", basename(fs))))
#   gs <- as.data.table(tstrsplit(
#     gsub(".synBlast.txt.gz", "", basename(fs)), "_vs_"))
#
#   ##############################################################################
#   # -- 2.3 subset genomes if genomeIDs provided
#   if(!is.null(genomeIDs)){
#     gs[,use := V1 %in% genomeIDs & V2 %in% genomeIDs]
#     if(sum(gs$use) == 0)
#       stop("those genomeIDs do not appear in the synHits directory\n")
#     fs <- fs[gs$use]
#     gs <- subset(gs, use)
#   }else{
#     genomeIDs <- unique(unlist(gs[,1:2]))
#   }
#
#   ##############################################################################
#   # -- 2.4 check that refGenome is among the genomeIDs
#   if(!is.null(refGenome)){
#     refGenome <- check_character(refGenome, onlySingleValue = T)
#     if(is.na(refGenome))
#       refGenome <- NULL
#   }
#   if(!is.null(refGenome)){
#     if(!refGenome %in% unlist(gs))
#       stop("provided refgenome is not among genomes available from synteny\n")
#   }
#
#   ##############################################################################
#   # -- 2.5 check regBed, make sure that all regions contain genes
#   if(!is.null(regBed)){
#     if(is.null(bed))
#       stop("there is a problem with the combBed file specification\n")
#     regBed <- data.table(regBed)
#     setkey(regBed, genome, chr, start, end)
#     setkey(bed, genome, chr, start, end)
#     bedInReg <- foverlaps(regBed, bed)
#     cnt <- bedInReg[,list(nGenes = uniqueN(ofID, na.rm = T)), by = "i.id"]
#     if(all(cnt$nGenes == 0))
#       stop("there is a problem with your bedReg object - none of the regions have any genes in the combBed.txt file\n")
#     if(any(cnt$nGenes == 0))
#       warning(sprintf("region(s) %s have no genes and will not be shown",
#                       paste(cnt$i.id[cnt$nGenes == 0], collapse = ", ")))
#     bedInReg <- subset(bedInReg, !i.id %in% cnt$i.id[cnt$nGenes == 0])
#     cnt <- subset(cnt, nGenes > 0)
#     regBed <- subset(regBed, id %in% cnt$i.id)
#     if(any(cnt$nGenes < blkSize)){
#       warning(sprintf("region(s) %s have fewer genes than supplied block size. Reducing min block size to %s",
#                       paste(cnt$i.id[cnt$nGenes < blkSize], collapse = ", "), min(cnt$nGenes)))
#       blkSize <- min(cnt$nGenes)
#     }
#   }
#
#   ##############################################################################
#   # -- 2.6 check if the synteny parameters are different than the gsParam
#   if(is.null(gsParam)){
#     sameParams <- TRUE
#   }else{
#     s <- gsParam$params$blkSize
#     r <- gsParam$params$blkRadius
#     g <- gsParam$params$nGaps
#     if(s == blkSize && r == blkRadius && g == nGaps){
#       sameParams <- TRUE
#     }else{
#       sameParams <- FALSE
#     }
#   }
#
#   ##############################################################################
#   # -- 2.7 check for interpolated syntenic position files
#   if(!is.null(pangenomeDir) && !is.null(refGenome) || !is.null(regBed)){
#     if(is.null(refGenome)){
#       rg <- regBed$genome
#     }else{
#       rg <- refGenome
#     }
#     if(!dir.exists(pangenomeDir))
#       stop("refGenome or regBed provided, but path to pangenome directory is not valid\n")
#     intf <- list.files(pangenomeDir , pattern = "integratedSynPos.txt$", full.names = T)
#     strs <- sprintf("_vs_%s.integratedSynPos.txt", rg)
#     hasit <- all(sapply(strs, function(x) any(grepl(x, basename(intf)))))
#     if(!hasit)
#       stop("there is a problem with the integratedSynPos.txt files - the reference genome(s) are not in there\n")
#   }
#
#   ##############################################################################
#   # -- 2.8 finally, use above info to determine function behavior
#   if(sameParams && is.null(refGenome) && is.null(regBed))
#     runType <- 1
#   if(sameParams && !is.null(refGenome) && is.null(regBed))
#     runType <- 2
#   if(!sameParams && is.null(refGenome) && is.null(regBed))
#     runType <- 3
#   if(!sameParams && !is.null(refGenome) && is.null(regBed))
#     runType <- 4
#   if(sameParams && is.null(refGenome) && !is.null(regBed))
#     runType <- 5
#   if(!sameParams && is.null(refGenome) && !is.null(regBed))
#     runType <- 6
#
#   ##############################################################################
#   # -- 2.9 split up the synHits files into chunks
#   if(verbose)
#     cat(sprintf("Found %s hits files across %s genomes",
#                 length(fs), uniqueN(unlist(gs[,1:2]))))
#   fspl <- split(fs, rep(1:length(fs), each = nCores)[1:length(fs)])
#
#
#   ##############################################################################
#   ##############################################################################
#   ##############################################################################
#   # 3 Re-call the blocks
#
#   ##############################################################################
#   # 3.1 -- runtype 1 (just the same as default)
#   if(runType == 1){
#     outblk <- rbindlist(lapply(fspl, function(chnki){
#       outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#         d <- subset(read_synHits(chnki[i]), isAnchor)
#         b <- calc_blkCoords(d, mirror = T)
#
#         u <- with(d[1,], paste(genome1, genome2))
#         b[,`:=`(refChr = NA, regID = NA, lab = i, isMirror = !paste(genome1, genome2) %in% u)]
#
#         return(b)
#       }))
#
#       outSv <- subset(outchnk, !isMirror)[,list(
#         nSVs = uniqueN(blkID) - uniqueN(paste(chr1, chr2))),
#         by = "lab"]
#       if(verbose)
#         with(outSv, cat(sprintf("%s: n. SVs = %s", lab, nSVs), sep = "\n"))
#
#       return(outchnk)
#     }))
#   }
#
#   ##############################################################################
#   # 3.2 -- runtype 2 (only a reference genome)
#   if(runType == 2){
#     intf <- file.path(pangenomeDir, sprintf(
#       "%s_vs_%s.integratedSynPos.txt", genomeIDs, refGenome))
#     intp <- rbindlist(mclapply(intf, mc.cores = nCores, read_intSynPos))
#     if(onlyAnchors)
#       intp <- subset(intp, isAnchor)
#     intp1 <- with(intp, data.table(
#       ofID1 = ofID, refChr1 = interpChr, refOrd1 = interpOrd))
#     intp2 <- with(intp, data.table(
#       ofID2 = ofID, refChr2 = interpChr, refOrd2 = interpOrd))
#     outblk <- rbindlist(lapply(fspl, function(chnki){
#       outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#         d <- subset(read_synHits(chnki[i]), isAnchor)
#         d <- merge(intp1, merge(intp2, d, by = "ofID2", allow.cartesian = T),
#                    by = "ofID1", allow.cartesian = T)
#         d <- subset(d, refChr1 == refChr2 & abs(refOrd1 - refOrd2) < maxInterpOrdDiff)
#         d[,blkID := sprintf("%sXXXrefGenomeChrXXX%s", blkID, refChr1)]
#         b <- calc_blkCoords(d, mirror = T)
#         b[,c("blkID", "refChr") := tstrsplit(blkID, "XXXrefGenomeChrXXX")]
#         u <- with(d[1,], paste(genome1, genome2))
#         b[,`:=`(regID = NA, lab = i, isMirror = !paste(genome1, genome2) %in% u)]
#         return(b)
#       }))
#
#       outSv <- subset(outchnk, !isMirror)[,list(
#         nSVs = uniqueN(paste(blkID, refChr)) - uniqueN(paste(chr1, chr2))),
#         by = "lab"]
#       if(verbose)
#         with(outSv, cat(sprintf("%s: n. SVs = %s", lab, nSVs), sep = "\n"))
#
#       return(outchnk)
#     }))
#   }
#
#   ##############################################################################
#   # 3.3 -- runtype 3, new params .. re-call blks
#   if(runType == 3){
#     outblk <- rbindlist(lapply(fspl, function(chnki){
#       outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#         d <- subset(read_synHits(chnki[i]), isAnchor)
#
#         for(j in 1:niter){
#           d[,`:=`(ordt1 = frank(ord1, ties.method = "dense"),
#                   ordt2 = frank(ord2, ties.method = "dense")),
#             by = c("chr1", "chr2", "blkID")]
#           d[,tmp := dbscan(frNN(
#             x = cbind(ordt1, ordt2),
#             eps = blkRadius),
#             minPts = blkSize)$cluster,
#             by = c("chr1", "chr2", "blkID")]
#           d <- subset(d, tmp != 0)
#         }
#         d[,blkID := as.numeric(as.factor(paste(chr1, chr2, blkID, tmp)))]
#         d[,blkID := paste(genome1, genome2,blkID)]
#
#         b <- calc_blkCoords(d, mirror = T)
#
#         u <- with(d[1,], paste(genome1, genome2))
#         b[,`:=`(refChr = NA, regID = NA, lab = i, isMirror = !paste(genome1, genome2) %in% u)]
#
#         return(b)
#       }))
#
#       outSv <- subset(outchnk, !isMirror)[,list(
#         nSVs = uniqueN(blkID) - uniqueN(paste(chr1, chr2))),
#         by = "lab"]
#       if(verbose)
#         with(outSv, cat(sprintf("%s: n. SVs = %s", lab, nSVs), sep = "\n"))
#
#       return(outchnk)
#     }))
#   }
#
#   ##############################################################################
#   # 3.4 -- runtype 4
#   if(runType == 4){
#     intf <- file.path(pangenomeDir, sprintf(
#       "%s_vs_%s.integratedSynPos.txt", genomeIDs, refGenome))
#     intp <- rbindlist(mclapply(intf, mc.cores = nCores, read_intSynPos))
#     if(onlyAnchors)
#       intp <- subset(intp, isAnchor)
#     intp1 <- with(intp, data.table(
#       ofID1 = ofID, refChr1 = interpChr, refOrd1 = interpOrd))
#     intp2 <- with(intp, data.table(
#       ofID2 = ofID, refChr2 = interpChr, refOrd2 = interpOrd))
#     outblk <- rbindlist(lapply(fspl, function(chnki){
#       outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#         d <- subset(read_synHits(chnki[i]), isAnchor)
#         d <- merge(intp1, merge(intp2, d, by = "ofID2", allow.cartesian = T),
#                    by = "ofID1", allow.cartesian = T)
#         d <- subset(d, refChr1 == refChr2 & abs(refOrd1 - refOrd2) < maxInterpOrdDiff)
#
#         for(j in 1:niter){
#           d[,`:=`(ordt1 = frank(ord1, ties.method = "dense"),
#                   ordt2 = frank(ord2, ties.method = "dense")),
#             by = c("chr1", "chr2", "blkID")]
#           d[,tmp := dbscan(frNN(
#             x = cbind(ordt1, ordt2),
#             eps = blkRadius),
#             minPts = blkSize)$cluster,
#             by = c("chr1", "chr2", "blkID")]
#           d <- subset(d, tmp != 0)
#         }
#         d[,blkID := as.numeric(as.factor(paste(chr1, chr2, blkID, tmp)))]
#         d[,blkID := paste(genome1, genome2,blkID)]
#
#         d[,blkID := sprintf("%sXXXrefGenomeChrXXX%s", blkID, refChr1)]
#         b <- calc_blkCoords(d, mirror = T)
#         b[,c("blkID", "refChr") := tstrsplit(blkID, "XXXrefGenomeChrXXX")]
#         u <- with(d[1,], paste(genome1, genome2))
#         b[,`:=`(regID = NA, lab = i, isMirror = !paste(genome1, genome2) %in% u)]
#         return(b)
#       }))
#
#       outSv <- subset(outchnk, !isMirror)[,list(
#         nSVs = uniqueN(paste(blkID, refChr)) - uniqueN(paste(chr1, chr2))),
#         by = "lab"]
#       if(verbose)
#         with(outSv, cat(sprintf("%s: n. SVs = %s", lab, nSVs), sep = "\n"))
#
#       return(outchnk)
#     }))
#   }
#
#   ##############################################################################
#   # 3.5 -- runtype 5
#   if(runType == 5){
#     ogv <- bedInReg$og; names(ogv) <- bedInReg$ofID
#     ogl <- split(bedInReg$og, bedInReg$i.id)
#     outblk <- rbindlist(lapply(fspl, function(chnki){
#       outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#         d <- subset(read_synHits(chnki[i]), isAnchor)
#         d[,`:=`(og1 = ogv[ofID1], og2 = ogv[ofID2])]
#         d <- subset(d, !is.na(og1) | !is.na(og2))
#         d <- subset(d, !(!is.na(og1) & !is.na(og2) & og1 != og2))
#
#         d <- rbindlist(lapply(names(ogl), function(j){
#           tmp <- subset(d, og1 %in% ogl[[j]] | og2 %in% ogl[[2]])
#           tmp[,regID := j]
#           return(tmp)
#         }))
#
#         d[,blkID := sprintf("%sXXXbedRegIDXXX%s", blkID, regID)]
#         b <- calc_blkCoords(d, mirror = T)
#         b[,c("blkID", "regID") := tstrsplit(blkID, "XXXbedRegIDXXX")]
#         u <- with(d[1,], paste(genome1, genome2))
#         b[,`:=`(refChr = NA)]
#       }))
#       return(outchnk)
#     }))
#
#
#     gids <- unique(regBed$genome)
#     refOut <- rbindlist(lapply(gids, function(rg){
#       refbedInReg <- subset(bedInReg, genome == rg)
#       refRegBed <- subset(regBed, genome == rg)
#       refOrdReg <- refbedInReg[,list(start = min(ord), end = max(ord)),
#                                by = c("genome", "chr", "i.id", "color")]
#       setkey(refOrdReg, genome, chr, start, end)
#       intf <- file.path(pangenomeDir, sprintf(
#         "%s_vs_%s.integratedSynPos.txt", genomeIDs, rg))
#
#       intp <- rbindlist(mclapply(intf, mc.cores = nCores, function(i){
#         x <- read_intSynPos(i)
#         regord <-
#         if(onlyAnchors)
#           x <- subset(x, isAnchor)
#         x[,`:=`(start = interpOrd, end = interpOrd, chr = interpChr, genome = interpGenome)]
#         setkey(x, genome, chr, start, end)
#         fo <- foverlaps(refOrdReg, x)
#         out <- with(fo, data.table(ofID = ofID, regID = i.id, interpOrd = interpOrd))
#         return(out)
#       }))
#       return(intp)
#     }))
#
#     intp1 <- with(refOut, data.table(
#       ofID1 = ofID, regID1 = regID, refOrd1 = interpOrd))
#     intp2 <- with(refOut, data.table(
#       ofID2 = ofID, regID2 = regID, refOrd2 = interpOrd))
#
#     outblk <- rbindlist(lapply(fspl, function(chnki){
#       outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#         d <- subset(read_synHits(chnki[i]), isAnchor)
#         d <- merge(intp1, merge(intp2, d, by = "ofID2", allow.cartesian = T),
#                    by = "ofID1", allow.cartesian = T)
#
#         return(b)
#       }))
#
#       return(outchnk)
#     }))
#
#
#
#   }
#
#   ##############################################################################
#   # 3.6 -- runtype 6
#   if(runType == 6){
#
#   }
#
#   if(!is.null(refGenome) && )
#
#
#
#   phaseIt <- !is.null(refGenome)
#   if(phaseIt){
#     bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
#     refg <- with(subset(bed, genome == refGenome), data.table(
#       refChr1 = chr, refChr2 = chr, og = og))
#     altg <- with(bed, data.table(ofID1 = ofID, ofID2 = ofID, og = og))
#     refChrDict <- merge(refg, altg, by = "og", allow.cartesian = TRUE)
#   }
#
#   outblk <- rbindlist(lapply(fspl, function(chnki){
#     outchnk <- rbindlist(mclapply(names(chnki), mc.cores = nCores, function(i){
#       d <- subset(read_synHits(chnki[i]), isAnchor)
#       if(phaseIt){
#         d <- merge(
#           d, refChrDict[,c("refChr1", "ofID1")],
#           by = "ofID1", allow.cartesian = T)
#         d <- merge(
#           d, refChrDict[,c("refChr2", "ofID2")],
#           by = "ofID2", allow.cartesian = T)
#         d <- subset(d, refChr1 == refChr2)
#         d[,`:=`(blkID = paste(blkID, refChr1), refChr = refChr1,
#                  refChr1 = NULL, refChr2 = NULL)]
#       }
#       for(j in 1:niter){
#         d[,`:=`(ordt1 = frank(ord1, ties.method = "dense"),
#                 ordt2 = frank(ord2, ties.method = "dense")),
#           by = c("chr1", "chr2", "blkID")]
#         d[,tmp := dbscan(frNN(
#           x = cbind(ordt1, ordt2),
#           eps = blkRadius),
#           minPts = blkSize)$cluster,
#           by = c("chr1", "chr2", "blkID")]
#         d <- subset(d, tmp != 0)
#       }
#       d[,blkID := as.numeric(as.factor(paste(chr1, chr2, blkID, tmp)))]
#       d[,blkID := paste(genome1, genome2,blkID)]
#       b <- calc_blkCoords(d, mirror = T)
#       u <- with(d[1,], paste(genome1, genome2))
#       b[,`:=`(lab = i, isMirror = !paste(genome1, genome2) %in% u)]
#       return(b)
#     }))
#     outSv <- subset(outchnk, !isMirror)[,list(
#       nSVs = uniqueN(blkID) - uniqueN(paste(chr1, chr2))),
#       by = "lab"]
#     if(verbose)
#       with(outSv, cat(sprintf("%s: n. SVs = %s", lab, nSVs), sep = "\n"))
#     return(outchnk)
#   }))
#
#
#   if(is.null(outFile)){
#     return(outblk)
#   }else{
#     fwrite(
#       outblk, file = outFile, sep = "\t", quote = FALSE, showProgress = FALSE)
#   }
# }


#' @title calculate syntenic block coordinates
#' @description
#' \code{calc_blkCoords} from a hits object, determine block coordinates,
#' orientation and membership
#' @rdname utils
#' @import data.table
#' @importFrom stats cor
#' @export
calc_blkCoords <- function(hits, mirror = FALSE){
  setDTthreads(1)

  # -- get the columns and complete observations for these
  hcols <- c("blkID", "start1", "start2", "end1", "end2", "ord1", "ord2",
             "chr1", "chr2", "genome1", "genome2", "ofID1", "ofID2")
  bhits <- subset(hits, complete.cases(hits[,hcols, with = F]))

  if(mirror){
    tmp <- data.table(bhits)
    setnames(tmp, gsub("2$", "3", colnames(tmp)))
    setnames(tmp, gsub("1$", "2", colnames(tmp)))
    setnames(tmp, gsub("3$", "1", colnames(tmp)))
    bhits <- rbind(bhits, tmp[,colnames(bhits), with = F])
    bhits <- subset(bhits, !duplicated(paste(ofID1, ofID2, blkID)))
  }

  # -- get the genome1 coordinates
  ofID1 <- start1 <- end1 <- ofID1 <- ord1 <- blkID <- NULL
  setkey(bhits, ord1)
  blks1 <- bhits[,list(
    startBp1 = min(start1), endBp1 = max(end1),
    startOrd1 = min(ord1), endOrd1 = max(ord1),
    firstGene1 = first(ofID1), lastGene1 = last(ofID1),
    nHits1 = uniqueN(ofID1)),
    by = c("blkID", "genome1","genome2", "chr1", "chr2")]

  # -- get the genome2 coordinates
  ofID2 <- start2 <- end2 <- ofID2 <- ord2 <- NULL
  setkey(bhits, ord2)
  blks2 <- bhits[,list(
    minBp2 = min(start2), maxBp2 = max(end2),
    minOrd2 = min(ord2), maxOrd2 = max(ord2),
    minGene2 = first(ofID2), maxGene2 = last(ofID2),
    nHits2 = uniqueN(ofID2),
    orient = ifelse(length(ord1) <= 1, "+",
                    ifelse(cor(jitter(ord1),
                               jitter(ord2)) > 0,"+", "-"))),
    by = c("blkID", "genome1","genome2", "chr1", "chr2")]

  # -- merge the two coordinates
  blks <- merge(blks1, blks2, by = c("genome1","genome2","chr1","chr2","blkID"))

  # -- fix the coordinates for inverted blocks
  orient <- NULL
  bgfor <- subset(blks, orient == "+")
  bgrev <- subset(blks, orient == "-")

  maxBp2 <- minBp2 <- maxOrd2 <- minOrd2 <- maxGene2 <- minGene2 <- NULL
  bgrev[,`:=`(startBp2 = maxBp2, endBp2 = minBp2,
              startOrd2 = maxOrd2, endOrd2 = minOrd2,
              firstGene2 = maxGene2, lastGene2 = minGene2)]
  bgfor[,`:=`(startBp2 = minBp2, endBp2 = maxBp2,
              startOrd2 = minOrd2, endOrd2 = maxOrd2,
              firstGene2 = minGene2, lastGene2 = maxGene2)]
  blks <- rbind(bgfor, bgrev)
  return(blks)
}
