#' @title Combine and annotate blast files
#' @description
#' \code{annotate_blast} Combines all orthogroup and exclustion data for each
#' gene
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#'
#' @import data.table
#' @importFrom grDevices pdf dev.off rgb
#' @importFrom parallel mclapply
#' @export
annotate_blast <- function(gsParam){

  bitScore <- ofID1 <- ofID2 <- og1 <- og2 <- og <- noAnchor1 <- ns1 <- ns2 <-
    noAnchor2 <- noAnchor <- rnd1 <- rnd2 <- ngene1 <- ngene2 <- sameOG <-
    ord1 <- ord2 <- ancOrd1 <- ancOrd2 <- n <- queryBlast <- targetBlast <-
    wt <- size1 <- size2 <- chunk <- lab <- query <- target <- NULL

  blNames <- c(
    "ofID1", "ofID2", "pid", "length", "mismatches","gapopenings",
    "queryStart", "queryEnd", "subjectStart", "subjectEnd", "Evalue", "bitScore")
  blNamesR <- c(
    "ofID2", "ofID1", "pid", "length", "mismatches","gapopenings",
    "subjectStart", "subjectEnd", "queryStart", "queryEnd", "Evalue", "bitScore")
  bdNames <- c(
    "chr", "start", "end", "id", "ofID", "ord", "genome", "og", "noAnchor", "isArrayRep")

  ##############################################################################
  # 1. Build the blast file metadata
  if(!"synteny" %in% names(gsParam))
    stop("must run set_syntenyParams prior to annotate_blast\n")

  # -- 1.1 read in the combined bed
  bedFile <- gsParam$synteny$combBed
  if(!file.exists(bedFile))
    stop("could not find combined bedfile in", bedFile)

  bed <- read_combBed(bedFile)
  nCores <- gsParam$params$nCores

  # -- 1.2 order by total size
  synMd <- data.table(gsParam$synteny$blast)
  synMd[,`:=`(size1 = ifelse(!is.na(queryBlast),
                             file.size(queryBlast), file.size(targetBlast)),
              size2 = ifelse(!is.na(targetBlast),
                             file.size(targetBlast), file.size(queryBlast)))]
  synMd[,wt := (size1 / 1e9) + (size2 / 1e9)]
  setorder(synMd, -wt)

  # -- 1.3 split into chunks
  synMd[,chunk := rep(1:.N, each = nCores)[1:.N]]
  synMd[,`:=`( wt = NULL, size1 = NULL, size2 = NULL)]
  synMd[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]
  synMdSpl <- split(synMd, by = "chunk")

  ##############################################################################
  # 2. Process each set
  # -- 2.1 set vectors of blast names

  # -- 2.2 make two bed files for the merge with the blast
  bed1 <- data.table(bed[,bdNames, with = F])
  setnames(bed1, paste0(names(bed1), "1"))
  bed2 <- data.table(bed[,bdNames, with = F])
  setnames(bed2, paste0(names(bed2), "2"))

  # ... for each chunk of nCores size
  outmd <- rbindlist(lapply(1:length(synMdSpl), function(chnki){

    if(nCores > 1)
      cat(sprintf(
        "\t# Chunk %s / %s (%s) ... \n",
        chnki, max(synMd$chunk), format(Sys.time(), "%X")))

    chnk <- data.table(synMdSpl[[chnki]])

    # ... for each row in each chunk
    out <- rbindlist(mclapply(1:nrow(chnk), mc.cores = nCores, function(i){
      # -- 2.3 read in the blast file(s)
      x <- chnk[i,]
      if(!is.na(x$queryBlast)){
        blf <- fread(
          x$queryBlast,
          showProgress = F,
          na.strings = c("", "NA"),
          col.names = blNames)
      }else{
        blf <- NULL
      }
      if(!is.na(x$targetBlast)){
        blr <- fread(
          x$targetBlast,
          showProgress = F,
          na.strings = c("", "NA"),
          col.names = blNamesR)
      }else{
        blr <- NULL
      }
      bl <- rbind(blf, blr, use.names = T)
      bl <- bl[,blNames, with = F]
      if(nrow(bl) == 0)
        stop("could not read in the blast files\n")

      # -- 2.2 drop duplicated hits, keeping the highest score
      setorder(bl, -bitScore)
      bl <- subset(bl, !duplicated(paste(ofID1, ofID2)))

      # -- 2.3 merge with bed information
      bl <- merge(bed1, merge(bed2, bl, by = "ofID2"), by = "ofID1")

      # -- 2.4 get anchor and og information
      bl[,sameOG := og1 == og2 & !is.na(og1) & !is.na(og2)]
      bl[,noAnchor := noAnchor1 | noAnchor2]
      bl[,`:=`(og1 = NULL, og2 = NULL, noAnchor1 = NULL, noAnchor2 = NULL)]

      bl[,`:=`(isAnchor = NA, inBuffer = NA, isSyntenic = NA,
               blkID = NA, regID = NA)]

      write_allBlast(bl, filepath = x$allBlast)
      x[,`:=`(nTotalHits = nrow(bl), nGlobOgHits = sum(bl$sameOG))]
      return(x)
    }))

    with(out, cat(sprintf("\t...%stotal hits = %s, same og = %s\n",
                          lab, nTotalHits, nGlobOgHits)))
    return(out)
  }))
  if(any(out$nGlobOgHits == 0)){
    stop(
      "some blast files have 0 hits in the same orthogroup.
      This is not to be expected and indicates that there is a
      problem with your orthofinder run.")
  }

  gsParam$synteny$blast <- outmd
  return(gsParam)
}

