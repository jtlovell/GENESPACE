#' @title convert_input2v1
#'
#' @description
#' \code{convert_input2v1} xxx
#'
#' @param existingDir A
#' @param v1Dir A
#' @details T...
#'
#' @return nothing
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#
#' @import data.table
#' @export
query_genespace <- function(gsParam,
                            regions){

  ##############################################################################
  # 1. Ensure that the input data looks ok
  if(!all(colnames(regions) %in% c("chr", "start", "end", "id", "genome", "dataType", "regID")))
    stop("regions input must be a six-column data.frame (or data.table) with the
         columns: chr, start, end, id, genome, and dataType. The genome
         and dataType columns MUST be specified, but the others are optional.
         If start or end are given, chr must also be given. If id is given,
         start, end and chr are ignored\n")

  if(!"regID" %in% colnames(regions) && nrow(regions) > 1)
    warning("regID is not a column in regions. Regions will be labeled by row number (reg1, reg2, ... regN)")
  if(any(duplicated(regions$regID))){
    warning("found duplicated region IDs. re-naming by row number (reg1, reg2, ... regN)")
    regions[,regID := NULL]
  }

  if(!"regID" %in% colnames(regions))
    regions[,regID := sprintf("region_%s", 1:.N)]

  if(nrow(regions) < 1)
    stop("regions input is empty\n")

  if(any(is.na(regions$genome)))
    stop("All elements of the regions genome column MUST be specified\n")

  if(any(is.na(regions$dataType)))
    stop("All elements of the regions dataType column MUST be specified\n")

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))

  if(!all(regions$genome %in% unique(bed$genome)))
    stop("Some genome entries in the regions input do not match the genome IDs\n")

  if(!all(regions$dataType %in% c("hits", "pangenome", "bed", "all")))
    stop("Some dataType entries in the regions input are not hits, pangenome, bed or all\n")

  if(any(is.na(regions$chr) & !is.na(regions$start)) ||
     any(is.na(regions$chr) & !is.na(regions$end)))
    stop("If start or end coordinates are specified in regions, must also give the chromosome\n")

  ##############################################################################
  # 2. Query the bed file
  regBed <- subset(regions, dataType %in% c("bed", "all"))
  if(nrow(regBed) > 0){
    bedList <- lapply(1:nrow(regBed), function(i){
      x <- regBed[i,]
      if(!is.na(x$id)){
        tmp <- subset(bed, genome == x$genome & id %in% x$id)
      }else{
        if(is.na(x$chr)){
          tmp <- subset(bed, genome == x$genome)
        }else{
          st <- ifelse(is.na(x$start), -Inf, x$start)
          en <- ifelse(is.na(x$end), Inf, x$end)
          tmp <- subset(
            bed, genome == x$genome & end >= st & start <= en & chr == x$chr)
        }
      }
      return(subset(bed, og %in% tmp$og))
    })
    names(bedList) <- regBed$regID
  }else{
    bedList <- NULL
  }

  ##############################################################################
  # 3. Query the hits files
  regHits <- subset(regions, dataType %in% c("hits", "all"))
  if(nrow(regHits) > 0){
    md <- data.table(gsParam$annotBlastMd)
    hitList <- lapply(1:nrow(regHits), function(i){
      x <- regHits[i,]
      hs <- subset(md, target == x$genome | query == x$genome & !is.na(annotBlastFile))
      if(nrow(hs) > 0){
        if(is.na(x$id)){
          st <- ifelse(is.na(x$start), -Inf, x$start)
          en <- ifelse(is.na(x$end), Inf, x$end)
          tmp <- subset(
            bed, genome == x$genome & end >= st & start <= en & chr == x$chr)
          ids <- tmp$id
        }else{
          ids <- x$id
        }
        hits <- rbindlist(lapply(hs$annotBlastFile, function(y)
          subset(read_synHits(y),
                 (id1 %in% ids & genome1 == x$genome) |
                   (id2 %in% ids & genome2 == x$genome))))
        return(hits)
      }
    })
    names(hitList) <- regHits$regID
  }else{
    hitList <- NULL
  }

  ##############################################################################
  # 4. Query the pangenome
  regPg <- subset(regions, dataType %in% c("pangenome", "all"))
  if(nrow(regHits) > 0){
    pgList <- lapply(1:nrow(regPg), function(i){
      x <- regPg[i,]
      pg <- fread(file.path(gsParam$paths$pangenome, sprintf(
        "%s_refPangenomeAnnot.txt", x$genome)))
      if(!is.na(x$id)){
        pgids <- unique(subset(pg, id %in% x$id)$pgID)
        pgw <- subset(pg, pgID %in% pgids)
      }else{
        st <- ifelse(is.na(x$start), -Inf, x$start)
        en <- ifelse(is.na(x$end), Inf, x$end)
        tmp <- subset(
          bed, genome == x$genome & end >= st & start <= en & chr == x$chr)
        rng <- range(tmp$ord)
        pgw <- subset(pg, pgOrd >= rng[1] & pgOrd <= rng[2] & pgChr == x$chr)
      }

      pgw[,flag := ifelse(isNSOrtho, "*", ifelse(!isArrayRep, "+", ""))]
      pgw[,id := sprintf("%s%s", id, flag)]
      pgw[,repGene := id[isRep][1], by = "pgID"]
      pgw <- dcast(pgw, pgID + pgGenome + pgChr + pgOrd + og + repGene ~ genome,
                   value.var = "id", fun.aggregate = list)
      return(pgw)
    })
    names(pgList) <- regPg$regID
  }else{
    pgList <- NULL
  }

  return(list(bedRegions = bedList,
              hitRegions = hitList,
              pangenomeRegions = pgList))
}
