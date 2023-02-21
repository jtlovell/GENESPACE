#' @title Query GENESPACE results
#' @name query_genespace
#' @description
#' \code{query_genespace} Functions to pull data stored in the /pangenome
#' and /syntenicHits directories
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param bed A data.table similar to the format of a bed file. At least two
#' columns must be specified: genome, chr. These must match genome-chromosome
#' combinations in the data. The user can also supply two additional columns:
#' start, end. These are the base-pair start and end coordinates of each region.
#' If start and end are not specified, they are set to 0 and Inf respectively
#' so that the entire chromosome is returned.
#' @param synOnly logical, given to query hits to return (or not) the synteny-
#' constrained hits only
#' @param transform logical, should the pangenome be transformed from the long
#' to the wide format?
#' @param showArrayMem logical, should all genes or only the array
#' representative genes (if FALSE) be shown?
#' @param showNSOrtho logical, should non-syntenic orthologs be included in the
#' output?
#' @param maxMem2Show integer specifying the maximum number of members to be
#' included in each entry-by-flag-by-genome combination.

#' @return a list of data.tables named $genome, $chr: $start-$end. One
#' data.table for each line in the bed file.

#' @title query_hits
#' @description
#' \code{query_hits} query_hits
#' @rdname query_genespace
#' @import data.table
#' @export
query_hits <- function(gsParam,
                       bed,
                       synOnly = TRUE){

  pass <- genome <- chr <- start <- end <- chr1 <- NULL

  cb <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  bed <- data.table(bed)
  if(nrow(bed) == 0)
    stop("bed must be a data.table or data.frame with >= 1 rows\n")
  if(!all(c("chr", "genome") %in% names(bed)))
    stop("bed must be a data.table or data.frame with columns named chr and genome\n")
  u <- unique(paste(cb$genome, cb$chr))
  bed[,pass := paste(genome, chr) %in% u]
  if(!all(bed$pass))
    stop("the following genome/chr combinations are not in the run:", subset(bed, !pass)[,1:2])
  bed[,pass := NULL]
  if(!"start" %in% names(bed))
    bed[,start := 0]
  if(!"end" %in% names(bed))
    bed[,end := Inf]

  spl <- split(bed, by = "genome")
  combout <- rbindlist(lapply(names(spl), function(i){
    inbed <- with(data.table(spl[[i]]), data.table(
      chr1 = chr, start1 = start, end1 = end,
      id = sprintf("%s, %s: %s-%s", genome, chr, start, end),
      key = c("chr1", "start1", "end1")))
    if(synOnly){
      rh <- read_refGenomeSynHits(gsParam = gsParam, refGenome = i)
    }else{
      rh <- read_refGenomeAllBlast(gsParam = gsParam, refGenome = i)
    }
    rh <- subset(rh, chr1 %in% inbed$chr)
    setkeyv(rh, c("chr1", "start1", "end1"))

    out <- foverlaps(inbed, rh)
    return(out)
  }))
  return(split(combout, by = "id"))
}

#' @title query_pangenome
#' @description
#' \code{query_pangenome} query_pangenome
#' @rdname query_genespace
#' @import data.table
#' @export
query_pangenes <- function(gsParam,
                            bed = NULL,
                            refGenome = NULL,
                            transform = TRUE,
                            showArrayMem = TRUE,
                            showNSOrtho = TRUE,
                            maxMem2Show = Inf){

  x <- flag <- id <- repGene <- pgRepID <- ofID <- og <- pgID <- index <-
    pass <- genome <- chr <- start <- end <- regID <- ord <- interpChr <-
    interpOrd <- flag <- id <- repGene <- pgRepID <- ofID <- og <- pgID <-
    index <- NULL


  if(is.null(bed)){
    if(is.null(refGenome))
      stop("must specify either a bed or refGenome\n")
    pgo <- fread(file.path(gsParam$paths$pangenes, sprintf(
      "%s_pangenes.txt.gz", refGenome)))

    if(!showArrayMem)
      pgo <- subset(pgo, flag != "array")

    if(!showNSOrtho)
      pgo <- subset(pgo, flag != "NSOrtho")

    if(transform){
      pgo[,flag := ifelse(flag == "NSOrtho", "*",
                          ifelse(flag == "array", "+", ""))]
      pgo[,id := sprintf("%s%s", id, flag)]
      pgo[,repGene := id[pgRepID == ofID][1], by = "pgID"]
      pgo[,og := og[flag == ""][1], by = "pgID"]
      if(is.finite(maxMem2Show)){
        setkey(pgo, pgID, flag)
        pgo[,index := 1:.N, by = c("pgID", "genome","flag")]
        pgo <- subset(pgo, index <= maxMem2Show)
      }
      pgo <- dcast(pgo, pgID + interpChr + interpOrd + og + repGene ~ genome,
                   value.var = "id", fun.aggregate = list)
    }
    return(pgo)
  }else{
    cb <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
    bed <- data.table(bed)
    if(nrow(bed) == 0)
      stop("bed must be a data.table or data.frame with >= 1 rows\n")
    if(!all(c("chr", "genome") %in% names(bed)))
      stop("bed must be a data.table or data.frame with columns named chr and genome\n")
    u <- unique(paste(cb$genome, cb$chr))
    bed[,pass := paste(genome, chr) %in% u]
    if(!all(bed$pass))
      stop("the following genome/chr combinations are not in the run:", subset(bed, !pass)[,1:2])
    bed[,pass := NULL]
    if(!"start" %in% names(bed))
      bed[,start := 0]
    if(!"end" %in% names(bed))
      bed[,end := Inf]
    bed[,regID := sprintf("%s, %s: %s-%s", genome, chr, start, end)]

    # -- convert bed to gene order positions
    setkey(cb, genome, chr, start, end)
    setkey(bed, genome, chr, start, end)
    conv <- foverlaps(cb, bed)
    bedo <- subset(conv, !is.na(regID))[,list(start = min(ord), end = max(ord)),
                                        by = c("genome", "chr", "regID")]

    out <- lapply(1:nrow(bedo), function(i){
      x <- bedo[i,]
      pg <- fread(file.path(gsParam$paths$pangenes, sprintf(
        "%s_pangenes.txt.gz", x$genome)))
      pgo <- subset(
        pg, interpChr == x$chr & interpOrd >= x$start & interpOrd <= x$end)

      if(!showArrayMem)
        pgo <- subset(pgo, flag != "array")

      if(!showNSOrtho)
        pgo <- subset(pgo, flag != "NSOrtho")

      if(transform){
        pgo[,flag := ifelse(flag == "NSOrtho", "*",
                            ifelse(flag == "array", "+", ""))]
        pgo[,id := sprintf("%s%s", id, flag)]
        pgo[,repGene := id[pgRepID == ofID][1], by = "pgID"]
        pgo[,og := og[flag == ""][1], by = "pgID"]
        if(is.finite(maxMem2Show)){
          setkey(pgo, pgID, flag)
          pgo[,index := 1:.N, by = c("pgID", "genome","flag")]
          pgo <- subset(pgo, index <= maxMem2Show)
        }
        pgo <- dcast(pgo, pgID + interpChr + interpOrd + og + repGene ~ genome,
                     value.var = "id", fun.aggregate = list)
      }
      return(pgo)
    })
    names(out) <- bedo$regID
    return(out)
  }
}
