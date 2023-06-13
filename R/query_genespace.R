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
#' @param refGenome character matching one of the genomeIDs if a bed object is
#' not provided
#' @param showUnPlacedPgs logical, when queried should un-located orthogroups be
#' shown?
#' @param onlyArrayReps logical, should only array representatives be considered
#' for CNV counts?
#' @param maxCopyNumber numeric, maximum copyNumber to tabulate. Orthogroups
#' with more than this number are combined into the largest CNV category

#' @return a list of data.tables named $genome, $chr: $start-$end. One
#' data.table for each line in the bed file.

#' @title query function for hits files
#' @description
#' \code{query_hits} extract syntenic (or all) hits that fall within a interval
#' in a genome.
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

#' @title reformat and subset pan-gene sets
#' @description
#' \code{query_pangenes} the primary engine to explore genespace output, this
#' lets you reformat the pan-genes to wide (entry - by - genome) matrix and
#' only look at specified positional bounds.
#' @rdname query_genespace
#' @import data.table
#' @export
query_pangenes <- function(gsParam,
                           bed = NULL,
                           refGenome = NULL,
                           transform = TRUE,
                           showArrayMem = TRUE,
                           showNSOrtho = TRUE,
                           maxMem2Show = Inf,
                           showUnPlacedPgs = FALSE){

  ##############################################################################
  parse_pangenes <- function(gsParam,
                             refGenome,
                             transform,
                             showArrayMem,
                             showNSOrtho,
                             maxMem2Show){
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
      pgo[,genome := factor(genome, levels = gsParam$genomeIDs)]
      pgo <- dcast(pgo, pgID + interpChr + interpOrd + og + repGene + pgRepID ~ genome,
                   value.var = "id", fun.aggregate = list)
      pgo <- merge(pgo, cbm, by = "pgRepID")
      pgo[,pgRepID := NULL]
      setcolorder(
        pgo,
        c("pgID", "interpChr", "interpOrd", "og", "repGene", "genome","chr",
          "start", "end", gsParam$genomeIDs))
    }

    setkey(pgo, pgID)
    return(pgo)
  }
  ##############################################################################

  x <- flag <- id <- repGene <- pgRepID <- ofID <- og <- pgID <- index <-
    pass <- genome <- chr <- start <- end <- regID <- ord <- interpChr <-
    interpOrd <- flag <- id <- repGene <- pgRepID <- ofID <- og <- pgID <-
    index <- NULL

  cb <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  cbm <- cb[,c("ofID", "genome", "chr", "start", "end")]
  setnames(cbm, "ofID", "pgRepID")
  if(is.null(bed)){
    if(is.null(refGenome))
      stop("must specify either a bed or refGenome\n")
    out <- parse_pangenes(
      gsParam = gsParam,
      refGenome = refGenome,
      transform = transform,
      showArrayMem = showArrayMem,
      showNSOrtho = showNSOrtho,
      maxMem2Show = maxMem2Show)
    if(!showUnPlacedPgs)
      out <- subset(out, !is.na(interpOrd))
    return(out)
  }else{
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

    out <- lapply(1:nrow(bed), function(i){
      x <- bed[i,]
      out <- parse_pangenes(
        gsParam = gsParam,
        refGenome = x$genome,
        transform = transform,
        showArrayMem = showArrayMem,
        showNSOrtho = showNSOrtho,
        maxMem2Show = maxMem2Show)
      tmp <- subset(
        out,
        genome == x$genome & chr == x$chr & end <= x$end & start >= x$start)
      outroi <- subset(out, pgID >= min(tmp$pgID) & pgID <= max(tmp$pgID))

      if(!showUnPlacedPgs)
        outroi <- subset(outroi, !is.na(interpOrd))

      return(outroi)
    })
    names(out) <- bed$regID

    return(out)
  }
}


#' @title reformat and subset pan-gene sets
#' @description
#' \code{query_cnv} the primary engine to explore genespace output, this
#' lets you reformat the pan-genes to wide (entry - by - genome) matrix and
#' only look at specified positional bounds.
#' @rdname query_genespace
#' @import data.table
#' @export
query_cnv <- function(gsParam,
                      bed = NULL,
                      onlyArrayReps = FALSE,
                      maxCopyNumber = Inf){

  count_cnv <- function(combBed, ofIDs){

    ofID <- ogType <- genome <- og <- copyNumber <- nGenes <- nOgs <- NULL
    cbl <- melt(
      combBed,
      measure.vars = c("globOG", "globHOG", "og"),
      id.vars = c("genome", "ofID", "chr", "start", "end"),
      variable.name = "ogType",
      value.name = "og")

    if(!is.null(ofIDs)){
      cbogs <- subset(cbl, ofID %in% ofIDs)[,c("ogType", "og")]
      cbogs <- subset(cbogs, !duplicated(cbogs))
      cbl <- merge(cbl, cbogs, by = c("ogType", "og"))
    }

    cbl[,ogType := ifelse(ogType == "og", "syntenicHOGs",
                          ifelse(ogType == "globHOG", "globalHOGs", "globalOGs"))]

    cbn <- cbl[,list(copyNumber = .N), by = c("genome", "og", "ogType")]
    cbn$copyNumber[cbn$copyNumber > maxCopyNumber] <- maxCopyNumber
    eg <- cbl[,CJ(genome = unique(genome), og = unique(og)), by = "ogType"]
    cbn <- merge(cbn, eg, by = colnames(eg), all = T)
    cbn$copyNumber[is.na(cbn$copyNumber)] <- 0
    cbo <- cbn[,list(nGenes = copyNumber*uniqueN(og), nOgs = uniqueN(og)),
               by = c("genome", "ogType", "copyNumber")]
    setkey(cbo, genome, ogType, copyNumber)
    cbo[,`:=`(percGenes = 100*(nGenes/sum(nGenes)),
              percOGs = 100*(nOgs/sum(nOgs))),
        by = c("genome", "ogType")]
    return(cbo)
  }

  isArrayRep <- pass <- genome <- chr <- start <- end <- regID <- genome <- NULL

  cb <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  if(onlyArrayReps){
    cb <- subset(cb, isArrayRep)
  }

  if(!is.null(bed)){
    bed <- data.table(bed)
    if(nrow(bed) == 0)
      stop("bed must be a data.table or data.frame with >= 1 rows\n")
    if(!all(c("chr", "genome") %in% names(bed)))
      stop("bed must be a data.table or data.frame with columns named chr and genome\n")
    u <- unique(paste(cb$genome, cb$chr))
    bed[,pass := paste(genome, chr) %in% u]

    if(!all(bed$pass))
      stop("the following genome/chr combinations are not in the run:",
           subset(bed, !pass)[,1:2])
    bed[,pass := NULL]
    if(!"start" %in% names(bed))
      bed[,start := 0]
    if(!"end" %in% names(bed))
      bed[,end := Inf]
    bed[,regID := sprintf("%s, %s: %s-%s", genome, chr, start, end)]


    out <- lapply(1:nrow(bed), function(i){
      x <- bed[i,]
      bogs <- subset(
        cb, genome == x$genome & chr == x$chr & end >= x$start & start <= x$end)
      outi <- count_cnv(combBed = cb, ofIDs = bogs$ofID)
      return(outi)
    })
    names(out) <- bed$regID
  }else{
    out <- count_cnv(combBed = cb, ofIDs = NULL)
  }
 return(out)
}
