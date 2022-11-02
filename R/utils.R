#' @title Generic internal functions used by genespace
#' @description
#' \code{utils} Convience functions for genespace, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param x single-value parameter, string, integer, numeric, list, vector
#' @param path file path character string
#' @param min if x is an integer or numeric, the minimum value allowed
#' @param max if x is an integer or numeric, the minimum value allowed
#' @param default if there is a problem with x, replace with this value
#' @param na.rm logical, should NA's be dropped
#' @param onlySingleValue logical, should long a single valuebe returned?
#' @param genomeIDs character vector of genomeIDs
#' @param to integer, top end value
#' @param which character specifying what to use
#' @param n integer, number of observations
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'

#' @title Check parameter value for integer values or vectors
#' @description
#' \code{check_integer} Checks and parses integer arguments to GENESPACE
#' functions. Replaces values of x not in range (min, max) with the minimum or
#' maximum values. If a single value and a value that is not coercible to an
#' integer is specified, returns the default value. If na.rm = TRUE and
#' onlySingleValue = FALSE, drops NAs.
#' from the vector
#' @rdname utils
#' @export
check_integer <- function(x,
                          min = -Inf,
                          max = Inf,
                          default = NA,
                          na.rm = FALSE,
                          onlySingleValue = length(x) <= 1){

  if(is.null(x))
    x <- NA

  suppressWarnings(min <- as.integer(min(x, na.rm = T)))
  suppressWarnings(max <- as.integer(max(x, na.rm = T)))
  suppressWarnings(na.rm <- as.logical(na.rm[1]))
  suppressWarnings(onlySingleValue <- as.logical(onlySingleValue[1]))

  if(is.na(min))
    min <- (-Inf)
  if(is.na(max))
    max <- (-Inf)
  if(is.na(onlySingleValue))
    onlySingleValue <- length(x) <= 1

  suppressWarnings(x <- as.integer(x))
  if(na.rm)
    x <- x[!is.na(x)]

  if(onlySingleValue || length(x) <= 1){
    x <- x[1]
    if(is.na(x) || is.null(x) || length(x) == 0){
      suppressWarnings(default <- as.integer(default[1]))
      x <- default
    }else{
      if(x > max) x <- max
      if(x < min) x <- min
    }
  }else{
    x[x > max & !is.na(x)] <- max
    x[x < min & !is.na(x)] <- min
  }

  return(x)
}

#' @title Check parameter value for numeric values or vectors
#' @description
#' \code{check_numeric} See check_integer. Same but for numeric values.
#' @rdname utils
#' @export
check_numeric <- function(x,
                          min = -Inf,
                          max = Inf,
                          default = NA,
                          na.rm = FALSE,
                          onlySingleValue = length(x) <= 1){
  if(is.null(x))
    x <- NA

  suppressWarnings(min <- as.numeric(min(x, na.rm = T)))
  suppressWarnings(max <- as.numeric(max(x, na.rm = T)))
  suppressWarnings(na.rm <- as.logical(na.rm[1]))
  suppressWarnings(onlySingleValue <- as.logical(onlySingleValue[1]))

  if(is.na(min))
    min <- (-Inf)
  if(is.na(max))
    max <- (-Inf)
  if(is.na(onlySingleValue))
    onlySingleValue <- length(x) <= 1

  suppressWarnings(x <- as.numeric(x))
  if(na.rm)
    x <- x[!is.na(x)]

  if(onlySingleValue || length(x) <= 1){
    x <- x[1]
    if(is.na(x) || is.null(x) || length(x) == 0){
      suppressWarnings(default <- as.numeric(default[1]))
      x <- default
    }else{
      if(x > max) x <- max
      if(x < min) x <- min
    }
  }else{
    x[x > max & !is.na(x)] <- max
    x[x < min & !is.na(x)] <- min
  }

  return(x)
}


#' @title Check parameter value for character values or vectors
#' @description
#' \code{check_character} See check_integer. Same but for character values.
#' @rdname utils
#' @export
check_character <- function(x,
                            default = NULL,
                            na.rm = FALSE,
                            onlySingleValue = length(x) <= 1){
  if(is.null(x))
    x <- NA
  if(all(is.na(x)))
    x <- rep("NA", length(x))
  suppressWarnings(na.rm <- as.logical(na.rm[1]))
  suppressWarnings(onlySingleValue <- as.logical(onlySingleValue[1]))

  if(is.na(onlySingleValue))
    onlySingleValue <- length(x) <= 1

  suppressWarnings(x <- as.character(x))
  if(na.rm)
    x <- x[!is.na(x)]

  if(onlySingleValue || length(x) <= 1){
    x <- x[1]
    if(is.na(x) || is.null(x) || length(x) == 0){
      suppressWarnings(default <- as.character(default[1]))
      if(is.na(default))
        stop(sprintf("NA value given, but default value %s could not be coerced to a character\n", default))
      x <- default
    }
  }
  x[x == "NA"] <- NA
  return(x)
}

#' @title Check parameter value for logical values or vectors
#' @description
#' \code{check_logical} See check_integer. Same but for logical values.
#' @rdname utils
#' @export
check_logical <- function(x,
                          default = NA,
                          na.rm = FALSE,
                          onlySingleValue = length(x) <= 1){

  if(is.null(x))
    x <- NA
  suppressWarnings(na.rm <- as.logical(na.rm[1]))
  suppressWarnings(onlySingleValue <- as.logical(onlySingleValue[1]))

  if(is.na(onlySingleValue))
    onlySingleValue <- length(x) <= 1

  suppressWarnings(x <- as.logical(x))
  if(na.rm)
    x <- x[!is.na(x)]
  if(length(x) == 0)
    x <- NA

  if(onlySingleValue || length(x) <= 1){
    x <- x[1]
    if(is.na(x) || is.null(x) || length(x) == 0){
      suppressWarnings(default <- as.logical(default[1]))
      x <- default
    }
  }

  return(x)
}

#' @title check_filePathParam
#' @description
#' \code{check_filePathParam} check_filePathParam
#' @rdname utils
#' @export
check_filePathParam <- function(path){
  if(is.null(path))
    path <- NA
  path <- path.expand(check_character(
    x = path, default = NA, na.rm = FALSE))
  chk <- dir.exists(path) || file.exists(path)
  if(!chk)
    path <- NA
  return(path)
}

#' @title Check only DNA
#' @description
#' \code{check_onlyDNA} Check only DNA
#' @rdname utils
#' @importFrom Biostrings readAAStringSet DNA_ALPHABET AAStringSet
#' @export
check_onlyDNA <- function(path){
  fa <- readAAStringSet(path)
  tmp <- gsub("[^A-Za-z]", "", paste(fa[1:100], collapse = ""))
  dnaa <- paste(DNA_ALPHABET, collapse = "|")
  tmp <- gsub(gsub("[^A-Za-z]", "", dnaa), "", tmp)
  return(nchar(tmp) == 0)
}


#' @title read_aaFasta
#' @description
#' \code{read_aaFasta} read_aaFasta
#' @rdname utils
#' @export
read_aaFasta <- function(path){
  chk <- tryCatch(
    {
      readAAStringSet(path)
    },
    error = function(err) {
      return(NA)
    }
  )
}

#' @title count the number of amino acids by gene
#' @description
#' \code{get_nAA} count the number of amino acids by gene
#' @rdname utils
#' @importFrom Biostrings readAAStringSet
#' @export
get_nAA <- function(path){
  y <- readAAStringSet(path)
  o <- width(y)
  names(o) <- names(y)
  return(o)
}

#' @title read_bed
#' @description
#' \code{read_bed} read_bed
#' @rdname utils
#' @import data.table
#' @importFrom stats complete.cases
#' @export
read_bed <- function(path){
  chk <- tryCatch(
    {
      suppressWarnings(suppressMessages(fread(
        path, verbose = FALSE, showProgress = FALSE, select = 1:4,
        colClasses = c("character", "numeric", "numeric", "character"),
        header = FALSE, col.names = c("chr", "start", "end", "id"))))
    },
    error = function(err) {
      return(NA)
    }
  )
  if(!is.data.table(chk))
    chk <- subset(chk, complete.cases(chk))

  return(chk)
}

#' @title read_bed
#' @description
#' \code{read_bed} read_bed
#' @rdname utils
#' @import data.table
#' @export
align_charLeft <- function(x){
  x <- check_character(x, default = "", na.rm = FALSE, onlySingleValue = FALSE)
  maxChar <- max(nchar(x))
  nbuff <- maxChar - nchar(x)
  buffBlank <- sapply(nbuff, function(z) paste(rep(" ", z), collapse = ""))
  return( x <- sprintf("%s%s", x, buffBlank))
}

#' @title read_bed
#' @description
#' \code{read_bed} read_bed
#' @rdname utils
#' @import data.table
#' @export
align_charRight <- function(x){
  x <- check_character(x, default = "", na.rm = FALSE, onlySingleValue = FALSE)
  maxChar <- max(nchar(x))
  nbuff <- maxChar - nchar(x)
  buffBlank <- sapply(nbuff, function(z) paste(rep(" ", z), collapse = ""))
  return(sprintf("%s%s", buffBlank, x))
}

#' @title Read orthofinder species IDs
#' @description
#' \code{read_orthofinderSpeciesIDs} Parses the SpeciesIDs.txt file into a
#' data.table and returns to R.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSpeciesIDs <- function(path){
  si <- fread(
    path,
    sep = ":",
    header = F,
    col.names = c("genomeNum", "genome"),
    colClasses = c("numeric", "character"))

  genome <- NULL
  si[,genome := gsub(".fa", "", genome, fixed = T)]
  sio <- si$genomeNum; names(sio) <- si$genome
  return(sio)
}

#' @title Read orthofinder sequence IDs
#' @description
#' \code{read_orthofinderSequenceIDs} Reads the sequence
#' IDs:gene name dictionary into memory.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSequenceIDs <- function(path){

  gi <- fread(
    path,
    header = F,
    sep = ":",
    col.names = c("ofID","id"),
    colClasses = c("character","character"))

  ofID <- NULL
  gi[,c("genomeNum", "geneNum") := tstrsplit(ofID, "_", type.convert = T)]
  return(gi)
}

#' @title Count the number of sequences in a fasta file
#' @description
#' \code{get_nSeqs} Counts the number of lines with ">" in a file. If the output
#' is not convertible to an integer, returns NA.
#' @rdname utils
#' @export
get_nSeqs <- function(path){
  nseq <- system2("grep", "-c '>'", path, stdout = TRUE, stderr = TRUE)
  suppressWarnings(nseq <- as.integer(nseq))
  return(nseq)
}


#' @title check_annotFiles
#' @description
#' \code{check_annotFiles}
#' @rdname utils
#' @import data.table
#' @export
check_annotFiles <- function(path, genomeIDs){

  # -- make sure the working directory is OK
  wd <- check_filePathParam(path)
  if(is.na(wd))
    stop("couldn't find working directory:", wd)

  # -- make sure the peptide files exist
  pepDir <- file.path(wd, "peptide")
  pepFiles <- file.path(pepDir, sprintf("%s.fa",genomeIDs))
  if(!all(file.exists(pepFiles))){
    mf <- basename(pepFiles[!file.exists(pepFiles)])
    stop(sprintf("Cannot find the following peptide files in %s:\n\t%s\n",
                 pepDir, paste(mf, collapse = "\n\t")))
  }

  # -- make sure the gff files exist
  bedDir <- file.path(wd, "bed")
  bedFiles <- file.path(bedDir, sprintf("%s.bed",genomeIDs))
  if(!all(file.exists(bedFiles))){
    mf <- basename(bedFiles[!file.exists(bedFiles)])
    stop(sprintf("Cannot find the following bed files in %s:\n\t%s\n",
                 bedDir, paste(mf, collapse = "\n\t")))
  }

  # -- check that the peptide files contain peptides and not DNA
  onlyDNA <- sapply(pepFiles, check_onlyDNA)
  if(any(onlyDNA)){
    mf <- basename(names(onlyDNA)[onlyDNA])
    warning(sprintf("The following fasta files appear to only have DNA sequence\n\t%s\nCheck these to make sure they contain peptides and not DNA\n",
                    paste(mf, collapse = "\n\t")))
  }

  # -- check exact match between bed and fasta headers
  labs <- align_charLeft(genomeIDs)
  chk <- rbindlist(lapply(1:length(genomeIDs), function(i){

    aa <- read_aaFasta(pepFiles[i])
    bd <- read_bed(bedFiles[i])
    un <- uniqueN(c(names(aa), bd$id))
    ui <- length(intersect(names(aa), bd$id))
    pf <- ifelse((ui/un) < .95, "FAIL", "PASS")
    cat(sprintf("\t%s: %s / %s geneIDs exactly match (%s)\n",
                labs[i], ui, un, pf))
    return(data.table(
      genomeID = genomeIDs[i],
      bedFile = bedFiles[i],
      peptideFile = pepFiles[i],
      nGenes = ui,
      annotationPass = pf == "PASS"))
  }))
  if(any(!chk$annotationPass))
    stop("some annotations do not have matching peptide headers and bed names\n")
  return(chk)
}


#' @title check MCScanX install
#' @description
#' \code{check_MCScanXhInstall} check that MCScanX_h can be called
#' @rdname utils
#' @export
check_MCScanXhInstall <- function(path){
  path <- check_filePathParam(path)
  if(is.na(path)){
    chk <- NA
  }else{
    pth <- check_filePathParam(file.path(path, "MCScanX_h"))
    if(is.na(pth)){
      chk <- NA
    }else{
      chk <- suppressWarnings(grepl(
        "prefix_fn",
        system2(pth, "-h", stdout = TRUE, stderr = FALSE)[1]))
      if(!chk)
        chk <- NA
    }
  }
  return(chk)
}

#' @title parse_ogs
#' @description
#' \code{parse_ogs} parse_ogs
#' @rdname utils
#' @import data.table
#' @export
parse_ogs <- function(path, genomeIDs){
  id <- genome <- Orthogroup <- NULL
  tmp <- fread(path, showProgress = F, verbose = F)
  tmp <- melt(
    tmp, id.vars = "Orthogroup", measure.vars = genomeIDs,
    variable.name = "genome", value.name = "id")
  tmp <- tmp[,list(id = strsplit(id, ",")[[1]]),
             by = c("Orthogroup", "genome")]
  tmp[,`:=`(genome = trimws(genome),
            id = trimws(id), Orthogroup = trimws(Orthogroup))]
  setnames(tmp, 1, "ogID")
  return(tmp)
}

#' @title parse_hogs
#' @description
#' \code{parse_hogs} parse_hogs
#' @rdname utils
#' @import data.table
#' @export
parse_hogs <- function(path, genomeIDs){
  id <- genome <- HOG <- NULL
  tmp <- fread(path, showProgress = F, verbose = F)
  tmp <- melt(
    tmp, id.vars = "HOG", measure.vars = genomeIDs,
    variable.name = "genome", value.name = "id")
  tmp <- tmp[,list(id = strsplit(id, ",")[[1]]), by = c("HOG", "genome")]
  tmp[,`:=`(genome = trimws(genome), id = trimws(id), HOG = trimws(HOG))]
  setnames(tmp, 1, "hogID")
  return(tmp)
}

#' @title parse_orthologues
#' @description
#' \code{parse_orthologues} parse_orthologues
#' @rdname utils
#' @import data.table
#' @export
parse_orthologues <- function(path){
  orthID <- id1 <- id2 <- NULL
  x <- fread(path, showProgress = F)
  refID <- colnames(x)[2]
  altID <- colnames(x)[3]
  setnames(x, c("og", "id1", "id2"))
  x1 <- subset(x, !grepl(",", paste(id1, id2)))
  x2 <- subset(x, grepl(",", paste(id1, id2)))
  x2[,orthID := 1:.N]
  x2r <- x2[,list(id1 = unique(strsplit(id1, ",")[[1]])),
            by = "orthID"]
  x2a <- x2[,list(id2 = unique(strsplit(id2, ",")[[1]])),
            by = "orthID"]
  x2 <- merge(x2r, x2a, by = "orthID", all = T, allow.cartesian = T)
  x1[,orthID := (1:.N)+max(x2$orthID)]
  x <- rbind(x1[,colnames(x2), with = F], x2)
  x[,`:=`(gen1 = refID, gen2 = altID,
          id1 = gsub(" ", "", id1), id2 = gsub(" ", "", id2))]
  return(x)
}

#' @title round_toInteger
#' @description
#' \code{round_toInteger} round_toInteger
#' @rdname utils
#' @import data.table
#' @export
round_toInteger <- function(x, to){
  round(x/to, 0) * to
}

#' @title add_rle
#' @description
#' \code{add_rle} add_rle
#' @rdname utils
#' @import data.table
#' @export
add_rle <- function(x, which = "n"){
  if(which == "n"){
    rep(rle(x)$lengths, rle(x)$lengths)
  }else{
    rep(1:length(rle(x)$lengths), rle(x)$lengths)
  }
}

#' @title gs_colors
#' @description
#' \code{gs_colors} gs_colors
#' @rdname utils
#' @import data.table
#' @importFrom grDevices colorRampPalette
#' @export
gs_colors <- function(n = 10){
  cols <- c("#C4645C", "#F5915C", "#FFC765",
            "#FCF8D5", "#BEF3F9", "#66B8FF", "#6666FF", "#9C63E1",
            "#F4BDFF")
  pal <- colorRampPalette(cols)
  return(pal(n))
}

#' @title clus_igraph
#' @description
#' \code{clus_igraph} Clus_igraph
#' @rdname utils
#' @importFrom igraph graph_from_data_frame clusters
#' @export
clus_igraph <- function(id1, id2){
  if(length(unique(id1)) == 1 & length(unique(id2)) == 2){
    return(rep(1, length(id1)))
  }else{
    clus <- clusters(graph_from_data_frame(
      data.frame(id1, id2),
      directed = F))$membership[id1]
    clus <- clus[!duplicated(names(clus))]
    return(clus)
  }
}

#' @title linear interpolation of missing values
#' @description
#' \code{interp_linear} linear interpolation of missing values in x
#' @rdname utils
#' @import data.table
#' @export
interp_linear <- function(x,
                          y){

  # -- convert to numeric, to ensure that NAs are correctly specified
  ord1 <- as.numeric(x)
  ord2 <- as.numeric(y)
  if(all(is.na(ord1)))
    stop("refOrd (achors) is all NAs. Can't proceed\n")
  if(any(is.na(ord2)))
    stop("found NAs in toInterpOrd (hits to interpolate) - not permitted\n")

  # -- if no NAs, just spit back ord2
  if(all(!is.na(ord1))){
    return(ord1)
  }else{
    # make into a data table with index and whether or not to interpolate
    dt <- data.table(
      x = ord1, y = ord2, index = 1:length(ord1), toInterp = is.na(ord1))
    dto <- data.table(dt)
    # order by anchor positions (ord1, x)
    setkey(dt, y)

    # -- find runs of NAs in y
    dt[,rl := add_rle(toInterp, which = "id")]

    # -- pull runs to infer (not first and last if they are NAs)
    interpThis <- subset(dt, !(toInterp & rl %in% c(1, max(rl))))

    # -- return original data if no bounding non-na runs and no internal NAs
    if(uniqueN(interpThis$rl[interpThis$toInterp]) < 1){
      return(ord1)
    }else{
      # -- get max right and min left values for each non-missing run
      minr <- with(subset(interpThis, !toInterp), tapply(x, rl, min))
      maxl <- with(subset(interpThis, !toInterp), tapply(x, rl, max))

      # -- linear interpolation of runs of NAs from bounding values
      out <- subset(interpThis, toInterp)
      out[,ip := seq(from = maxl[as.character(rl-1)],
                     to = minr[as.character(rl+1)],
                     length.out = .N+2)[-c(1, .N+2)],
          by = "rl"]

      # -- fill NAs and return
      dto$x <- as.numeric(dto$x)
      dto$x[out$index] <- out$ip
      return(dto$x)
    }
  }
}

#' @title flag_boundingNAs
#' @description
#' \code{flag_boundingNAs} flag_boundingNAs
#' @rdname utils
#' @export
flag_boundingNAs <- function(x){
  z <- rep(FALSE, length(x))
  y <- add_rle(is.na(x), which = "id")
  if(is.na(x)[1])
    z[y == 1] <- TRUE

  if(is.na(x)[length(x)])
    z[y == max(y)] <- TRUE

  return(z)
}

#' @title calc_blkCoords
#' @description
#' \code{calc_blkCoords} calc_blkCoords
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
  ofID1 <- start1 <- end1 <- ofID1 <- ord1 <- NULL
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


#' @title check if a vector is coercible to R colors
#' @description
#' \code{are_colors} check if a vector is coercible to R colors
#' @rdname utils
#' @importFrom grDevices col2rgb
#' @export
are_colors <- function(col) {
  sapply(col, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}


#' @title scale a vector between a range
#' @description
#' \code{scale_between} scale a vector between a range
#' @rdname utils
#' @export
scale_between <- function(x, min, max, scale1toMean = TRUE){
  if(length(unique(x)) > 1){
    return((x - min(x)) / (max(x) - min(x)) * (max - min) + min)
  }else{
    if(scale1toMean){
      return(mean(c(min, max)))
    }else{
      return(max)
    }
  }
}

#' @title read_combBed
#' @description
#' \code{read_combBed} read_combBed
#' @rdname utils
#' @export
read_combBed <- function(filepath){
  bedNames <- c(
    "chr", "start", "end", "id", "ofID", "pepLen", "ord", "genome", "arrayID",
    "isArrayRep", "globOG", "globHOG", "synOG", "inblkOG", "noAnchor", "og")
  cl <- c("numeric", "character", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, bedNames))
    stop("combBed.txt file is malformed\n")
  bedClass <- cl[c(2,1,1,2,2,1,1,2,2,3,2,2,2,2,3,2)]
  bed <- fread(
    filepath, na.strings = c("", "NA"), select = bedNames,
    colClasses = bedClass, showProgress = F)
  return(bed)
}

#' @title write_combBed
#' @description
#' \code{write_combBed} write_combBed
#' @rdname utils
#' @export
write_combBed <- function(x, filepath){
  bedNames <- c(
    "chr", "start", "end", "id", "ofID", "pepLen", "ord", "genome", "arrayID",
    "isArrayRep", "globOG", "globHOG", "synOG", "inblkOG", "noAnchor", "og")
  if(!all(bedNames %in% colnames(x)))
    stop("bed data.table is malformed\n")
  x <- x[,bedNames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}


#' @title read_synHits
#' @description
#' \code{read_synHits} read_synHits
#' @rdname utils
#' @export
read_synHits <- function(filepath){
  hnames <- c(
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "pid", "length", "mismatches", "gapopenings", "queryStart", "queryEnd",
    "subjectStart", "subjectEnd", "Evalue", "bitScore", "sameOg", "noAnchor",
    "isAnchor", "inBuffer", "regID", "blkID", "lgBlkID")
  cl <- c("numeric", "character", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, hnames))
    stop("synhits.txt file is malformed\n")
  hc <- cl[c(2,2,1,1,2,1,2,3,2,2,1,1,2,1,2,3,1,1,1,1,1,1,1,1,1,1,3,3,3,3,2,2,2)]
  hits <- fread(
    filepath, na.strings = c("", "NA"), select = hnames,
    colClasses = hc, showProgress = F)
  return(hits)
}

#' @title write_synBlast
#' @description
#' \code{write_synBlast} write_synBlast
#' @rdname utils
#' @export
write_synBlast <- function(x, filepath){
  hnames <- c(
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "pid", "length", "mismatches", "gapopenings", "queryStart", "queryEnd",
    "subjectStart", "subjectEnd", "Evalue", "bitScore", "sameOg", "noAnchor",
    "isAnchor", "inBuffer", "regID", "blkID", "lgBlkID")
  if(!all(hnames %in% colnames(x)))
    stop("synhits data.table is malformed\n")
  x <- x[,hnames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}


#' @title add transparency to a color
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col,
                      alpha = 1){

  if(missing(col) || !all(are_colors(col)))
    stop("Colors are misspecified\n")
  if(length(alpha) != 1 || alpha > 1 || alpha < 0)
    stop("alpha is misspecified\n")

  return(apply(sapply(col, col2rgb)/255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha)))
}

#' @title read_intSynPos
#' @description
#' \code{read_intSynPos} read_intSynPos
#' @rdname utils
#' @export
read_intSynPos <- function(filepath){
  hnames <- c("genome", "ofID", "chr", "ord", "og", "interpGenome",
              "interpChr", "interpOrd", "isAnchor")
  hclass <- c("character", "character", "character", "numeric", "character",
              "character", "character", "numeric", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, hnames))
    stop("integratedSynPos file is malformed\n")
  synpos <- fread(
    filepath, na.strings = c("", "NA"), select = hnames,
    colClasses = hclass, showProgress = F)
  return(synpos)
}
