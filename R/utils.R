#' @title Generic internal functions used by genespace
#' @description
#' \code{utils} Convience functions for genespace, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param gsParam list of genespace parameters
#' @param path file path character string
#' @param col character or integer specifying a color
#' @param alpha numeric (0-1) specifying the transparency
#' @param pattern regular expression to search for
#' @param recursive logical, should the search be recursive
#' @param which character specifying which method to use
#' @param x vector of observations
#' @param y vector of observations
#' @param id1 character vector of ids with length matching id2
#' @param id2 character vector of ids with length matching id1
#' @param min numeric, length 1 specifying the minumum value in the scaled data
#' @param max numeric, length 1 specifying the maximum value in the scaled data
#' @param raw logical, should a raw vector of peptide widths be returned?
#' @param scale1toMean logical, if single value, should it be the mean?
#' @param width integer, the number of characters to return in string.
#' @param to what should the value be rounded to?
#' @param refGenome character string specifying which genome is the reference
#' @param blFile file path to the blast-like text file
#' @param ofID1 orthofinder ID of the first gene
#' @param ofID2 orthofinder ID of the second gene
#' @param onlyIDScore logical, should only the geneIDs and score be returned?
#' @param onlyCheckRun logical, should the run just be checked?
#' @param hits data.table containing annotated blast-format pairwise hits
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'
#' @title Order files by date of last modification
#' @description
#' \code{order_filesByMtime} Order files by date of last modification
#' @rdname utils
#' @export
order_filesByMtime <- function(path = getwd(),
                               pattern = "*",
                               recursive = F){
  if (length(path) == 1) {
    allFiles <- list.files(
      path = path,
      full.names = T,
      pattern = pattern,
      recursive = recursive)
  }else{
    allFiles <- path
  }

  details <- file.info(allFiles, extra_cols = F)
  details <- details[rev(with(details, order(as.POSIXct(mtime)))), ]
  return(rownames(details))
}

#' @title Check if orthofinder is installed
#' @description
#' \code{check_orthofinderInstall} Check if orthofinder is installed
#' @rdname utils
#' @export
check_orthofinderInstall <- function(path){
  grepl("OrthoFinder",
        system2(path, "-h", stdout = TRUE)[2])
}

#' @title Check logical argument
#' @description
#' \code{check_logicalArg} Ensure a logical (length 1) argument is coerced
#' correctly. If cannot be coerced, return error.
#' @rdname utils
#' @importFrom parallel detectCores
#' @export
check_logicalArg <- function(x){
  x <- as.logical(x)
  if(is.null(x)) stop(x, "must be logical\n")
  if(is.na(x)) stop(x, "must be logical\n")
  if(length(x) == 0) stop(x, "must be logical\n")
  if(length(x) > 1) x <- x[1]
  return(x)
}

#' @title convert vector to RLE
#' @description
#' \code{add_rle} Convert vector into run length equivalents
#' @rdname utils
#' @export
add_rle <- function(x, which = "n"){
  if (which == "n") {
    rep(rle(x)$lengths, rle(x)$lengths)
  }else{
    rep(1:length(rle(x)$lengths), rle(x)$lengths)
  }
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
    return(clusters(graph_from_data_frame(
      data.frame(id1, id2),
      directed = F))$membership[id1])
  }
}

#' @title count the number of amino acids by gene
#' @description
#' \code{get_nAA} count the number of amino acids by gene
#' @rdname utils
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
get_nAA <- function(path, raw = FALSE){
  setDTthreads(1)
  if(!raw){
    pepF <- list.files(path, pattern = "^Species", full.names = T)
    pepF <- pepF[grep(".fa$", pepF)]
    peps <- rbindlist(lapply(pepF, function(x){
      y <- readAAStringSet(x)
      return(data.table(ofID = names(y),
                        nAA = width(y)))
    }))
  }else{
    y <- readAAStringSet(path)
    o <- width(y)
    names(o) <- names(y)
    return(o)
  }
}

#' @title calculate the mode
#' @description
#' \code{find_modalValue} find the most commmon value in a series
#' @rdname utils
#' @export
find_modalValue <- function(x){
  tab <- table(x)
  return(names(tab[order(-tab)])[1])
}

#' @title scale a vector between a range
#' @description
#' \code{scale_between} scale a vector between a range
#' @rdname utils
#' @export
scale_between <- function(x, min, max){
  (x - min(x)) / (max(x) - min(x)) * (max - min) + min
}

#' @title flatten a list
#' @description
#' \code{flatten_list} convert a list into a vector while keeping names
#' @rdname utils
#' @export
flatten_list <- function(x){
  y <- unlist(x)
  names(y) <- rep(names(x), sapply(x, length))
  return(y)
}

#' @title
#' @description
#' \code{scale_between} ...
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

#' @title pull_strWidth
#' @description
#' \code{pull_strWidth} pull_strWidth
#' @rdname utils
#' @export
pull_strWidth <- function(x, width){
  y <- substr(x, 1, width)
  if(nchar(y) < width)
    y <- sprintf("%s%s", y, paste(rep(" ", width-nchar(y)), collapse = ""))
  return(y)
}

#' @title round_toInteger
#' @description
#' \code{round_toInteger} round_toInteger
#' @rdname utils
#' @export
round_toInteger <- function(x, to){
  round(x / to, 0) * to
}

#' @title Read orthofinder species IDs
#' @description
#' \code{read_orthofinderSpeciesIDs} Parses the SpeciesIDs.txt file into a
#' data.table and returns to R.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSpeciesIDs <- function(path){
  setDTthreads(1)
  genome <- NULL
  si <- fread(
    file.path(path, "SpeciesIDs.txt"),
    sep = ":",
    header = F,
    col.names = c("genomeNum", "genome"),
    colClasses = c("numeric", "character"))
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
  setDTthreads(1)
  ofID <- NULL
  gi <- fread(
    file.path(path, "SequenceIDs.txt"),
    header = F,
    sep = ":",
    col.names = c("ofID","id"),
    colClasses = c("character","character"))
  gi[,c("genomeNum","geneNum") := tstrsplit(ofID, "_", type.convert = T)]
  return(gi)
}

#' @title choose_mostRecentOF
#' @description
#' \code{choose_mostRecentOF} choose_mostRecentOF
#' @rdname utils
#' @import data.table
#' @export
choose_mostRecentOF <- function(path){
  setDTthreads(1)
  p <- file.path(path, "Orthofinder")
  ps <- order_filesByMtime(p)
  pschk <- lapply(ps, function(x){
    chks <- c(dir.exists(file.path(x, "Orthologues")),
              dir.exists(file.path(x, "Gene_Duplication_Events")),
              file.exists(file.path(x, "Gene_Duplication_Events","Duplications.tsv")))
    if(all(chks))
      return(x)
  })
  wh <- min(which(!sapply(pschk, is.null)))
  if(length(wh) > 0){
    return(pschk[[wh]])
  }else{
    stop("Cannot find orthogues and duplications, has the full orthofinder pipe been run?\n")
  }
}

#' @title parse orthogroups file into a data.table
#' @description
#' \code{parse_ogs} wide to long format conversion for orthogroups.tsv
#' @rdname utils
#' @export
parse_ogs <- function(gsParam){
  setDTthreads(1)
  id <- genome <- Orthogroup <- NULL
  ogtsv <- file.path(gsParam$paths$orthogroupsDir, "Orthogroups.tsv")
  og <- fread(ogtsv, showProgress = F, verbose = F)
  og <- melt(
    og, id.vars = "Orthogroup", variable.name = "genome", value.name = "id")
  og <- og[,list(id = strsplit(id, ",")[[1]]), by = c("Orthogroup", "genome")]
  og[,`:=`(genome = trimws(genome), id = trimws(id), Orthogroup = trimws(Orthogroup))]
  setnames(og, 1, "ogID")
  return(og)
}

#' @title parse_orthologues
#' @description
#' \code{parse_orthologues} parse_orthologues
#' @rdname utils
#' @import data.table
#' @export
parse_orthologues <- function(gsParam, refGenome){
  orthID <- NULL
  setDTthreads(1)
  od <- file.path(gsParam$paths$orthologuesDir,
                  sprintf("Orthologues_%s", refGenome))
  odf <- list.files(od, full.names = T, pattern = "__v__")
  ogo <- rbindlist(lapply(odf, function(i){
    x <- fread(i, showProgress = F)
    refID <- colnames(x)[2]
    altID <- colnames(x)[3]
    setnames(x, c("og", "id1", "id2"))
    id1 <- id2 <- NULL
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
  }))
  return(ogo)
}

#' @title Read orthofinder blast file
#' @description
#' \code{read_blast} Reads in a single pairwise orthofinder-formaatted blast
#' file
#' @rdname utils
#' @import data.table
#' @export
read_blast <- function(blFile = NULL,
                       ofID1 = NULL,
                       ofID2 = NULL,
                       path = NULL,
                       onlyIDScore = TRUE){
  setDTthreads(1)
  if(is.null(blFile)){
    blFile <- file.path(path, sprintf("Blast%s_%s.txt.gz", ofID1, ofID2))
  }
  if(!file.exists(blFile))
    stop("cannot find ", blFile, "\n")

  if(!onlyIDScore){
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE)
    g1 <- strsplit(bl$V1[1], "_")[[1]][1]
    g2 <- strsplit(bl$V2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[, c(2, 1, 3:6, 8, 7, 10, 9, 11, 12)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl), with = F]
      bl <- rbind(bl, tmp)
      V12 <- NULL
      setorder(bl, -V12)
      bl <- subset(bl, !duplicated(bl[, c(1:2)]))
    }
  }else{
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE,
      select = c(1,2,12),
      col.names = c("ofID1", "ofID2", "score"))
    g1 <- strsplit(bl$ofID1[1], "_")[[1]][1]
    g2 <- strsplit(bl$ofID2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[,c(2,1,3)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[, colnames(bl), with = F]
      bl <- rbind(bl, tmp)
      score <- NULL
      setorder(bl, -score)
      bl <- subset(bl, !duplicated(bl[, c(1:2)]))
    }
  }

  return(bl)
}

#' @title find_orthofinderResults
#' @description
#' \code{find_orthofinderResults} find_orthofinderResults
#' @rdname utils
#' @import data.table
#' @importFrom R.utils gunzip
#' @export
find_orthofinderResults <- function(gsParam, onlyCheckRun = F){
  ogsFile <- order_filesByMtime(
    path = gsParam$paths$orthofinder,
    pattern = "Orthogroups.tsv",
    recursive = T)
  if(onlyCheckRun){
    return(length(ogsFile) > 0)
  }else{
    if(length(ogsFile) > 1)
      warning("Found multiple orthofinder runs, only using the most recent\n")
    if(length(ogsFile) == 0)
      stop("Can't find the 'orthogroups.tsv' file\n\tHave you run orthofinder yet?\n")

    ofResDir <- dirname(dirname(ogsFile[1]))

    pfile <- file.path(ofResDir, "Gene_Duplication_Events")
    if(!dir.exists(pfile)){
      paralogsDir <- NA
    }else{
      paralogsDir <- pfile
    }

    orthfile <- file.path(ofResDir, "Orthologues")
    if(!dir.exists(orthfile)){
      orthologuesDir <- NA
    }else{
      orthologuesDir <- orthfile
    }

    blsFile <- order_filesByMtime(
      path = gsParam$paths$orthofinder,
      pattern = "diamondDBSpecies0.dmnd",
      recursive = T)

    blastDir <- dirname(blsFile[1])

    gsParam$paths$blastDir <- blastDir
    gsParam$paths$orthogroupsDir <- dirname(ogsFile[1])
    gsParam$paths$paralogsDir <- paralogsDir
    gsParam$paths$orthologuesDir <- orthologuesDir

    return(gsParam)
  }
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
             "chr1", "chr2", "gen1", "gen2", "ofID1", "ofID2")
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
    by = c("blkID", "gen1","gen2", "chr1", "chr2")]

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
    by = c("blkID", "gen1","gen2", "chr1", "chr2")]

  # -- merge the two coordinates
  blks <- merge(blks1, blks2, by = c("gen1","gen2","chr1","chr2","blkID"))

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

#' @title clus_dbscan
#' @description
#' \code{clus_dbscan} clus_dbscan
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
clus_dbscan <- function(hits,
                        radius,
                        blkSize){
  setDTthreads(1)

  # -- split up the hits
  x <- split(hits, by = c("gen1","gen2","chr1","chr2"))

  ord1 <- ord2 <- chr1 <- chr2 <- NULL
  x <- rbindlist(lapply(x,  function(y){
    y[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = radius),
                       minPts = blkSize)$cluster]
    return(y)
  }))

  # -- subset the blocks, then add unique IDs
  blkID <- NULL
  x <- subset(x, blkID > 0)
  x[,blkID := sprintf("blk_%s_%s_%s", chr1, chr2, blkID)]

  return(x)
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

#' @title drop_unusedPeptides
#' @description
#' \code{drop_unusedPeptides} drop_unusedPeptides
#' @rdname utils
#' @import data.table
#' @export
drop_unusedPeptides <- function(gsParam){
  f <- list.files(path = dirname(gsParam$paths$peptide[1]), full.names = F)
  fi <- basename(gsParam$paths$peptide)
  if(any(!f %in% fi)){
    fo <- f[!f %in% fi]
    for(i in fo)
      file.remove(file.path(dirname(gsParam$paths$peptide[1]), i))
  }
}

#' @title linear interpolation of missing values
#' @description
#' \code{interp_linear} linear interpolation of missing values in x
#' @rdname utils
#' @import data.table
#' @export
interp_linear <- function(refOrd,
                          toInterpOrd){

  # -- convert to numeric, to ensure that NAs are correctly specified
  dt <- x <- y <- NULL
  ord1 <- as.numeric(refOrd)
  ord2 <- as.numeric(toInterpOrd)
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
