#' @title Generic internal functions used by genespace
#' @description
#' \code{utils} Convience functions for genespace, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param gsParam list of genespace parameters
#' @param path file path character string
#' @param pattern regular expression to search for
#' @param recursive logical, should the search be recursive
#' @param which character specifying which method to use
#' @param x vector of observations
#' @param id1 character vector of ids with length matching id2
#' @param id2 character vector of ids with length matching id1
#' @param min numeric, length 1 specifying the minumum value in the scaled data
#' @param max numeric, length 1 specifying the maximum value in the scaled data
#' @param dt data.table with the first two columns as id1, id2.
#' @param raw logical, should a raw vector of peptide widths be returned?
#' @param scale1toMean logical, if single value, should it be the mean?
#' @param width integer, the number of characters to return in string.
#' @param x x position of the scale bar
#' @param to what should the value be rounded to?
#' @param byGrpCol column which serves as the by factor
#' @param windowSize integer specifying the window size
#' @param yCol character specifying the column name containing the y values
#' @param fun function to pass to sliding window
#' @param gffFiles vector of file paths pointing to the gff files
#' @param refGenome character string specifying which genome is the reference
#' @param nCores integer, the number of parallel processes to run
#' @param blFile file path to the blast-like text file
#' @param ofID1 orthofinder ID of the first gene
#' @param ofID2 orthofinder ID of the second gene
#' @param onlyIDScore logical, should only the geneIDs and score be returned?
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param blastDir file.path to the location of the blast results.
#' @param onlyCheckRun logical, should the run just be checked?
#' @param genomeIDs an optional vector of genomeIDs to consider. If not
#' specified (default) taken from gsParam$genomeIDs$genomeIDs
#' @param hits data.table containing annotated blast-format pairwise hits
#' @param radius numeric of length 1 specifying the eps dbscan parameter; the
#' search radius within which to count clustered density-based xy points.
#' @param blkSize integer of length 1 specifying the minimum size for a syntenic
#' block and the -s 'size' MCScanX parameter
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#' @param minRbhScore integer of length 1, see set_syntenyParams
#' @param genome1 character string specifying first of two genomeIDs
#' @param genome2 character string specifying second of two genomeIDs
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param synBuff integer of length 1 specifying the maximum euclidean distance
#' from an 'anchor' so that it can be considered syntenic
#' @param selfOnly logical, should only self hits be considered
#' @param overwrite logical, should the results be overwrittem?
#' @param maskHits data.table of hits that should be excluded
#' @param synParam data.table with synteny parameters. See set_syntenyParams.
#' @param selfRegionMask integer, the radius around self hits that should be
#' masked
#' @param type either 'primary' or 'secondary' depending on the scale of
#' inference
#' @param dropInterleavesSmallerThan integer, the minimum block size to retain
#' after splitting overlapping blocks
#' @param minPropDup numeric (0-1) specifying the minimum proportion of
#' duplicated hits to allow two overlapping blocks to not be split
#' @param maxIter integer, the maximum number of block splitting interations
#' @param nhits integer, the number of hits to retain
#' @param blks data.table containing the block coordinates
#' @param useBlks logical, for cross compatibility with plot_hits
#' @param allowRBHinOg logical, should reciprocal best blast hits be used when
#' finalizing block coordinates?
#' @param verbose logical, should updates be printed to the console?
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
        system(paste(path, "-h"),
               intern = T)[2])
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
#' @import data.table
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

#' @title complete subgraphs
#' @description
#' \code{complete_subgraph} given a hits data.table, expand unique orthofinder
#' IDs among all unique elements in a subgraph
#' @rdname utils
#' @export
complete_subgraph <- function(dt){
  clusters <- clusters(graph_from_data_frame(
    dt,
    directed = F))
  clusters <- with(clusters, data.table(
    id = names(membership),
    group = membership))
  out <- data.table(merge(
    clusters,
    clusters,
    by = "group",
    allow.cartesian = T))
  setnames(out, c("cluster", "x", "y"))
  return(out)
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


#' @title checks parameter for minimum length of a peptide
#' @description
#' \code{check_minPepLen} checks for init_genespace
#' @rdname utils
#' @export
check_minPepLen <- function(x){
  x <- as.integer(x)[1]
  if(is.na(x) || is.null(x))
    stop("minPepLen must be an integer\n")
  if(x < 0)
    x <- 0
  return(x)
}

#' @title checks parameter for interleave size between blocks
#' @description
#' \code{check_dropInterSize} checks for init_genespace
#' @rdname utils
#' @export
check_dropInterSize <- function(x){
  x <- as.integer(x)[1]
  if(is.na(x) || is.null(x))
    stop("dropInterleavesSmallerThan must be an integer\n")
  if(x < 1)
    x <- 1
  return(x)
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

#' @title calc_partialSw
#' @description
#' \code{calc_partialSw} calc_partialSw
#' @rdname utils
#' @export
calc_partialSw <- function(dt,
                           byGrpCol,
                           windowSize,
                           yCol,
                           fun = function(x) mean(x, na.rm = T)){
  yDataFakeCol <- byGrpFakeCol <- swLeft <- swRight <- NULL
  dt <- data.table(dt)
  # rename variables
  if("yDataFakeCol" %in% colnames(dt) & yCol != "yDataFakeCol")
    stop("yDataFakeCol in colnames but not ycol, rename please\n")
  if("byGrpFakeCol" %in% colnames(dt) & byGrpCol != "byGrpFakeCol")
    stop("byGrpFakeCol in colnames but not byGrpCol, rename please\n")

  if(yCol != "yDataFakeCol")
    dt[,yDataFakeCol := dt[[yCol]]]
  if(byGrpCol != "byGrpFakeCol")
    dt[,byGrpFakeCol := dt[[byGrpCol]]]

  # apply sliding windows
  dt[,`:=`(swLeft = frollapply(yDataFakeCol, n = windowSize, align = "left", FUN = fun),
           swRight = frollapply(yDataFakeCol, n = windowSize, align = "right", FUN = fun)),
     by = "byGrpFakeCol"]

  # replace NAs with the alternative means
  dt$swLeft[is.na(dt$swLeft)] <- dt$swRight[is.na(dt$swLeft)]
  dt$swRight[is.na(dt$swRight)] <- dt$swLeft[is.na(dt$swRight)]

  # get mean of left and right aligned windows
  dt[,`:=`(sw = (swLeft + swRight)/2,
           swLeft = NULL, swRight = NULL)]

  # drop renamed variables
  if(yCol != "yDataFakeCol")
    dt[,yDataFakeCol := NULL]
  if(byGrpCol != "byGrpFakeCol")
    dt[,byGrpFakeCol := NULL]
  return(dt)
}

#' @title load gff annotations
#' @description
#' \code{read_gff} reads a genespace-formatted gff-like annotation file into
#' memory
#' @rdname utils
#' @export
read_gff <- function(gffFiles){
  genome <- NULL
  gff <- rbindlist(lapply(names(gffFiles), function(i){
    x <- fread(
      gffFiles[[i]],
      key = c("chr", "start","end","strand"))
    x[,genome := i]
    return(x)
  }))
  gffCols <- c("ord","chr","genome","start", "end", "strand","id")
  if(!all(gffCols %in% colnames(gff)))
    stop(paste(gffCols, collapse = ", "),
         " must all be column names in gff\n")
  return(gff)
}

#' @title Read orthofinder species IDs
#' @description
#' \code{read_orthofinderSpeciesIDs} Parses the SpeciesIDs.txt file into a
#' data.table and returns to R.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSpeciesIDs <- function(path){
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
#' @importFrom parallel mclapply
#' @export
parse_orthologues <- function(gsParam, refGenome, nCores = 1){
  id1 <- id2 <- orthID <- NULL
  od <- file.path(gsParam$paths$orthologuesDir,
                  sprintf("Orthologues_%s", refGenome))
  odf <- list.files(od, full.names = T, pattern = "__v__")
  ogo <- rbindlist(mclapply(odf, mc.cores = nCores, function(i){
    x <- fread(i, showProgress = F)
    refID <- colnames(x)[2]; altID <- colnames(x)[3]
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
  V12 <- score <- NULL
  if(is.null(blFile)){
    blFile <- file.path(path, paste0("Blast", ofID1, "_", ofID2,".txt.gz"))
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
      tmp <- data.table(bl[,c(2,1,3:6,8,7,10,9,11,12)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl),with = F]
      bl <- rbind(bl, tmp)
      setorder(bl, -V12)
      bl <- subset(bl, !duplicated(bl[,c(1:2)]))
    }
  }else{
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE,
      select = c(1,2,12),
      col.names = c("ofID1","ofID2","score"))
    g1 <- strsplit(bl$ofID1[1], "_")[[1]][1]
    g2 <- strsplit(bl$ofID2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[,c(2,1,3)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl),with = F]
      bl <- rbind(bl, tmp)
      setorder(bl, -score)
      bl <- subset(bl, !duplicated(bl[,c(1:2)]))
    }
  }

  return(bl)
}

#' @title add orthofinder ID to a gff object
#' @description
#' \code{add_ofID2gff} read the orthofinder species and gene IDs and merge
#' these with the gff-like data.table
#' @rdname utils
#' @export
add_ofID2gff <- function(gff,
                         blastDir){
  id <- ofID <- genomeNum <- genome <- NULL
  specIDs <- read_orthofinderSpeciesIDs(blastDir)
  gv <- names(specIDs); names(gv) <- as.character(specIDs)
  seqIDs <- read_orthofinderSequenceIDs(blastDir)
  seqIDs[,genome :=  gv[as.character(genomeNum)]]
  idv <- seqIDs$ofID; names(idv) <- with(seqIDs, paste(genome, id))
  gff[,ofID := idv[paste(genome, id)]]
  return(gff)
}

#' @title add peptide length to a gff object
#' @description
#' \code{add_pepLen2gff} read the peptide lengths and merge with the gff
#' @rdname utils
#' @export
add_pepLen2gff <- function(gff,
                           gsParam){
  pepLen <- id <- ofID <-  NULL
  spl <- split(gff, by = "genome")
  naa <- rbindlist(lapply(spl, function(x){
    x[,pepLen := get_nAA(gsParam$paths$peptide[x$genome[1]], raw = T)[id]]
    return(x[,c("ofID","pepLen")])
  }))
  nao <- naa$pepLen; names(nao) <- naa$ofID
  gff[,pepLen := nao[ofID]]
  return(gff)
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

#' @title prep_ofDbFromPeptide
#' @description
#' \code{prep_ofDbFromPeptide} prep_ofDbFromPeptide
#' @rdname utils
#' @import data.table
#' @importFrom R.utils gunzip
#' @export
prep_ofDbFromPeptide <- function(gsParam){

  # clean up leftover peptides if necessary
  drop_unusedPeptides(gsParam)

  # check for and remove orthofinder directory if needed
  if(dir.exists(gsParam$paths$orthofinder)){
    unlink(gsParam$paths$orthofinder, recursive = T)
  }

  # convert to orthofinder
  com <- sprintf(
    "%s -f %s -t %s -op -o %s 1>/dev/null 2>&1",
    gsParam$paths$orthofinderCall,
    dirname(gsParam$paths$peptide[1]),
    gsParam$params$nCores,
    gsParam$paths$orthofinder)
  system(com)
  return(com)
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
calc_blkCoords <- function(hits){
  ord1 <- ord2 <- start1 <- start2 <- end1 <- end2 <- ofID1 <- ofID2 <-minBp2 <- NULL
  blkID <- orient <- maxBp2 <- maxBp1 <- maxOrd2 <- minOrd2 <-maxGene2 <- minGene2 <- NULL

  bhits <- subset(hits, !is.na(blkID))
  setkey(bhits, ord1)
  blks1 <- bhits[,list(startBp1 = min(start1, na.rm = T), endBp1 = max(end1, na.rm = T),
                       startOrd1 = min(ord1, na.rm = T), endOrd1 = max(ord1, na.rm = T),
                       firstGene1 = first(ofID1), lastGene1 = last(ofID1),
                       nHits1 = uniqueN(ofID1)),
                 by = c("blkID", "gen1","gen2", "chr1", "chr2")]
  setkey(bhits, ord2)
  blks2 <- bhits[,list(minBp2 = min(start2, na.rm = T), maxBp2 = max(end2, na.rm = T),
                       minOrd2 = min(ord2, na.rm = T), maxOrd2 = max(ord2, na.rm = T),
                       minGene2 = first(ofID2), maxGene2 = last(ofID2, 1),
                       nHits2 = uniqueN(ofID2),
                       orient = ifelse(length(ord1) <= 1, "+",
                                       ifelse(cor(jitter(ord1),
                                                  jitter(ord2)) > 0,"+", "-"))),
                 by = c("blkID")]
  blks <- merge(blks1, blks2, by = "blkID")

  bgfor <- subset(blks, orient == "+")
  bgrev <- subset(blks, orient == "-")
  bgrev[,`:=`(startBp2 = maxBp2, endBp2 = minBp2,
              startOrd2 = maxOrd2, endOrd2 = minOrd2,
              firstGene2 = maxGene2, lastGene2 = minGene2)]
  bgfor[,`:=`(startBp2 = minBp2, endBp2 = maxBp2,
              startOrd2 = minOrd2, endOrd2 = maxOrd2,
              firstGene2 = minGene2, lastGene2 = maxGene2)]
  blks <- rbind(bgfor, bgrev)
  return(blks)
}


#' @title finalize_blocks
#' @description
#' \code{finalize_blocks} finalize_blocks
#' @rdname utils
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
finalize_blocks <- function(hits,
                            synParam,
                            gsParam,
                            type = "primary",
                            verbose = FALSE,
                            minRbhScore = 50){

  isOg <- maxScr1 <- maxScr2 <- score <- og <- ofID1 <- ofID2 <- ord1 <- ord2 <- NULL
  isCollin <- blkID <- bord1 <- bord2 <- gen1 <- gen2 <- nblk <- regID <- NULL
  regAnchor <- blkIDn <- chr1 <- chr2 <- isAnchor <- NULL
  nhits1 <- nSecondHits1 <- nhits2 <- nSecondHits2 <- blkSizeSecond <- NULL
  synBuffSecond <- nGapsSecond <- NULL
  ##############################################################################
  # 1. Pull anchor hits and re-rank
  ##############################################################################
  sp <- data.table(synParam[1,])
  if(type == "Second")
    sp[,`:=`(nhits1 = (nhits1-1 + nSecondHits1),
             nhits2 = (nhits2-1 + nSecondHits2),
             blkSize = blkSizeSecond,
             nGaps = nGapsSecond,
             synBuff = synBuffSecond)]

  dropInterleavesSmallerThan <- sp$dropInterleavesSmallerThan
  blkSize <- sp$blkSize
  nGaps <- sp$nGaps
  synBuff <- sp$synBuff
  path2mcscanx <- gsParam$paths$mcscanxCall
  nCores <- gsParam$params$nCores

  anchs <- subset(hits, regAnchor & !is.na(regID))

  # -- calculate rank order within regions
  anchs[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")),
        by = "regID"]

  ##############################################################################
  # 2. remove non-colinear hits within regions
  ##############################################################################
  # -- choose representative hits within each region
  spl <- split(anchs, by = "regID")
  colAnch <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){
    suppressWarnings(x[,isCollin := !is.na(run_mcscanx(
      hits = x,
      gsParam = gsParam,
      blkSize = blkSize,
      nGaps = nGaps,
      path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)])])
    if(!"isCollin" %in% colnames(x))
      x[,isCollin := FALSE]
    return(x)
  }), use.names=TRUE)

  # -- drop non colinear hits and re-rank
  colAnch <- subset(colAnch, isCollin)
  colAnch[,isCollin := NULL]
  colAnch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                ord2 = frank(ord2, ties.method = "dense")),
          by = "regID"]

  ##############################################################################
  # 3. cluster within regions into blks
  ##############################################################################
  # -- initial clustering, allowing to drop some non-clustering points
  spl <- split(colAnch, by = "regID")
  iblks <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){
    x[,blkID := dbscan(frNN(cbind(ord1, ord2),
                            eps = sqrt((blkSize^2)*2)),
                       minPts = blkSize)$cluster]
    return(x)
  }))
  iblks <- subset(iblks, blkID != 0)

  # -- clustering within blks, again dropping blocks too small
  iblks[,blkID := sprintf("%s_%s_%s_%s",gen1, gen2, regID, blkID)]
  iblks[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")),
        by = "blkID"]
  iblks[,nblk := dbscan(frNN(cbind(ord1, ord2),
                             eps = sqrt(2)+.1),
                        minPts = blkSize)$cluster,
        by = "blkID"]
  fblks <- subset(iblks, blkID != 0)

  # -- final clustering within blks
  fblks[,blkID := sprintf("%s_%s", blkID, nblk)]
  fblks[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")),
        by = "blkID"]
  fblks[,nblk := dbscan(frNN(cbind(ord1, ord2),
                             eps = sqrt(2)+.1),
                        minPts = 0)$cluster,
        by = "blkID"]

  # -- split interleaved blocks
  outBlks <- fblks[,colnames(hits), with = F]

  blv <- outBlks$blkID; names(blv) <- with(outBlks, paste(ofID1, ofID2))
  blAnc <- with(outBlks, paste(ofID1, ofID2))

  ##############################################################################
  # 4. format and return blk regs
  ##############################################################################

  # get block regions and hits therein
  hits[,blkID := blv[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(hits, !is.na(blkID)))
  inReg <- get_hitsInBlks(blks = bc, hits = hits, nCores = nCores)
  blv <- inReg$blkID; names(blv) <- with(inReg, paste(ofID1, ofID2))

  # -- get hits in buffer from block anchors
  inReg[,isAnchor := paste(ofID1, ofID2) %in% blAnc]
  inBuff <- find_hitsInBuff(
    hits = inReg, nCores = nCores, synBuff = synBuff)
  buffv <- with(inBuff, paste(ofID1, ofID2))

  # -- return a new hits obj
  hits[,`:=`(blkID = blv[paste(ofID1, ofID2)],
             blkAnchor = paste(ofID1, ofID2) %in% blAnc,
             blkBuffer = paste(ofID1, ofID2) %in% buffv)]
  return(hits)
}

#' @title parse_blast4synteny
#' @description
#' \code{parse_blast4synteny} parse_blast4synteny
#' @rdname utils
#' @import data.table
#' @export
parse_blast4synteny <- function(gsParam,
                                genome1,
                                genome2,
                                gff,
                                selfOnly){

  score <- ofID1 <- ofID2 <- isOg <- og1 <- scrRank1 <- scrRank2 <- NULL
  chr1 <- chr2 <- ord1 <- ord2 <- NULL
  # -- get vectors from the gff
  gv <- gff$genome; cv <- gff$chr; ogv <- gff$globOG
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(gv) <- names(cv) <- names(ogv) <- gff$ofID
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  # -- read in the blast hits
  ofsp <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  if(genome1 == genome2){
    bl <- read_blast(
      path = gsParam$paths$blastDir,
      ofID1 = ofsp[genome1],
      ofID2 = ofsp[genome1])
  }else{
    bl <- rbind(read_blast(
      path = gsParam$paths$blastDir,
      ofID1 = ofsp[genome1],
      ofID2 = ofsp[genome2]),
      with(read_blast(
        path = gsParam$paths$blastDir,
        ofID1 = ofsp[genome2],
        ofID2 = ofsp[genome1]),
        data.table(
          ofID1 = ofID2, ofID2 = ofID1, score = score)))
  }

  # -- choose only the best scoring non-duplicated hits
  setorder(bl, -score)
  bl <- subset(bl, !duplicated(paste(ofID1, ofID2)))

  # -- add annotation information in
  bl[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
           start1 = sv[ofID1], start2 = sv[ofID2],
           end1 = ev[ofID1], end2 = ev[ofID2],
           chr1 = cv[ofID1], chr2 = cv[ofID2],
           ord1 = ov[ofID1], ord2 = ov[ofID2],
           og1 = ogv[ofID1], og2 = ogv[ofID2],
           scrRank1 = 1, scrRank2 = 1,
           isOg = ogv[ofID1] == ogv[ofID2])]
  if(selfOnly)
    bl <- subset(bl, chr1 == chr2)
  bl <- subset(bl, !is.na(ord1) & !is.na(ord2))

  # -- get the number of unique hits by chr (for later culling)
  bl[,`:=`(nChr1 = uniqueN(ofID1[isOg]),
           nChr2 = uniqueN(ofID2[isOg])),
     by = c("gen1","gen2","chr1","chr2")]
  bl[,`:=`(og = ifelse(isOg, og1, NA), og1 = NULL, og2 = NULL)]

  # -- optionally add the rank scores by geneID
  if(!selfOnly){
    setorder(bl, ofID1, -score)
    bl[,scrRank1 := 1:.N, by = "ofID1"]
    setorder(bl, ofID2, -score)
    bl[,scrRank2 := 1:.N, by = "ofID2"]
  }
  return(bl)
}

#' @title find_hitsInBuff
#' @description
#' \code{find_hitsInBuff} find_hitsInBuff
#' @rdname utils
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
find_hitsInBuff <- function(hits,
                            nCores,
                            synBuff){

  inBuffer <- ord1 <- ord2 <- isAnchor <- NULL

  if(!"isAnchor" %in% colnames(hits))
    stop("column name isAnchor must be in hits\n")
  splHit <- split(hits, by = c("gen1","gen2","chr1","chr2"))
  out <- rbindlist(mclapply(splHit, mc.cores = nCores, function(x){
    if(any(x$isAnchor)){
      nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synBuff))
      wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
      x[,inBuffer := 1:.N %in% wh]
      return(x)
    }
  }), fill = T)
  if(!"inBuffer" %in% colnames(out)){
    return(NULL)
  }else{
    return(subset(out, inBuffer))
  }
}

#' @title clus_dbscan
#' @description
#' \code{clus_dbscan} clus_dbscan
#' @rdname utils
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
clus_dbscan <- function(hits,
                        radius,
                        blkSize,
                        nCores){
  blkID <- ord1 <- ord2 <- chr1 <- chr2 <- NULL

  x <- split(hits, by = c("gen1","gen2","chr1","chr2"))
  x <- rbindlist(mclapply(x, mc.cores = nCores, function(y){
    y[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = radius),
                       minPts = blkSize)$cluster]
    return(y)
  }))
  x <- subset(x, blkID > 0)
  x[,blkID := sprintf("blk_%s_%s_%s", chr1, chr2, blkID)]
  return(x)
}

#' @title find_globalAnchors
#' @description
#' \code{find_globalAnchors} find_globalAnchors
#' @rdname utils
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
find_globalAnchors <- function(hits,
                               gsParam,
                               synParam,
                               type = "primary"){

  add_rbh2og <- function(hits, minRbhScore = 50){
    ofID1 <- ofID2 <- NULL
    maxScr1 <- maxScr2 <- score <- og <- NULL
    hits[,maxScr1 := max(score), by = c("blkID", "ofID1")]
    hits[,maxScr2 := max(score), by = c("blkID", "ofID2")]

    outrb <- rbindlist(lapply(split(hits, by = "blkID"), function(x){
      xog <- subset(x, !is.na(og))
      xrbh <- subset(x, maxScr1 == score & maxScr2 == score & score > minRbhScore &
                       !ofID1 %in% xog$ofID1 & !ofID2 %in% xog$ofID2)
      xrbh[,og := "RBH"]
      return(rbind(xrbh, xog))
    }))
    outrb$og[outrb$og == "RBH"] <- sprintf("RBH%s",1:sum(outrb$og == "RBH"))
    return(outrb)
  }

  nhits1 <- nSecondHits1 <- nhits2 <- nSecondHits2 <- blkSizeSecond <- NULL
  nGapsSecond <- synBuffSecond <- isMasked <- ord1 <- ord2 <- isArrayRep <- NULL
  og <- ofID1 <- score <- scrRank1 <- ofID2 <- scrRank2 <- blkID <- NULL
  isAnchor <- inBuffer <- onlyOgAnchorsSecond <- NULL
  # -- get synteny parameters in line
  sp <- data.table(synParam[1,])
  if(type != "primary")
    sp[,`:=`(nhits1 = (nhits1-1 + nSecondHits1),
             nhits2 = (nhits2-1 + nSecondHits2),
             blkSize = blkSizeSecond,
             nGaps = nGapsSecond,
             synBuff = synBuffSecond,
             onlyOgAnchors = onlyOgAnchorsSecond)]

  # -- drop masked hits and re-rank order
  hits <- subset(hits, !isMasked)
  hits[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]

  # -- drop non-array rep genes
  anch <- subset(hits, isArrayRep)

  # -- drop non-og hits, if necessary
  if(sp$onlyOgAnchors)
    anch <- subset(anch, !is.na(og))

  # -- subset to top n hits
  setorder(anch, ofID1, -score)
  anch[,scrRank1 := 1:.N, by = "ofID1"]
  setorder(anch, ofID2, -score)
  anch[,scrRank2 := 1:.N, by = "ofID2"]
  anch <- subset(anch, scrRank1 <= sp$nhits1 & scrRank2 <= sp$nhits2)

  # -- re-rank
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]

  # -- initial anchor building
  anchu <- run_mcscanx(
    hits = anch,
    gsParam = gsParam,
    blkSize = sp$blkSize,
    nGaps = sp$nGaps,
    path2mcscanx = gsParam$paths$mcscanxCall)
  anch <- subset(anch, paste(ofID1, ofID2) %in% names(anchu))

  # -- initial anchor clustering
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anch <- clus_dbscan(
    hits = anch,
    radius = sp$synBuff,
    blkSize = sp$blkSize,
    nCores = gsParam$params$nCores)
  anchv <- anch$blkID; names(anchv) <- with(anch, paste(ofID1, ofID2))

  # -- find all hits in regions in the array rep data
  anch <- subset(hits, isArrayRep)
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anch[,blkID := anchv[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(anch, !is.na(blkID)))
  anch <- get_hitsInBlks(
    blks = bc,
    hits = anch,
    nCores = gsParam$params$nCores)
  anch[,isAnchor := paste(ofID1, ofID2) %in% names(anchv)]
  anch <- subset(find_hitsInBuff(
    hits = anch,
    nCores = gsParam$params$nCores,
    synBuff = sp$synBuff), inBuffer)

  # -- add in RBHs and re cluster
  anch <- add_rbh2og(
    hits = subset(anch, !is.na(blkID) & inBuffer))
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anchu <- run_mcscanx(
    hits = anch,
    gsParam = gsParam,
    blkSize = sp$blkSize,
    nGaps = sp$nGaps,
    path2mcscanx = gsParam$paths$mcscanxCall)
  anch <- subset(anch, paste(ofID1, ofID2) %in% names(anchu))

  # -- final anchor clustering
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anch <- clus_dbscan(
    hits = anch,
    radius = sp$synBuff,
    blkSize = sp$blkSize,
    nCores = gsParam$params$nCores)
  anchv <- anch$blkID; names(anchv) <- with(anch, paste(ofID1, ofID2))
  hits[,blkID := anchv[paste(ofID1, ofID2)]]

  # -- split overlapping blocks
  anch2 <- split_overlaps(
    hits = subset(hits, !is.na(blkID)),
    dropInterleavesSmallerThan = 2,
    minPropDup = .05,
    verbose = F)
  anchv2 <- anch2$blkID; names(anchv2) <- with(anch2, paste(ofID1, ofID2))

  # -- final assignment
  hits[,blkID := anchv2[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(hits, !is.na(blkID)))
  out <- get_hitsInBlks(
    blks = bc,
    hits = hits,
    nCores = gsParam$params$nCores)
  out[,isAnchor := paste(ofID1, ofID2) %in% names(anchv)]
  out <- find_hitsInBuff(
    hits = out,
    nCores = gsParam$params$nCores,
    synBuff = sp$synBuff)
  return(out)
}

#' @title pipe_synteny
#' @description
#' \code{pipe_synteny} pipe_synteny
#' @rdname utils
#' @import data.table
#' @export
pipe_synteny <- function(gsParam,
                         gff,
                         maskHits,
                         synParam,
                         type = "primary"){

  gen1 <- gen2 <- isArrayRep <- ofID1 <- ofID2 <- isMasked <- regBuffer <- NULL
  blkIDn <- blkID <- chr1 <- chr2 <- inBuffer <- isAnchor <- lgRegionID <- NULL

  verbose <- gsParam$params$verbose
  # -- read in the blasts
  bl <- with(synParam, parse_blast4synteny(
    gsParam = gsParam,
    genome1 = genome1,
    genome2 = genome2,
    gff = gff,
    selfOnly = FALSE))
  bl <- subset(bl, !is.na(gen1) & !is.na(gen2))
  if(verbose)
    cat(sprintf("%s / ", nrow(bl)))

  # -- flag array reps
  arrayReps <- gff$ofID[gff$isArrayRep]
  bl[,isArrayRep := ofID1 %in% arrayReps & ofID2 %in% arrayReps]

  # -- optionally remove masked regions
  bl[,isMasked := FALSE]
  if(!is.null(maskHits)){
    u <- with(maskHits, paste(ofID1, ofID2))
    bl$isMasked[with(bl, paste(ofID1, ofID2)) %in% u] <- TRUE
  }

  # -- get the syntenic region hits
  regs <- find_globalAnchors(
    hits = bl,
    gsParam = gsParam,
    synParam = synParam,
    type = type)
  regs[,blkIDn := as.numeric(as.factor(blkID)), by = c("chr1", "chr2")]
  regs[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, blkIDn)]

  if(verbose)
    cat(sprintf("%s (%s) / ", sum(regs$inBuffer), uniqueN(regs$blkID, na.rm = T)))

  # -- get the syntenic blocks
  regs[,`:=`(regID = blkID, regBuffer = inBuffer, regAnchor = isAnchor)]
  blks <- finalize_blocks(
    hits = regs,
    synParam = synParam,
    gsParam = gsParam,
    verbose = T,
    type = type)

  out <- merge(
    bl,
    blks[,c("ofID1","ofID2","regID","blkID","regAnchor", "blkAnchor","blkBuffer","regBuffer")],
    by = c("ofID1", "ofID2"), all = T)

  h <- clus_dbscan(
    subset(out, regBuffer),
    radius = synParam$synBuff+1,
    blkSize = 1,
    nCores = gsParam$params$nCores)
  hb <- h$blkID;  names(hb) <- with(h, paste(ofID1, ofID2))
  out[,lgRegionID := hb[paste(ofID1, ofID2)]]
  # -- reformat and return
  if(verbose)
    cat(sprintf("%s (%s)\n", sum(blks$blkBuffer), uniqueN(blks$blkID, na.rm = T)))
  for(i in c("regAnchor", "blkAnchor","blkBuffer","regBuffer","lgRegionID"))
    out[[i]][is.na(out[[i]])] <- FALSE
  return(out)
}

#' @title pull_nonSelfReg
#' @description
#' \code{pull_nonSelfReg} pull_nonSelfReg
#' @rdname utils
#' @import data.table
#' @export
pull_nonSelfReg <- function(hits,
                            selfRegionMask,
                            nhits){
  ofID1 <- ofID2 <- ord1 <- ord2 <- NULL
  u <- with(pull_selfRegion(hits, synBuff = selfRegionMask),
            unique(paste(ofID1, ofID2)))
  x <- subset(hits, !paste(ofID1, ofID2) %in% u)
  x[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
          ord2 = frank(ord2, ties.method = "dense"))]
  return(x)
}

#' @title pull_selfRegion
#' @description
#' \code{pull_selfRegion} pull_selfRegion
#' @rdname utils
#' @import data.table
#' @export
pull_selfRegion <- function(hits,
                            synBuff){
  chr1 <- chr2 <- ord1 <- ord2 <- ofID1 <- ofID2 <- NULL
  buff <- sqrt(synBuff^2 * 2) + 1
  x <- subset(
    hits,
    chr1 == chr2 &
      abs(ord1 - ord2) <= synBuff*2)
  x[,`:=`(blkID = sprintf("blk_%s_%s_self", chr1, chr2),
          isAnchor = ofID1 == ofID2)]
  return(x)
}

#' @title get_hitsInBlks
#' @description
#' \code{get_hitsInBlks} get_hitsInBlks
#' @rdname utils
#' @import data.table
#' @importFrom parallel mclapply
#' @export
get_hitsInBlks <- function(blks, hits, nCores){
  ord1 <- ord2 <- blkID <- NULL
  splh <- split(hits, by = c("chr1","chr2","gen1","gen2"))
  splb <- split(blks, by = c("chr1","chr2","gen1","gen2"))
  splh <- splh[names(splb)]
  hits <- rbindlist(mclapply(1:nrow(blks), mc.cores = 1, function(i){
    y <- blks[i,]
    x <- splh[[with(y, paste(chr1, chr2, gen1, gen2, sep = "."))]]

    x <- subset(x, ord1 >= y$startOrd1 & ord1 <= y$endOrd1 &
                  ord2 >= y$minOrd2 & ord2 <= y$maxOrd2)
    x[,blkID := y$blkID]
    return(x)
  }))
  return(hits)
}

#' @title split_overlaps
#' @description
#' \code{split_overlaps} split_overlaps
#' @rdname utils
#' @export
split_overlaps <- function(hits,
                           dropInterleavesSmallerThan,
                           minPropDup = 0.05,
                           maxIter = 20,
                           verbose){
  ############################################################################
  ############################################################################

  ############################################################################
  # -- function to split blocks by RLEs
  split_ovlBlks <- function(blk1, blk2, hits, dropInterleavesSmallerThan){
    blkID <- ofID1 <- ofID2 <- ord1 <- rl1 <- ord2 <- rl2 <- n <- NULL
    y <- subset(hits, blkID %in% c(blk1, blk2))
    y <- subset(y, !duplicated(ofID1))
    y <- subset(y, !duplicated(ofID2))
    setkey(y, ord1)
    y[,rl1 := add_rle(blkID)]
    setkey(y, ord2)
    y[,rl2 := add_rle(blkID)]
    y[,blkID := sprintf("%s_%s%s", blkID, rl1, rl2)]
    y[,n := .N, by = blkID]
    y <- subset(y, n >= dropInterleavesSmallerThan)
    y[,`:=`(rl1 = NULL, rl2 = NULL, n = NULL)]
    return(y)
  }
  ############################################################################
  # -- function to count and find overlapping blocks
  count_ovlHits <- function(hits, minPropDup){
    blkID <- i.blkID <- blk1 <- blk2 <- r1 <- r2 <- propDup <- n <- hasPriority <- NULL
    # -- calculate block coordinates
    hitstmp <- data.table(hits)
    bc <- calc_blkCoords(hitstmp)
    bc1 <- with(bc, data.table(
      blkID = blkID, chr = chr1, start = startOrd1, end = endOrd1,
      key = c("chr","start","end")))
    bc2 <- with(bc, data.table(
      blkID = blkID, chr = chr2, start = minOrd2, end = maxOrd2,
      key = c("chr","start","end")))

    # -- find non-self overlapping regions
    fo1 <- subset(foverlaps(bc1, bc1), blkID != i.blkID)
    fo2 <- subset(foverlaps(bc2, bc2), blkID != i.blkID)

    # -- add overlapping info into an output object
    out <- data.table(
      blk1 = c(fo1$blkID, fo2$blkID),
      blk2 = c(fo1$i.blkID, fo2$i.blkID))

    # -- choose the non-repeated overlaps
    u <- unique(c(out$blk1, out$blk2))
    un <- order(u); names(un) <- u
    out[,`:=`(r1 = un[blk1], r2 = un[blk2])]
    out[,u := ifelse(r1 > r2,  paste(blk1, blk2), paste(blk2, blk1))]
    out <- subset(out, !duplicated(u))
    out[,`:=`(r1 = NULL, r2 = NULL, u = NULL)]

    # -- count the number of duplicated hits between each pair of blocks
    if(nrow(out) > 0){
      tmp <- split(subset(hitstmp, blkID %in% unlist(out)), by = "blkID")
      out[,propDup := sapply(1:nrow(out), function(j){
        tb1 <- tmp[[out$blk1[j]]]
        tb2 <- tmp[[out$blk2[j]]]
        tot <- uniqueN(c(tb1$ofID1, tb1$ofID2, tb2$ofID1, tb2$ofID2))
        dup <- sum(tb1$ofID1 %in% tb2$ofID1, tb1$ofID2 %in% tb2$ofID2,
                   tb2$ofID1 %in% tb1$ofID1, tb2$ofID2 %in% tb1$ofID2)
        return(dup/tot)
      })]
      out <- subset(out, propDup <= minPropDup)
      if(nrow(out) > 0){
        # -- get the total number of hits in each block
        out[,`:=`(n1 = sapply(tmp[blk1], nrow),
                  n2 = sapply(tmp[blk2], nrow))]

        # -- prioritize by the pairs with the largest small block
        out[,n := apply(.SD, 1, min), .SDcols = c("n1", "n2")]

        setorder(out, -n)
        out[,n := NULL]
        v <- matrix(duplicated(as.character(unlist(t(out[,1:2])))), ncol = 2, byrow = T)
        out[,hasPriority := !apply(v, 1, any)]
        return(out)
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  }
  ############################################################################
  ############################################################################
  # -- iterate through trying to split overlapping blocks
  hasPriority <- blkID <- NULL
  anchs <- data.table(hits)
  nov <- 1
  iter <- 0
  while(nov > 0 && iter <= maxIter){
    iter <- iter + 1
    if(verbose)
      cat(sprintf("\titer %s: %s hits / %s blks\n",
                  iter, nrow(anchs), uniqueN(anchs$blkID)))
    ovlBlks <- count_ovlHits(anchs, minPropDup = minPropDup)
    if(!is.null(ovlBlks)){
      ovlBlks <- subset(ovlBlks, hasPriority)
      nov <- nrow(ovlBlks)
      if(is.null(nov))
        nov <- 0
      if(nov > 0){
        splAnch <- rbindlist(lapply(1:nrow(ovlBlks), function(j)
          split_ovlBlks(
            blk1 = ovlBlks$blk1[j],
            blk2 = ovlBlks$blk2[j],
            hits = anchs,
            dropInterleavesSmallerThan)))
        if(!is.null(splAnch)){
          anchs <- rbind(
            splAnch,
            subset(anchs, !blkID %in% unique(unlist(ovlBlks[,c(1:2)]))))
        }else{
          nov <- 0
        }
      }else{
        nov <- 0
      }
    }else{
      nov <- 0
    }
  }
  if(verbose)
    cat("\nfound no additional overlapping blocks\n")
  return(anchs)
}

#' @title add number of anchors to gff
#' @description
#' \code{annotate_gff} add number of anchors to gff
#' @rdname utils
#' @importFrom parallel mclapply
#' @export
annotate_gff <- function(gsParam,
                         genomeIDs = NULL,
                         overwrite = FALSE){

  genome <- ord  <- arrayID <- isArrayRep <- med <- pepLen <- rnk <- ofID <- NULL
  medbp <- start <- dist2med <- dist2bp <- og <- n <- chr <- globOG <- NULL

  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(file.exists(gffFile) & !overwrite){
    tmp <- fread(gffFile, na.strings = c("-", "NA", ""), showProgress = F)
    if(all(c("pepLen", "globOG", "arrayID","isArrayRep") %in% colnames(tmp)))
      stop("annotated gff file exists and !overwrite, so not running ...\n")
  }
  if(!is.data.table(gsParam$params$synteny))
    stop("Must run set_syntenyParams first!\n")

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs

  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  synBuff <- max(gsParam$params$synteny$synBuff)
  if(verbose)
    cat("Loading annotations ...\n")
  # -- get paths to the orthofinder run
  if(is.na(gsParam$paths$orthogroupsDir)){
    if(verbose)
      cat("\tIndexing location of orthofinder results ... ")
    gsParam <- find_orthofinderResults(gsParam)
    if(verbose)
      cat("Done!\n")
  }

  ##############################################################################
  # 1. Load the gff and add global metadata
  ##############################################################################
  # -- read in the gff
  if(verbose)
    cat("\tReading the gffs ... ")
  gff <- read_gff(gsParam$paths$gff)
  gff <- add_ofID2gff(gff, gsParam$paths$blastDir)
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)

  # -- add peptide length
  if(verbose)
    cat("Done!\n\tPulling gene lengths ... ")
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)

  # -- add global orthogroupsto the gff
  if(verbose)
    cat("Done!\n\tParsing global orthogroups ... ")
  ogs <- parse_ogs(gsParam)
  gff <- merge(gff, ogs, by = c("genome","id"), all.x = T)
  gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
  setnames(gff, "ogID", "globOG")

  ##############################################################################
  # 2. Define collinear arrays
  ##############################################################################
  # -- build arrays from method in gsParam
  if(verbose)
    cat("Done!\nDefining collinear orthogroup arrays ... \n")
  gff <- add_arrays2gff(gsParam = gsParam, gff = gff)
  gff[,arrayID := paste(arrayID, og)]
  gff[,n := .N, by = "arrayID"]
  gff$arrayID[gff$n == 1] <- NA
  gff[,`:=`(n = NULL, isArrayRep = TRUE, og = NULL)]

  # -- choose the array reps
  if(verbose)
    cat("\tChoosing array representative genes ... ")
  gffa <- subset(gff, !is.na(arrayID))
  setkey(gffa, genome, ord)
  gffn <- subset(gff, is.na(arrayID))
  gffa[,arrayID := sprintf(
    "%s_%s_%s",
    genome, chr, as.numeric(factor(arrayID, levels = unique(arrayID))))]
  gffa[,med := as.numeric(median(ord)), by = "arrayID"]
  gffa[,medbp := as.numeric(median(start)), by = "arrayID"]
  gffa[,`:=`(dist2med = abs(med - ord),
             dist2bp = abs(medbp - start))]
  setorder(gffa, arrayID, dist2med, dist2bp, -pepLen)
  gffa[,rnk := 1:.N, by = "arrayID"]
  gffa[,`:=`(isArrayRep = rnk == 1, rnk = NULL, dist2med = NULL, dist2bp = NULL,
             med = NULL, medbp = NULL)]
  gff <- rbind(gffa, gffn)
  setkey(gff, genome, ord)

  if(verbose)
    cat("Done!\n")
  gff[,`:=`(synOG = NA, inBlkOG = NA, combOG = NA, og = globOG, refCoord = NA)]
  fwrite(gff, file = gffFile, sep = "\t", quote = F, showProgress = F)
  return(gsParam)
}

#' @title add array representative to a gff object
#' @description
#' \code{add_arrayRep2gff} choose most central gene by orthogroup
#' @rdname utils
#' @importFrom parallel mclapply
#' @export
add_arrayRep2gff <- function(gff,
                             gsParam){
  synArr <- chr <- NULL
  # -- count peptides
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)
  genomeIDs <- unique(gff$genome)
  nCores <- gsParam$params$nCores

  # -- count number of orthologs
  di <- dir.exists(gsParam$paths$orthologuesDir)
  dl <- length(list.files(gsParam$paths$orthologuesDir)) > 1

  nGenome <- nGenes <- gen2 <- id1 <- gen1 <- genome <- id <- og <- NULL
  if(di && dl){
    ogcnt <- rbindlist(mclapply(genomeIDs, mc.cores = nCores, function(i){
      ogs <- parse_orthologues(gsParam = gsParam, refGenome = i)
      ogn <- ogs[,list(nGenome = uniqueN(gen2),
                       nGenes = .N), by = c("gen1","id1")]
      return(ogn)
    }))
    setorder(ogcnt, -nGenome, -nGenes)
    ogcnt <- subset(ogcnt, !duplicated(paste(gen1, id1)))
    nog <- ogcnt$nGenome; ng <- ogcnt$nGenes
    names(ng) <- names(nog) <- with(ogcnt, paste(gen1, id1))
    gff[,`:=`(nGenomeOrthologs = nog[paste(genome, id)],
              nTotalOrthologs = ng[paste(genome, id)])]
    gff$nGenomeOrthologs[is.na(gff$nGenomeOrthologs)] <- 0
    gff$nTotalOrthologs[is.na(gff$nTotalOrthologs)] <- 0
  }else{
    gff[,`:=`(nGenomeOrthologs = 0,
              nTotalOrthologs = 0)]
  }

  # -- split into single and multiple member arrays
  gffi <- data.table(gff)
  d2h <- gsParam$params$maxDistBtwPgHits
  gff[,synArr := as.integer(as.factor(paste(genome, chr, og)))]
  gff[,nog := .N, by = "synArr"]
  g1 <- subset(gff, nog == 1)
  g1[,nog := NULL]
  g2 <- subset(gff, nog > 1)

  # -- calculate the maximum distance between genes in an array
  maxJump <- ord <- NULL
  g2[,maxJump := max(diff(ord[order(ord)])), by = "synArr"]
  g2r <- subset(g2, maxJump > d2h)
  g2 <- subset(g2, maxJump <= d2h)

  # -- cluster genes in arrays with big jumps
  clus <- NULL
  if(nrow(g2r) > 1){
    g2[,maxJump := NULL]
    g2r[,clus := dbscan(frNN(cbind(ord, ord), eps = d2h), minPts = 0)$cluster,
        by = "synArr"]
    g2r[,synArr := paste(synArr, clus)]
    g2 <- rbind(g2,  g2r[,colnames(g2),with = F])
  }

  # -- calculate distance to the median
  dist2median <- nGenomeOrthologs <- pepLen <- nTotalOrthologs <- ofID <- NULL
  g2[,dist2median := abs(as.numeric(median(ord, na.rm = T)) - ord),
     by = c("synArr","genome","chr")]

  # -- order and rank genes, choosing representatives for each array
  setorder(g2, genome, chr, synArr, -nGenomeOrthologs, -nTotalOrthologs,
           dist2median, -pepLen, ord)
  arep <- rbind(g1, g2[,colnames(g1), with = F])
  sar <- as.numeric(as.factor(arep$synArr)); names(sar) <- arep$ofID
  arep <- subset(arep, !duplicated(synArr))
  gffi[,`:=`(synArray = sar[ofID],
             isArrayRep = ofID %in% arep$ofID)]
  return(gffi)
}

#' @title add_arrays2gff
#' @description
#' \code{add_arrays2gff} add_arrays2gff
#' @rdname utils
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
add_arrays2gff <- function(gsParam,
                           gff){
  arrayID <- genome <- chr <- globOG <- n <- rng <- ord <- clus <- og <- NULL
  collinearOG <- NULL

  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  synBuff <- max(gsParam$params$synteny$synBuff)

  # -- make global arrays from orthogroups
  gff[,arrayID := sprintf("%s_%s_%s", genome, chr, globOG)]
  gff[,n := .N, by = "arrayID"]

  # -- combine 1x ogs with ogs in regions < synBuff
  g1 <- subset(gff, n == 1)
  g2 <- subset(gff, n > 1)
  g2[,rng := diff(range(ord)),  by = "arrayID"]
  g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
  g2 <- subset(g2, rng > synBuff)

  # -- combine above with ogs without max gap < synBuff
  g2[,rng := max(diff(ord[order(ord)])), by = "arrayID"]
  g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
  g2 <- subset(g2, rng > synBuff)

  # -- split ogs with gaps
  g2[,clus := dbscan(frNN(cbind(ord, ord), eps = synBuff), minPts = 1)$cluster,
     by = "arrayID"]
  g2[,arrayID := sprintf("%s_%s", arrayID, clus)]
  gff <- rbind(g1, g2[,colnames(g1), with = F])

  # -- NA out arrays with just one member
  gff[,n := .N, by = "arrayID"]
  gff[,og := arrayID]
  gff$arrayID[gff$n == 1] <- NA
  gff[,n := NULL]
  setkey(gff, genome, ord)

  # -- print updates and number of global orthogroups
  if(verbose){
    cat("\tUsing collinear orthogroups for array identity:\n")
    nu <- lapply(split(subset(gff, !is.na(arrayID)), by = "genome"), function(x)
      cat(sprintf("\t%s: %s genes in %s collinear arrays\n",
                  x$genome[1], nrow(x), uniqueN(x$arrayID))))
  }

  return(gff)
}

#' @title add_synOg2gff
#' @description
#' \code{add_synOg2gff} add_synOg2gff
#' @rdname utils
#' @export
add_synOg2gff <- function(gff,
                          hits = NULL,
                          gsParam,
                          genomeIDs,
                          allowRBHinOg,
                          useBlks){
  blkBuffer <- isSelf <- ofID1 <- ofID2 <- og <- ofID <- synOG <- blkID <- NULL
  # -- find hits files
  if(is.null(hits)){
    eg <- CJ(genomeIDs, genomeIDs)
    fs <- file.path(gsParam$paths$results,
                    sprintf("%s_%s_synHits.txt.gz",eg[[1]], eg[[2]]))
    fs <- fs[file.exists(fs)]

    # -- read all hits
    nCores <- gsParam$params$nCores
    hts <- rbindlist(mclapply(fs, mc.cores = nCores, function(i){
      if(useBlks){
        x <- fread(
          i, select = c("ofID1","ofID2","og","blkBuffer","blkAnchor","blkID"),
          na.strings = c("","NA"))
      }else{
        x <- fread(
          i, select = c("ofID1","ofID2","og","regBuffer","regAnchor","regID"),
          na.strings = c("","NA"),
          col.names = c("ofID1","ofID2","og","blkBuffer","blkAnchor","blkID"))
      }
      x <- subset(x, blkBuffer)
      x[,isSelf := any(ofID1 == ofID2), by = "blkID"]
      x <- subset(x, !isSelf & !is.na(og))
      if(!allowRBHinOg)
        x <- subset(x, !grepl("RBH", og))
      return(x[,c("ofID1", "ofID2", "blkAnchor", "blkID")])
    }))
  }else{
    hts <- data.table(hits)
    hts <- subset(hts, !is.na(og))
    if(!allowRBHinOg)
      hts <- subset(hts, !grepl("RBH", og))
    if(useBlks){
      hts <- hts[,c("ofID1", "ofID2", "blkBuffer", "blkAnchor", "blkID")]
    }else{
      hts <- with(hts, data.table(
        ofID1 = ofID1, ofID2 = ofID2,
        blkBuffer = regBuffer, blkAnchor = regAnchor, blkID = regID))
    }
    hts <- subset(hts, blkBuffer)
    hts[,blkBuffer := NULL]
  }

  # -- convert to syntenic orthogroups
  ic <- with(subset(hts, !is.na(blkID)), clus_igraph(
    id1 = c(ofID1, ofID2), id2 = c(ofID2, ofID1)))
  gff[,synOG := ic[ofID]]
  nmis <- sum(is.na(gff$synOG))
  mol <- max(gff$synOG, na.rm = T)
  gff$synOG[is.na(gff$synOG)] <- (mol + 1):(mol + nmis)
  gff[,synOG := as.integer(factor(synOG, levels = unique(synOG)))]
  return(gff)
}


#' @title combine_inblkSynOG
#' @description
#' \code{combine_inblkSynOG} combine_inblkSynOG
#' @rdname utils
#' @import data.table
#' @export
combine_inblkSynOG <- function(genomeIDs,
                               gff,
                               gsParam){

  ofID <- ofID1 <- ofID2 <- clus <- combOG <- inBlkOG <- synOG <- NULL
  if(gsParam$params$verbose)
    cat("Combining synteny-constrained and inblock orthogroups ...\n")

  if(gsParam$params$verbose)
    cat(sprintf("\tsyn OGs: %s, inblk OGs: %s",
                uniqueN(gff$synOG, na.rm = T), uniqueN(gff$inBlkOG, na.rm = T)))
  if(all(is.na(gff$inBlkOG)))
    gff[,inBlkOG := synOG]
  inblk <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "inBlkOG"]
  syn <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "synOG"]
  u <- with(inblk, paste(ofID1, ofID2))
  syn <- subset(syn, !paste(ofID1, ofID2) %in% u & ofID1 != ofID2)
  tmp <- rbind(syn[,c("ofID1", "ofID2")], inblk[,c("ofID1", "ofID2")])
  tmp[,clus := clus_igraph(ofID1, ofID2)]
  ov <- with(tmp, c(clus, clus)); names(ov) <- with(tmp, c(ofID1, ofID2))
  gff[,combOG := ov[ofID]]
  nmis <- sum(is.na(gff$combOG))
  mol <- max(gff$combOG, na.rm = T)
  gff$combOG[is.na(gff$combOG)] <- (mol + 1):(mol + nmis)
  gff[,combOG := as.integer(factor(combOG, levels = unique(combOG)))]

  if(gsParam$params$verbose)
    cat(sprintf(", combined OGs: %s\n",
                uniqueN(gff$combOG)))
  return(gff)
}

#' @title pull_nonSynOrthologs
#' @description
#' \code{pull_nonSynOrthologs} pull_nonSynOrthologs
#' @rdname utils
#' @import data.table
#' @export
pull_nonSynOrthologs <- function(gsParam,
                                 gff){

  gen1 <- og2 <- og1 <- orthIDs <- ofID <- id2 <- gen2 <- id1 <- NULL
  idv <- gff$ofID; names(idv) <- with(gff, paste(genome, id))
  ogv <- gff$combOG; names(ogv) <- gff$ofID

  orths <- rbindlist(lapply(unique(gff$genome), function(i){
    x <- parse_orthologues(
      gsParam = gsParam,
      refGenome = i,
      nCores = gsParam$params$nCores)
    x[,`:=`(ofID = idv[paste(gen1, id1)], orthIDs = idv[paste(gen2, id2)])]
    x[,`:=`(og1 = ogv[ofID], og2 = ogv[orthIDs])]
    x <- subset(x, og1 != og2)
    return(x)
  }))
  return(orths[,c("ofID", "orthIDs", "orthID")])
}

#' @title pull_blkAnchors
#' @description
#' \code{pull_blkAnchors} pull_blkAnchors
#' @rdname utils
#' @import data.table
#' @export
pull_blkAnchors <- function(gsParam,
                            gff,
                            refGenome){
  isSelf <- blkAnchor <- ofID2 <- ofID1 <- ord1 <- ord2 <- g1 <- g2 <- NULL
  genomeIDs <- unique(gff$genome)

  # -- get the hit files
  pfs <- CJ(g1 = genomeIDs, g2 = genomeIDs)
  pfs[,fs := file.path(gsParam$paths$results,
                       sprintf("%s_%s_synHits.txt.gz", g1, g2))]
  pfs <- subset(pfs, file.exists(fs))
  fs <- pfs$fs

  # -- block coords from hits
  out <- rbindlist(mclapply(fs, mc.cores = gsParam$params$nCores, function(i){
    x <- fread(i,
               select = c("ofID1", "ofID2","blkID","blkAnchor","gen1"),
               showProgress = F,
               na.strings = c("NA", "-", ""))
    if(x$gen1 != refGenome){
      setnames(x, c("ofID2", "ofID1","blkID","blkAnchor","gen1"))
      x <- x[,c("ofID1", "ofID2","blkID","blkAnchor","gen1")]
    }
    x[,isSelf := any(ofID1 %in% ofID2), by = "blkID"]
    return(subset(x, blkAnchor & !isSelf)[,1:3])
  }))

  gv <- gff$genome; ov <- gff$ord; cv <- gff$chr
  names(gv) <- names(ov) <- names(cv) <- gff$ofID
  out[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
            chr1 = cv[ofID1], chr2 = cv[ofID2],
            ord1 = ov[ofID1], ord2 = ov[ofID2])]
  bc <- out[,list(start1 = min(ord1, na.rm = T), end1 = max(ord1, na.rm = T),
                  start2 = min(ord2, na.rm = T), end2 = max(ord2, na.rm = T)),
            by = c("gen1","gen2", "chr1", "chr2", "blkID")]
  return(list(anchors = out, coords = bc))
}
