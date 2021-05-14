#' @title Generic internal functions used by genespace
#' @description
#' \code{utils} Convience functions for genespace, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
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
#' @param gffFiles named vector of file paths pointing to formatted gff3-like
#' files
#' @param gff data.table parsed from a gff3-formatted annotation
#' @param blastDir path to the directory containing the initial
#' blast/orthofinder run.
#' @param genomeIDs character vector of genomeIDs
#' @param blFile file.path pointing to a blast8-formatted text file
#' @param ofID1 column of orthofinderIDs
#' @param ofID2 column of orthofinderIDs
#' @param onlyIDScore logical, should only three columns, ofID1, ofID2 and score
#' be returned?
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run.
#'
#' @note \code{utils} is a generic name for the functions documented.
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

#' @title Check if a fasta is peptide
#' @description
#' \code{check_isPeptideFasta} Load the first few entries of a fasta and
#' make sure it a peptide. Faster than reading in the whole thing.
#' @rdname utils
#' @importFrom parallel detectCores
#' @export
check_isPeptideFasta <- function(x){
  tmp <- readLines(x, 1000)
  tmp <- tmp[!(grepl("^>",tmp) |
                 grepl("#",tmp, fixed = T) |
                 grepl("*",tmp, fixed = T))][1]
  isPeptide <- nchar(gsub("A|T|G|C|N|X", "", toupper(tmp))) > 0
  if(!isPeptide)
    stop(x,"appears to not be a peptide fasta file\n")
}

#' @title Check number of cores specified
#' @description
#' \code{check_nCores} Ensure that user-specified number of parallel processes
#' is within the range of the system.
#' @rdname utils
#' @importFrom parallel detectCores
#' @export
check_nCores <- function(nCores){
  if(length(nCores > 1)) nCores <- nCores[1]
  nCores <- as.integer(nCores)
  if(is.null(nCores)) nCores <- detectCores()/2
  if(is.na(nCores)) nCores <- detectCores()/2
  if(nCores < 1) nCores <- 1
  if(nCores > detectCores()){
    warning(sprintf(
      "user specified %s cores, but only %s available\n\tSetting nCores to %s",
      nCores, detectCores(), detectCores()))
    nCores <- detectCores()
  }
  return(nCores)
}

#' @title Check that orthofinder is installed
#' @description
#' \code{check_orthofinderInstall} Check that orthofinder is installed
#' @rdname utils
#' @export
check_orthofinderInstall <- function(path){
  grepl("OrthoFinder",
        system(paste(path, "-h"),
               intern = T)[2])
}

#' @title Check that diamond is installed
#' @description
#' \code{check_diamondInstall} Check that diamond is installed
#' @rdname utils
#' @export
check_diamondInstall <- function(path){
  grepl("diamond",
        system(paste(path, "help"),
               intern = T,ignore.stderr = T)[1])
}

#' @title Check that MCScanX_h is installed
#' @description
#' \code{check_MCScanXhInstall} Check that MCScanX_h is installed
#' @rdname utils
#' @export
check_MCScanXhInstall <- function(path){
  suppressWarnings(grepl("prefix_fn",
                         system(paste(path, "-h"),
                                intern = T,ignore.stderr = T)[1]))
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
get_nAA <- function(path){
  pepF <- list.files(path, pattern = "^Species", full.names = T)
  pepF <- pepF[grep(".fa$", pepF)]
  peps <- rbindlist(lapply(pepF, function(x){
    y <- readAAStringSet(x)
    return(data.table(ofID = names(y),
                      nAA = width(y)))
  }))
  return(peps)
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

#' @title add orthofinder ID to a gff object
#' @description
#' \code{add_ofID2gff} read the orthofinder species and gene IDs and merge
#' these with the gff-like data.table
#' @rdname utils
#' @export
add_ofID2gff <- function(gff, blastDir){
  id <- ofID <- genomeNum <- genome <- NULL
  specIDs <- read_orthofinderSpeciesIDs(blastDir)
  gv <- names(specIDs); names(gv) <- as.character(specIDs)
  seqIDs <- read_orthofinderSequenceIDs(blastDir)
  seqIDs[,genome :=  gv[as.character(genomeNum)]]
  idv <- seqIDs$ofID; names(idv) <- with(seqIDs, paste(genome, id))
  gff[,ofID := idv[paste(genome, id)]]
  return(gff)
}
