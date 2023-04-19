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
#' @param id1 character, first id
#' @param id2 character, second id
#' @param y numeric, y values
#' @param hits data.table of syntenic hits
#' @param col color
#' @param scale1toMean logical, should values be scaled to 1?
#' @param filepath file.path
#' @param alpha numeric, transparency
#' @param gsParam genespace parameters, see init_genespace.
#' @param refGenome character string specifying the reference genome
#' @param bed data.table containing the combined bed object
#' @param synBuff see init_genespace.
#' @param maxIter integer, the maximum number of iterations to use
#' @param reorder logical, should the gene rank position be re-ordered?
#' @param isAnchor logical, is a hit an anchor?
#' @param radius numeric, the 2d search radius.
#' @param blkID vector of block IDs
#' @param ladderize logical, should the tree be ladderized?
#' @param treFile file.path to the tree file.
#' @param verbose logical, should updates be printed to the console?
#' @param ... additional parameters passed on to other functions
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'

#' @title startup messages
#' @description
#' \code{.onAttach} startup messages
#' @rdname utils
#' @export
.onAttach <- function(...) {
  packageStartupMessage(paste(strwrap(
    "GENESPACE v1.2.3: synteny and orthology constrained
    comparative genomics\n",
    indent = 0, exdent = 8), collapse = "\n"))
}

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

#' @title check that a file.path is valid
#' @description
#' \code{check_filePathParam} QC of user-specified parameter
#' @rdname utils
#' @export
check_filePathParam <- function(filepath){
  if(is.null(filepath))
    filepath <- NA
  filepath <- path.expand(check_character(
    x = filepath, default = NA, na.rm = FALSE))
  chk <- dir.exists(filepath) || file.exists(filepath)
  if(!chk)
    filepath <- NA
  return(filepath)
}

#' @title Check if a sequence is only DNA
#' @description
#' \code{check_onlyDNA} QC to ensure the peptides are actually peptides
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


#' @title read peptide fasta
#' @description
#' \code{read_aaFasta} read fasta-formatted peptide sequences
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

#' @title read bed file
#' @description
#' \code{read_bed} read and check a raw bed file with four columns.
#' @rdname utils
#' @import data.table
#' @importFrom stats complete.cases
#' @export
read_bed <- function(filepath){
  chk <- tryCatch(
    {
      suppressWarnings(suppressMessages(fread(
        filepath, verbose = FALSE, showProgress = FALSE, select = 1:4,
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

#' @title left justify
#' @description
#' \code{align_charLeft} for a vector of character strings, add " " to the right
#' side so they all align to the left when printed
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

#' @title right justify
#' @description
#' \code{align_charRight} for a vector of character strings, add " " to the left
#' side so they all align to the right when printed
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
read_orthofinderSpeciesIDs <- function(filepath){
  si <- fread(
    filepath,
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
read_orthofinderSequenceIDs <- function(filepath){
  V1 <- ofID <- NULL
  gi <- data.table(readLines(filepath))
  gi[,c("ofID", "id") := tstrsplit(V1, ": ")]
  gi[,c("genomeNum", "geneNum") := tstrsplit(ofID, "_", type.convert = T)]
  gi[,V1 := NULL]
  return(gi)
}

#' @title Count the number of sequences in a fasta file
#' @description
#' \code{get_nSeqs} Counts the number of lines with ">" in a file. If the output
#' is not convertible to an integer, returns NA.
#' @rdname utils
#' @export
get_nSeqs <- function(filepath){
  nseq <- system2("grep", "-c '>'", filepath, stdout = TRUE, stderr = TRUE)
  suppressWarnings(nseq <- as.integer(nseq))
  return(nseq)
}


#' @title check genespace-formatted annotations
#' @description
#' \code{check_annotFiles} ensure the annotations match correctly
#' @rdname utils
#' @import data.table
#' @export
check_annotFiles <- function(filepath, genomeIDs){

  # -- make sure the working directory is OK
  wd <- check_filePathParam(filepath)
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
    if(any(grepl(":", names(aa), fixed = T)))
      stop("some genes have ':' in the names. This string cannot be in gene names as it is used as the dictionary separator by OrthoFinder")
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
check_MCScanXhInstall <- function(filepath){
  path <- check_filePathParam(filepath)
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

#' @title parse OGs
#' @description
#' \code{parse_ogs} read and parse orthofinder orthogroups.tsv files
#' @rdname utils
#' @import data.table
#' @export
parse_ogs <- function(filepath, genomeIDs){
  id <- genome <- Orthogroup <- NULL
  tmp <- fread(filepath, showProgress = FALSE, header = TRUE, check.names = FALSE)
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

#' @title parse HOGs
#' @description
#' \code{parse_hogs} read and parse orthofinder phylogenetically hierarchical
#' orthogroup (N0.tsv) files
#' @rdname utils
#' @import data.table
#' @export
parse_hogs <- function(filepath){
  id <- genome <- HOG <- hogID <- OG <- orthofinderInternalData_XXX_OG <-
    orthofinderInternalData_XXX_hogID <- orthofinderInternalData_XXX_HOG <- NULL
  d <- fread(filepath, showProgress = FALSE, header = TRUE, check.names = FALSE)
  setnames(d, 1:3, sprintf("orthofinderInternalData_XXX_%s", colnames(d)[1:3]))
  sd <- colnames(d)[-(1:3)]
  d[,orthofinderInternalData_XXX_hogID := paste(
    orthofinderInternalData_XXX_HOG, orthofinderInternalData_XXX_OG)]
  m <- melt(
    d, id.vars = "orthofinderInternalData_XXX_hogID",
    measure.vars = sd,
    value.name = "id", variable.name = "genome")
  m <- subset(m, id != "")
  tmp <- m[,list(id = strsplit(id, ",")[[1]]),
           by = c("orthofinderInternalData_XXX_hogID", "genome")]
  tmp[,id := trimws(id)]
  setnames(tmp, "orthofinderInternalData_XXX_hogID", "hogID")
  return(tmp)
}

#' @title parse orthologs
#' @description
#' \code{parse_orthologues} read and parse orthofinder orthologs
#' @rdname utils
#' @import data.table
#' @export
parse_orthologues <- function(filepath){
  orthID <- id1 <- id2 <- NULL
  x <- fread(filepath, showProgress = F, header = TRUE, check.names = FALSE)
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

#' @title round to the nearest integer
#' @description
#' \code{round_toInteger} flexible rounding to any integer.
#' @rdname utils
#' @import data.table
#' @export
round_toInteger <- function(x, to){
  round(x/to, 0) * to
}

#' @title convert vector to RLE
#' @description
#' \code{add_rle} run-length equivalent conversion, either as the length of the
#' runs or the unique run ids.
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

#' @title genespace colors
#' @description
#' \code{gs_colors} get a set of colors from the genespace palette
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

#' @title clustering via igraph
#' @description
#' \code{clus_igraph} cluster connected subgraphs from pairwise observations
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

#' @title read combBed file
#' @description
#' \code{read_combBed} ensures consistent combBed IO
#' @rdname utils
#' @export
read_combBed <- function(filepath){
  bedNames <- c(
    "chr", "start", "end", "id", "ofID", "pepLen", "ord", "genome", "arrayID",
    "isArrayRep", "globOG", "globHOG",  "noAnchor", "og")
  cl <- c("numeric", "character", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, bedNames))
    stop("combBed.txt file is malformed\n")
  bedClass <- cl[c(2,1,1,2,2,1,1,2,2,3,2,2,3,2)]
  bed <- fread(
    filepath, na.strings = c("", "NA"), select = bedNames,
    colClasses = bedClass, showProgress = F)
  return(bed)
}

#' @title write combBed file
#' @description
#' \code{write_combBed} ensures consistent combBed IO
#' @rdname utils
#' @export
write_combBed <- function(x, filepath){
  bedNames <- c(
    "chr", "start", "end", "id", "ofID", "pepLen", "ord", "genome", "arrayID",
    "isArrayRep", "globOG", "globHOG", "noAnchor", "og")
  if(!all(bedNames %in% colnames(x)))
    stop("bed data.table is malformed\n")
  x <- x[,bedNames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}


#' @title read allBlast file
#' @description
#' \code{read_allBlast} ensures consistent allBlast IO
#' @rdname utils
#' @export
read_allBlast <- function(filepath, ...){
  hnames <- c(
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "pid", "length", "mismatches", "gapopenings", "queryStart", "queryEnd",
    "subjectStart", "subjectEnd", "Evalue", "bitScore", "sameOG", "noAnchor",
    "isAnchor", "isSyntenic", "regID", "blkID")
  cl <- c("numeric", "character", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, hnames))
    stop("allBlast.txt file is malformed\n")
  hc <- cl[c(2,2,1,1,2,1,2,3,2,2,1,1,2,1,2,3,1,1,1,1,1,1,1,1,1,1,3,3,3,3,2,2)]
  hits <- fread(
    filepath, na.strings = c("", "NA"), select = hnames,
    colClasses = hc, showProgress = F, ...)
  return(hits)
}

#' @title write allBlast file
#' @description
#' \code{write_allBlast} ensures consistent allBlast IO
#' @rdname utils
#' @export
write_allBlast <- function(x, filepath){
  hnames <- c(
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "pid", "length", "mismatches", "gapopenings", "queryStart", "queryEnd",
    "subjectStart", "subjectEnd", "Evalue", "bitScore", "sameOG", "noAnchor",
    "isAnchor", "isSyntenic", "regID", "blkID")
  if(!all(hnames %in% colnames(x)))
    stop("allBlast data.table is malformed\n")
  x <- x[,hnames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}

#' @title read synHits file
#' @description
#' \code{read_synHits} ensures consistent synHit IO
#' @rdname utils
#' @export
read_synHits <- function(filepath, ...){
  hnames <- c(
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "bitScore", "sameOG", "noAnchor", "isAnchor", "regID", "blkID")
  cl <- c("numeric", "character", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, hnames))
    stop("synhits.txt file is malformed\n")
  hc <- cl[c(2,2,1,1,2,1,2,3,2,2,1,1,2,1,2,3,1,3,3,3,2,2)]
  hits <- fread(
    filepath, na.strings = c("", "NA"), select = hnames,
    colClasses = hc, showProgress = F, ...)
  return(hits)
}

#' @title write synHits file
#' @description
#' \code{write_synHits} ensures consistent synHit IO
#' @rdname utils
#' @export
write_synHits <- function(x, filepath){
  hnames <- c(
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "bitScore", "sameOG", "noAnchor", "isAnchor", "regID", "blkID")
  if(!all(hnames %in% colnames(x)))
    stop("synHits data.table is malformed\n")
  x <- x[,hnames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}


#' @title add transparency
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col, alpha = 1){

  if(missing(col) || !all(are_colors(col)))
    stop("Colors are misspecified\n")
  if(length(alpha) != 1 || alpha > 1 || alpha < 0)
    stop("alpha is misspecified\n")

  return(apply(sapply(col, col2rgb)/255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha)))
}

#' @title read interpolated syntenic position
#' @description
#' \code{read_intSynPos} utility to read interpolated syntenic position files
#' @rdname utils
#' @export
read_intSynPos <- function(filepath){
  hnames <- c("blkID", "genome1", "chr1", "genome2", "chr2", "ofID1", "ofID2",
              "ord1", "ord2")
  hclass <- c("character", "character", "character", "character", "character",
              "character", "character", "numeric", "numeric")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, hnames))
    stop("integratedSynPos file is malformed\n")
  synpos <- fread(
    filepath, na.strings = c("", "NA"), select = hnames,
    colClasses = hclass, showProgress = F)
  return(synpos)
}

#' @title write interpolated syntenic position
#' @description
#' \code{write_intSynPos} utility to write interpolated syntenic position files
#' @rdname utils
#' @export
write_intSynPos <- function(x, filepath){
  hnames <- c("blkID", "genome1", "chr1", "genome2", "chr2", "ofID1", "ofID2",
              "ord1", "ord2")
  x <- x[,hnames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}

#' @title Get the version of the orthofinder install
#' @description
#' \code{get_orthofinderVersion} Checks that orthofinder is installed and if so,
#' returns the installed version.
#' @rdname utils
#' @export
get_orthofinderVersion <- function(filepath){
  # -- check if orthofinder is callable
  path <- path.expand(filepath)
  wh <- Sys.which(as.character(path))
  isThere <- basename(wh) == "orthofinder"
  # -- if orthofinder is installed, check version
  if(isThere){
    ver <- system2(path, "-h", stdout = TRUE)[2]
    ver <- strsplit(ver, " ")[[1]][3]
    vern <- strsplit(ver, ".", fixed = T)[[1]]
    vern <- as.numeric(sprintf("%s.%s%s", vern[1], vern[2], vern[3]))
    return(vern)
  }else{
    return(NA)
  }
}

#' @title Get the version of the DIAMOND install
#' @description
#' \code{get_diamondVersion} Checks that DIAMOND is installed and if so,
#' returns the installed version.
#' @rdname utils
#' @export
get_diamondVersion <- function(filepath){
  path <- path.expand(filepath)
  chk <- tryCatch(
    {
      system2(path, "help",
              stdout = TRUE, stderr = FALSE)[1]
    },
    error = function(err) {
      return(NA)
    }
  )
  if(!is.na(chk)){
    ver <- strsplit(gsub(" |diamond|v", "", chk),"(", fixed = T)[[1]][1]
    vern <- strsplit(ver, ".", fixed = T)[[1]]
    vern <- as.numeric(sprintf("%s.%s%s", vern[1], vern[2], vern[3]))
    chk <- vern
  }
  return(chk)
}

#' @title ggplot2 theme for genespace
#' @description
#' \code{theme_genespace} specifies publication-style themes. Col here is the
#' color of the panel.background.
#' @rdname utils
#' @export
theme_genespace <- function(col = "black"){
  theme(panel.background = element_rect(fill = col),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(
          color = rgb(1, 1, 1, .2), size = .2, linetype = 2),
        panel.spacing = unit(.1, "mm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(
          angle = 90, family = "Helvetica", size = 5),
        strip.text.y.left = element_text(
          angle = 0, family = "Helvetica", size = 5),
        axis.title = element_text(family = "Helvetica", size = 6),
        plot.title = element_text(family = "Helvetica", size = 7))
}

#' @title download example data
#' @description
#' \code{download_exampleData} downloads chicken and human annotations from NCBI
#' @rdname utils
#' @importFrom utils download.file
#' @export
download_exampleData <- function(filepath){

  path <- path.expand(check_character(filepath))
  if(!dir.exists(path)){
    if(!dir.exists(dirname(path))){
      stop(sprintf("%s and its parent directory do not exist. provide a path to a valid directory."))
    }else{
      dir.create(path)
    }
  }

  hpath <- file.path(path, "human")
  if(!dir.exists(hpath))
    dir.create(hpath)
  cat(sprintf("Downloading human data to %s ... ", hpath))
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz",
    destfile = file.path(hpath, "GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz"))
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
    destfile = file.path(hpath, "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"))


  cpath <- file.path(path, "chicken")
  if(!dir.exists(cpath))
    dir.create(cpath)
  cat(sprintf("Done!\nDownloading chicken data to %s ...", cpath))
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz",
    destfile = file.path(cpath, "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz"))
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz",
    destfile = file.path(cpath, "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz"))
  cat("Done!")
  return(path)
}

#' @title linear interploation of missing positions
#' @description
#' \code{interp_approx} use approx to interpolate missing positions based on
#' the positions of the bounder x/y coordinates.
#' @rdname utils
#' @importFrom utils download.file
#' @importFrom stats approx
#' @export
interp_approx <- function(x, y){
  y <- as.numeric(y)
  if(all(is.na(y)))
    stop("must have some non-NA values in y")
  x <- as.numeric(x)
  if(any(is.na(x)))
    stop("x cannot contain NAs")
  if(any(is.na(y))){
    newy <- approx(
      x = x,
      y = y,
      ties = mean,
      xout = x[is.na(y)])$y
    y[is.na(y)] <- newy
  }
  return(y)
}

#' @title read syntenic hits for a genome
#' @description
#' \code{read_refGenomeSynHits} read in all syntenic hits files involving a
#' single reference genome and, where necessary, invert the hits so that the
#' reference genome is always the query (genome1).
#' @rdname utils
#' @importFrom utils download.file
#' @export
read_refGenomeSynHits <- function(gsParam,
                                  refGenome){

  query <- target <- ofID1 <- ofID2 <- NULL

  if(!"synteny" %in% names(gsParam))
    gsParam <- set_syntenyParams(gsParam)

  blMd <- data.table(gsParam$synteny$blast)
  nCores <- gsParam$params$nCores

  tnames <- c(
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "bitScore", "sameOG", "noAnchor", "isAnchor", "regID", "blkID")

  isq <- subset(blMd, query == refGenome)
  ist <- subset(blMd, target == refGenome & query != target)

  synh <- rbindlist(mclapply(isq$synHits, mc.cores = nCores, read_synHits))

  if(nrow(ist) > 0){
    synt <- rbindlist(mclapply(ist$synHits, mc.cores = nCores, read_synHits))
    synt <- synt[,tnames, with = F]
    setnames(synt, colnames(synh))
    synh <- rbind(synh, synt)
  }
  synh <- subset(synh, !duplicated(paste(ofID1, ofID2)))

  return(synh)
}

#' @title read all hits for a genome
#' @description
#' \code{read_refGenomeAllBlast} read in all hits files involving a
#' single reference genome and, where necessary, invert the hits so that the
#' reference genome is always the query (genome1).
#' @rdname utils
#' @importFrom utils download.file
#' @export
read_refGenomeAllBlast <- function(gsParam,
                                  refGenome){

  query <- target <- ofID1 <- ofID2 <- NULL

  if(!"synteny" %in% names(gsParam))
    gsParam <- set_syntenyParams(gsParam)

  blMd <- data.table(gsParam$synteny$blast)
  nCores <- gsParam$params$nCores

  tnames <- c(
    "ofID2", "chr2", "start2", "end2", "id2", "ord2", "genome2", "isArrayRep2",
    "ofID1", "chr1", "start1", "end1", "id1", "ord1", "genome1", "isArrayRep1",
    "pid", "length", "mismatches", "gapopenings", "queryStart", "queryEnd",
    "subjectStart", "subjectEnd", "Evalue", "bitScore", "sameOG", "noAnchor",
    "isAnchor", "isSyntenic", "regID", "blkID")

  isq <- subset(blMd, query == refGenome)
  ist <- subset(blMd, target == refGenome & query != target)

  synh <- rbindlist(mclapply(isq$allBlast, mc.cores = nCores, read_allBlast))

  if(nrow(ist) > 0){
    synt <- rbindlist(mclapply(ist$allBlast, mc.cores = nCores, read_allBlast))
    synt <- synt[,tnames, with = F]
    setnames(synt, colnames(synh))
    synh <- rbind(synh, synt)
  }
  synh <- subset(synh, !duplicated(paste(ofID1, ofID2)))

  return(synh)
}

#' @title subset the bed object within blocks
#' @description
#' \code{get_bedInBlk} splits the bed file so that two entries (query and target
#' ) match the physical bounds of blocks in the hits object.
#' @rdname utils
#' @importFrom utils download.file
#' @export
get_bedInBlk <- function(hits, bed){

  isAnchor <- blkID <- ord1 <- ord2 <- isArrayRep <- genome1 <- chr1 <-
    start1 <- end1 <- i.start1 <- i.end1 <- genome2 <- chr2 <- start2 <-
    end2 <- i.start2 <- i.end2 <- ofID1 <- ofID2 <- NULL

  anch <- subset(hits, isAnchor & !is.na(blkID))
  anchu1 <- with(anch, paste(ofID1, blkID))
  anchu2 <- with(anch, paste(ofID2, blkID))
  bc <- anch[,list(
    start1 = min(ord1), start2 = min(ord2),
    end1 = max(ord1), end2 = max(ord2)),
    by = c("genome1","chr1", "genome2", "chr2", "blkID")]

  b1 <- with(subset(bed, isArrayRep), data.table(
    genome1 = genome, chr1 = chr, start1 = ord,
    end1 = ord, ofID1 = ofID, ord1 = ord, og = og))

  b2 <- with(subset(bed, isArrayRep), data.table(
    genome2 = genome, chr2 = chr, start2 = ord,
    end2 = ord, ofID2 = ofID, ord2 = ord, og = og))

  setkey(bc, genome1, chr1, start1, end1)
  setkey(b1, genome1, chr1, start1, end1)
  fo1 <- foverlaps(bc, b1)
  fo1[,`:=`(start1 = i.start1, end1 = i.end1, i.start1 = NULL, i.end1 = NULL,
            start2 = NULL, end2 = NULL, chr2 = NULL, genome2 = NULL)]

  setkey(bc, genome2, chr2, start2, end2)
  setkey(b2, genome2, chr2, start2, end2)
  fo2 <- foverlaps(bc, b2)
  fo2[,`:=`(start2 = i.start2, end2 = i.end2, i.start2 = NULL, i.end2 = NULL,
            start1 = NULL, end1 = NULL, chr1 = NULL, genome1 = NULL)]
  fo1 <- subset(fo1, !paste(ofID1, blkID) %in% anchu1)
  fo2 <- subset(fo2, !paste(ofID2, blkID) %in% anchu2)
  out <- rbind(fo1, fo2, fill = T)
  out <- out[,c("ofID1", "ofID2", "ord1", "ord2", "blkID"), with = F]
  anchout <- subset(hits, isAnchor & !is.na(blkID))
  anchout <- anchout[,c("ofID1", "ofID2", "ord1", "ord2", "blkID"), with = F]

  out <- rbind(out, anchout)
  # out <- subset(out, !paste(ofID1, blkID) %in% hasAnch)
  # out <- subset(out, !duplicated(paste(ofID2, blkID)))
  setorder(out, blkID, ord1, ord2, na.last = T)

  out <- merge(bc, out, by = "blkID")
  out[,`:=`(start1 = NULL, start2 = NULL, end1 = NULL, end2 = NULL)]
  return(out)
}


#' @title calculate collinear arrays
#' @description
#' \code{add_array2bed} add array ID to the bed file
#' @rdname utils
#' @import data.table
#' @export
add_array2bed <- function(bed, synBuff, maxIter = 10, reorder = TRUE){

  # -- we can make arrays for everything except that we want to exclude the
  # arrays that are huge and problematic
  genome <- chr <- id <- arrayID <- nOGPlaces <- tord <- ofID <- NULL
  tmp <- data.table(bed)
  tmp[,tord := as.numeric(ord)]

  # -- set up the iteration
  cnt <- 1
  diffn <- 1
  tmp[,arrayID := sprintf(
    "tmp%s", as.integer(as.factor(paste(genome, chr, id))))]
  while(cnt <= maxIter && diffn > 0){

    # -- for each iteration, calculate clusters by the size of jump between
    # genes larger than the synBuffer
    cnt <- cnt + 1
    initn <- uniqueN(tmp$arrayID)
    genome <- chr <- og <- ord <- jumpLeft <- clus <- n <- NULL
    setkey(tmp, genome, chr, og, tord)
    if(reorder)
      tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
    tmp[,jumpLeft := c(synBuff + 1, diff(tord)),
        by = c("genome", "chr", "og")]
    tmp[,clus := as.integer(jumpLeft > synBuff), by = c("genome", "chr")]
    tmp[,clus := cumsum(clus), by = c("genome", "chr", "og")]
    tmp[,`:=`(
      arrayID = sprintf("tmp%s",
                        as.integer(as.factor(paste(genome, chr, og, clus)))),
      jumpLeft = NULL, clus = NULL)]

    # -- get the new order of genes based on array ID
    tmp[,n := .N, by = "arrayID"]
    tmp1 <- subset(tmp, n == 1)
    tmp2 <- subset(tmp, n > 1)
    tmp2[,tord := as.numeric(tord)]
    tmp2[,tord := mean(as.numeric(tord)), by = "arrayID"]
    tmp <- rbind(tmp1, tmp2)
    tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
    newn <- uniqueN(tmp$arrayID)
    diffn <- initn - newn
  }

  # -- relabel
  ogID <- NULL
  lab <- gsub(" ", "0",
              align_charRight(
                as.numeric(factor(tmp$arrayID,
                                  levels = unique(tmp$arrayID)))))
  arrv <- sprintf("Arr%s", lab); names(arrv) <- tmp$ofID
  bed[,arrayID := arrv[ofID]]

  tmp <- subset(bed, is.na(arrayID))
  wh <- with(tmp,
             as.numeric(as.factor(paste(genome, chr, og))))
  lab <- gsub(" ", "0", align_charRight(wh))
  wh <- which(is.na(bed$arrayID))
  bed$arrayID[wh] <- sprintf("NoArr%s", lab)
  return(bed)
}

#' @title calculate the array representatives
#' @description
#' \code{add_arrayReps2bed} add array representative genes to the combined bed
#' object
#' @rdname utils
#' @import data.table
#' @export
add_arrayReps2bed <- function(bed){
  n <- ord <- medOrd <- medDiff <- pepLen <- arrayID <- isRep <-
    isArrayRep <- ofID <- NULL
  bed[,n := .N, by = "arrayID"]
  tmp <- subset(bed, n > 1)
  tmp[,medOrd := as.numeric(median(ord)), by = "arrayID"]
  tmp[,medDiff := as.numeric(abs(medOrd - ord))]
  setorder(tmp, arrayID, medDiff, -pepLen)
  tmp$noAnchor[duplicated(tmp$arrayID)] <- TRUE
  tmp[,isRep := !duplicated(arrayID)]
  bed[,isArrayRep := ofID %in% tmp$ofID[tmp$isRep] | n == 1]
  tmp <- subset(bed, isArrayRep)
  tmp[,ord := frank(ord, ties.method = "dense"), by = "genome"]
  di <- tmp$ord; names(di) <- tmp$ofID
  return(bed)
}

#' @title write pan-gene file
#' @description
#' \code{write_pangenes} utility to correctly write in long-formatted pan=gene
#' text files.
#' @rdname utils
#' @export
write_pangenes <- function(x, filepath){
  hnames <- c(
    "ofID", "pgID", "interpChr", "interpOrd", "pgRepID", "genome", "og",
    "flag", "id", "chr", "start", "end", "ord")
  if(!all(hnames %in% colnames(x)))
    stop("pangenes data.table is malformed\n")
  x <- x[,hnames, with = F]
  fwrite(x, file = filepath, showProgress = F, quote = F, sep = "\t")
}

#' @title read pan-gene text file
#' @description
#' \code{read_pangenes} utility to correctly read in long-formatted pangene text
#' files.
#' @rdname utils
#' @export
read_pangenes <- function(x, filepath, ...){
  hnames <- c(
    "ofID", "pgID", "interpChr", "interpOrd", "pgRepID", "genome", "og",
    "flag", "id", "chr", "start", "end", "ord")
  cl <- c("numeric", "character", "logical")
  chk <- strsplit(readLines(filepath, 1), "\t")[[1]]
  if(!identical(chk, hnames))
    stop("pangenes.txt file is malformed\n")
  hc <- cl[c(2,1,2,1,2,2,1,2,2,2,1,1,1)]
  pangenes <- fread(
    filepath, na.strings = c("", "NA"), select = hnames,
    colClasses = hc, showProgress = F, ...)
  return(pangenes)
}

#' @title find nearest neighbor hits
#' @description
#' \code{find_nnHit} find anchor hits that are nearest to non-anchor xy position
#' @rdname utils
#' @import data.table
#' @export
find_nnHit <- function(x, y, isAnchor, radius){
  ##############################################################################
  # -- parameter checking
  if(length(x) != length(y))
    stop("x and y must be vectors of the same length\n")
  if(length(x) < 2)
    stop("x and y must have lengths > 1\n")
  x <- as.numeric(x)
  if(any(is.na(x)))
    stop("x values must all be numeric or integer\n")
  y <- as.numeric(y)
  if(any(is.na(y)))
    stop("y values must all be numeric or integer\n")
  if(sum(isAnchor) < 1)
    stop("at least one TRUE value must be given for isAnchor (logical vector of observations are anchors)\n")
  isAnchor <- as.logical(isAnchor)
  if(any(is.na(isAnchor)))
    stop("values given to isAnchor (logical vector of isAnchor observations are anchors) must all be coercible to logical\n")
  radius <- as.numeric(radius[1])
  if(is.na(radius) || is.null(radius))
    stop("radius must be a single numeric value > 0\n")
  radius <- as.numeric(radius[1])
  if(radius <= 0)
    stop("radius must be a single numeric value > 0\n")
  ##############################################################################
  # -- function to do it
  # -- get fixed radius nearest neighbors
  xy <- data.frame(x, y, isAnchor)
  fn <- frNN(
    x = subset(xy, isAnchor)[c("x", "y")],
    query = subset(xy, !isAnchor)[c("x", "y")],
    eps = radius)
  # -- get index of nearest neighbor hits
  wh <- sapply(fn$id, function(x) x[1])
  return(wh)
}

#' @title flag proximate hits
#' @description
#' \code{flag_hitsInRadius} given a vector of anchors, pulls xy positions within
#' radius of anchors using dbscan
#' @rdname utils
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
flag_hitsInRadius <- function(x, y, isAnchor, radius){
  ##############################################################################
  # -- parameter checking
  if(length(x) != length(y))
    stop("x and y must be vectors of the same length\n")
  if(sum(isAnchor) == 0)
    return(rep(FALSE, length(x)))
  if(sum(!isAnchor) == 0)
    return(rep(TRUE, length(x)))
  if(length(x) < 2)
    stop("x and y must have lengths > 1\n")
  x <- as.numeric(x)
  if(any(is.na(x)))
    stop("x values must all be numeric or integer\n")
  y <- as.numeric(y)
  if(any(is.na(y)))
    stop("y values must all be numeric or integer\n")

  isAnchor <- as.logical(isAnchor)
  if(any(is.na(isAnchor)))
    stop("values given to isAnchor (logical vector of isAnchor observations are anchors) must all be coercible to logical\n")

  radius <- as.numeric(radius[1])
  if(is.na(radius) || is.null(radius))
    stop("radius must be a single numeric value > 0\n")
  radius <- as.numeric(radius[1])
  if(radius <= 0)
    stop("radius must be a single numeric value > 0\n")

  ##############################################################################
  # -- function to do it
  # -- get fixed radius nearest neighbors
  if(all(isAnchor) || all(!isAnchor)){
    return(isAnchor)
  }else{
    xy <- data.frame(x = x, y = y, isAnchor = isAnchor, inBuffer = isAnchor)
    fn <- frNN(
      x = subset(xy, isAnchor)[c("x", "y")],
      query = subset(xy, !isAnchor)[c("x", "y")],
      eps = radius)

    # -- subset to positions with any nearest neighbors in radius
    hasAnch <- fn$id[sapply(fn$id, length) > 0]
    wh <- as.integer(names(hasAnch))
    # -- set these observations to true and return
    if(length(wh) > 0)
      xy$inBuffer[wh] <- TRUE
    return(xy$inBuffer)
  }
}

#' @title flag hits within blocks
#' @description
#' \code{flag_hitsInBlk} finds hits within the bounding coordinates of blocks
#' @rdname utils
#' @import data.table
#' @export
flag_hitsInBlk <- function(x, y, blkID){

  isAnchor <- start1 <- start2 <- end1 <- end2 <-
    i.blkID <- blkID1 <- blkID2 <- id <- NULL

  ##############################################################################
  # -- parameter checking
  if(length(x) != length(y))
    stop("x and y must be vectors of the same length\n")
  if(length(x) < 2)
    stop("x and y must have lengths > 1\n")
  x <- as.numeric(x)
  if(any(is.na(x)))
    stop("x values must all be numeric or integer\n")
  y <- as.numeric(y)
  if(any(is.na(y)))
    stop("y values must all be numeric or integer\n")

  if(length(blkID) != length(x))
    stop("x, y and blkID must have the same length\n")

  xy <- data.table(
    start1 = x, end1 = x, start2 = y, end2 = y,
    isAnchor = !is.na(blkID), id = 1:length(x), blkID = blkID)

  if(all(is.na(blkID)) || all(!is.na(blkID))){
    return(xy$isAnchor)
  }else{
    blks <- subset(xy, isAnchor)[,list(
      start1 = min(start1), end1 = max(start1),
      start2 = min(start2), end2 = max(start2)),
      by = "blkID"]

    b1 <- blks[,c("start1" ,"end1", "blkID")]
    b2 <- blks[,c("start2" ,"end2", "blkID")]
    setkey(b1, start1, end1)
    setkey(b2, start2, end2)
    setkey(xy, start1, end1)
    t1 <- foverlaps(b1, xy)
    t1[,`:=`(i.start1 = NULL, i.end1 = NULL, blkID1 = i.blkID, i.blkID = NULL)]
    setkey(t1, start2, end2)
    t2 <- foverlaps(b2, t1)
    t2[,`:=`(i.start2 = NULL, i.end2 = NULL, blkID2 = i.blkID, i.blkID = NULL)]
    out <- subset(t2, blkID1 == blkID2)
    bv <- out$blkID1; names(bv) <- with(out, paste(start1, start2))
    xy[,blkID := bv[paste(start1, start2)]]
    setkey(xy, id)
    return(xy$blkID)
  }
}

#' @title get ordered tip labels from a phylogenetic tree
#' @description
#' \code{get_orderedTips} Respect ordering of tree when ladderized
#' @rdname utils
#' @export
get_orderedTips <- function(treFile, ladderize = TRUE, genomeIDs){
  if(requireNamespace("ape", quietly = T)){
    tre <- ape::read.tree(treFile)
    if(all(genomeIDs %in% tre$tip.label)){
      if(ladderize)
        tre <- ape::ladderize(tre)
      is_tip <- tre$edge[,2] <= length(tre$tip.label)
      ordered_tips <- tre$edge[is_tip, 2]
      outlabs <- tre$tip.label[ordered_tips]
      outlabs <- outlabs[outlabs %in% genomeIDs]
      return(outlabs)
    }else{
      cat(" some genomeIDs are not in the tree file ")
      return(genomeIDs)
    }
  }else{
    cat("the ape library is not available, not using the tree to order genomes ")
    return(genomeIDs)
  }
}

#' @title get all pairwise combination of hits between two genomes
#' @description
#' \code{pull_pairwise} Builds new pairwise files in /syntenicHits
#' @rdname utils
#' @import data.table
#' @export
pull_pairwise <- function(gsParam, verbose){

  sameOG <- genome <- NULL

  bl <- data.table(gsParam$synteny$blast)
  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))

  for(i in 1:nrow(bl)){
    qu <- bl$query[i]
    ta <- bl$target[i]
    blf <- bl$synHits[i]
    blo <- file.path(gsParam$paths$syntenicHits,
                     sprintf("%s_vs_%s.pairwise.txt.gz", qu, ta))
    h <- subset(read_synHits(blf), sameOG)[,c("ofID1", "ofID2", "regID", "blkID")]
    a1 <- with(subset(bed, genome == qu), data.table(
      ofID1 = ofID, arrayID1 = arrayID))
    a2 <- with(subset(bed, genome == ta), data.table(
      ofID2 = ofID, arrayID2 = arrayID))
    ho <- merge(a1, merge(a2, h, by  = "ofID2", allow.cartesian = T, all.y = T),
                by = "ofID1", allow.cartesian = T, all.y = T)
    b1 <- with(subset(bed, genome == qu), data.table(
      genome1 = genome, id1 = id, arrayID1 = arrayID))
    b2 <- with(subset(bed, genome == ta), data.table(
      genome2 = genome, id2 = id, arrayID2 = arrayID))
    ho[,`:=`(ofID1 = NULL, ofID2 = NULL)]
    ho2 <- merge(b2, ho, by  = "arrayID2", allow.cartesian = T, all.y = T)
    ho2 <- subset(ho2, !duplicated(ho2))
    ho2 <- merge(b1, ho2, by = "arrayID1", allow.cartesian = T, all.y = T)
    ho2 <- subset(ho2, !duplicated(ho2))
    ho2[,`:=`(arrayID1 = NULL, arrayID2 = NULL)]
    if(verbose)
      cat(sprintf("\t%s vs. %s: synHits = %s, pairwiseHits = %s\n", qu, ta, nrow(h), nrow(ho2)))
    fwrite(ho2, file = blo, sep = "\t", quote = F, verbose = F, showProgress = F)
  }
}
