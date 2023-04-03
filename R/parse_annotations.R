#' @title Accessory function to help build GENESPACE input files
#' @description
#' \code{parse_annotations} Peptide and gff3 gene annotation matching and
#' conversion to .bed format. Speeds up downstream compute and catches problems
#' with annotation files. This is NOT required for GENESPACE, but does help
#' get the input files in order. There are many other methods to convert
#' gff3 --> bed and match the names with fasta headers.
#'
#' @name parse_annotations
#' @param rawGenomeRepo file path to the location of gff3 and fasta annotations
#' @param genomeDirs character vector giving exact matches to subdirectories
#' in rawGenome repo.
#' @param genomeIDs character vector of length equal to genomeDirs. By default,
#' takes values from genomeDirs, but, if specified, re-names the files
#' accordingly. Useful if you want to shorten the names of genomes in your
#' GENESPACE run.
#' @param gffString regular expression of length 1 specifying the string to
#' search for in rawGenomeRepo/genomeDirs that will exactly match the gff3-
#' formatted annotation. Default is any text or .gz file ending in gff3, gff.
#' @param faString same as gffString but for the fasta-formatted peptide
#' annotation. Default is any text or .gz file ending in fa, fasta or faa.
#' @param genespaceWd file.path of length 1 specifying the GENESPACE working
#' directory. Will make two subdirectories: /bed and /peptide for the parsed
#' annotations
#' @param presets character string: "none", "phytozome" or "ncbi" which sets
#' the below parameters to parse phytozome or ncbi-formatted annotations
#' correctly. See details.
#' @param gffIdColumn character, specifying the field name in the gff3
#' attributes column.
#' @param minPepLen numeric, specifying the shortest peptide (in daltons) to be
#' kept
#' @param dropDuplicates logical, should only one of a set of duplicated
#' peptide sequences be kept?
#' @param removeNonAAs logial, should "." and "-" characters be stripped from
#' the amino acids?
#' @param headerEntryIndex integer specifying the field index in the fasta
#' header which contains the gene ID information to match with the gff.
#' @param headerSep character used as a field delimiter in the fasta header.
#' @param gffStripText regular expression of length 1 specifying a gsub command
#' to remove text from the gff ID.
#' @param headerStripText like gffStripText, but for the fasta header
#' @param chrIdDictionary a named vector where the names are the values in the
#' first ("seqnames") gff3 column and the element names in the vector are
#' the values to replace.
#' @param troubleShoot logical, should the raw and parsed files be printed?
#' @param overwrite logical, should existing files be overwritten?
#' @param path2fasta deprecated, kept to maintain backwards compatibility
#' @param path2gff deprecated, kept to maintain backwards compatibility
#' @param genomeID single genomeID to consider
#' @param convertSpecialCharacters Character string with a non-special character
#' of length 1. Replaces special characters (punctionation other
#' than ".", "-", and "_") if they are present in the gene IDs.
#' @param ... additional arguments passed on
#'
#' @details parse_annotations assumes that you have a 'rawGenomeRepo' directory
#' that contains a subdirectory for each genome to parse. These subdirectory
#' names are given in "genomeDirs". So, if rawGenomeRepo = "~/Destop/genomeRepo"
#' and genomeDirs = c("human", "mouse"), then parse_annotations assumes there
#' are two directories: ~/Desktop/genomeRepo/human and ../mouse. Each of these
#' dicectories must contain a gff3-formatted gene annotation and a fasta-
#' formatted peptide annotation. These annotation files can be further nested
#' in the subdirectories, but each must be named with a uniquely findable
#' "gffString" and "faString". If multiple (or no) files match these strings,
#' an error will be returned.
#'
#' Given differences in how gff3 and fasta headers are constructed, there are a
#' number of parameters to choose how the files should be matched. Unless the
#' files come from phytozome or NCBI, these need to be chosen manually. If you
#' have differently named or formatted files, you can run parse_annotations
#' several times with different parameters and paths to the files.
#'
#' For each genomeDir, the pair of files are read in, parsed, matched, then
#' written to $genespaceWd/$genomeID/peptide and $genespaceWd/$genomeID/bed
#' respectively. By default, the genomeID is the same as the genomeDir, but
#' this can be customized.
#'
#' Presets: In the case of "ncbi", this assumes that you have downloaded the
#' 'translated_cds' peptide file and gene.gff3 annotation. Given this file,
#' it uses the following parameters: gffIdColumn <- "gene"; headerEntryIndex <- 2;
#' headerSep <- " "; gffStripText <- ""; headerStripText <- "gene=|\\[|\\]"
#' Present 'ncbi"  also builds a chrIdDictionary to re-name sequenceIDs based on
#' entries labeled "chromosome" in the third gff3 column.
#'
#' Phytozome presets: gffIdColumn <- "Name"; headerEntryIndex <- 4;
#' headerSep <- " "; gffStripText <- ""; headerStripText <- "locus="
#'
#' @return a data.table containing the file paths to the raw and parsed
#' annotations.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#'
#' @title match gff and fasta transcript annotations
#' @description
#' \code{parse_annotations} parse gff into a bed format with one entry per
#' primary transcript, and ensure that the peptide fasta headers match the
#' name column
#' @rdname parse_annotations
#' @import data.table
#' @export
parse_annotations <- function(rawGenomeRepo,
                              genomeDirs,
                              genomeIDs = genomeDirs,
                              gffString = "gff$|gff3$|gff3\\.gz$|gff\\.gz",
                              faString = "fa$|fasta$|faa$|fa\\.gz$|fasta\\.gz|faa\\.gz",
                              genespaceWd,
                              minPepLen = 0,
                              dropDuplicates = FALSE,
                              removeNonAAs = FALSE,
                              presets = "none",
                              gffIdColumn = "Name",
                              headerEntryIndex = 4,
                              headerSep = " ",
                              gffStripText = "",
                              headerStripText = "locus=",
                              convertSpecialCharacters = "_",
                              chrIdDictionary = NULL,
                              troubleShoot = FALSE,
                              overwrite = FALSE){

  if(!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  # -- make sure you can write the output directory
  genespaceWd <- check_character(genespaceWd)
  if(!dir.exists(genespaceWd)){
    if(!dir.exists(dirname(genespaceWd))){
      stop("could not find genespace directory nor its parent directory")
    }else{
      dir.create(genespaceWd)
    }
  }

  # -- make sure the repo directory is good
  rawGenomeRepo <- check_character(rawGenomeRepo)
  if(!dir.exists(rawGenomeRepo))
    stop("could not find rawGenomeRepo directory")
  rawDis <- file.path(rawGenomeRepo, genomeDirs)
  if(any(!dir.exists(rawDis)))
    stop("could not find the following raw genome directories: \n\t",
         paste(basename(rawDis[!dir.exists(rawDis)]), collapse = "\n\t"))

  # -- get parameters in order
  genomeDirs <- check_character(genomeDirs, na.rm = T)
  genomeIDs <- check_character(genomeIDs, na.rm = T)
  gffStripText <- check_character(gffStripText, na.rm = T, onlySingleValue = T)
  headerStripText <- check_character(headerStripText, na.rm = T, onlySingleValue = T)
  gffIdColumn <- check_character(gffIdColumn, na.rm = T, onlySingleValue = T)
  faString <- check_character(
    faString, na.rm = T, default = "gff", onlySingleValue = T)
  gffString <- check_character(
    gffString, na.rm = T, default = "gff", onlySingleValue = T)
  presets <- check_character(
    presets, na.rm = T, default = "none", onlySingleValue = T)
  headerSep <- check_character(
    headerSep, na.rm = T, default = "none", onlySingleValue = T)
  headerEntryIndex <- check_integer(
    headerEntryIndex, min = 1, na.rm = T, onlySingleValue = T)

  if(length(genomeDirs) != length(genomeIDs))
    stop("genomeDirs (subdirectory in rawGenomeRepo) and genomeIDs must be the same length")

  outPaths <- rbindlist(lapply(1:length(genomeIDs), function(i){
    # -- get input file locations
    dirID <- genomeDirs[i]
    genID <- genomeIDs[i]
    gf <- list.files(
      file.path(rawGenomeRepo, dirID),
      pattern = gffString, recursive = T, full.names = T)
    fa <- list.files(
      file.path(rawGenomeRepo, dirID),
      pattern = faString, recursive = T, full.names = T)

    if(length(fa) == 0)
      stop(sprintf("could not find fasta file with faString = %s in %s. Files available:\n\t%s",
           faString, file.path(rawGenomeRepo, dirID),
           paste(list.files(file.path(rawGenomeRepo, dirID), recursive = T), collapse = "\n\t")))
    if(length(gf) == 0)
      stop(sprintf("could not find gff file with gffString = %s in %s. Files available:\n\t%s",
           gffString, file.path(rawGenomeRepo, dirID),
           paste(list.files(file.path(rawGenomeRepo, dirID), recursive = T), collapse = "\n\t")))
    if(length(fa) > 1)
      stop(sprintf("found multiple fasta files with faString = %s:\n\t%s",
                   faString, file.path(rawGenomeRepo, dirID),
                   paste(gf, collapse = "\n\t")))
    if(length(gf) > 1)
      stop(sprintf("found multiple gff files with gffString = %s:\n\t%s",
                   gffString, file.path(rawGenomeRepo, dirID),
                   paste(fa, collapse = "\n\t")))
    # do the matching
    mtch <- match_fasta2gff(
      path2fasta = fa, path2gff = gf, genespaceWd = genespaceWd,
      genomeID = genID, troubleShoot = troubleShoot,
      headerEntryIndex = headerEntryIndex, presets = presets,
      gffIdColumn = gffIdColumn, headerSep = headerSep,
      gffStripText = gffStripText, headerStripText = headerStripText,
      chrIdDictionary = chrIdDictionary,
      convertSpecialCharacters = convertSpecialCharacters)

    return(data.table(
      gffFileIn = gf, faFileIn = fa,
      bedFileOut = mtch$bed, faFileOut = mtch$fasta))
  }))
  return(outPaths)
}


#' @title match a single peptide fasta to a gff
#' @description
#' \code{match_fasta2gff} engine for reading, parsing and writing annotation
#' files
#' @rdname parse_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width readAAStringSet AAStringSet DNA_ALPHABET vcountPattern replaceAt vmatchPattern
#' @importFrom utils head
#' @export
match_fasta2gff <- function(path2fasta,
                            path2gff,
                            genespaceWd,
                            genomeID,
                            presets,
                            gffIdColumn,
                            headerEntryIndex,
                            headerSep,
                            minPepLen,
                            dropDuplicates,
                            removeNonAAs,
                            gffStripText,
                            headerStripText,
                            chrIdDictionary,
                            convertSpecialCharacters,
                            troubleShoot){
  if(!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")
  # -- ensure the place to write to is valid
  wd <- path.expand(genespaceWd)
  path2fasta <- path.expand(path2fasta)
  path2gff <- path.expand(path2gff)
  if(!dir.exists(wd))
    stop("wd must be an existing directory\n")
  # -- ensure input files are valid
  if(!file.exists(path2fasta))
    stop("specified path2fasta does not exist\n")
  if(!file.exists(path2gff))
    stop("specified path2gff does not exist\n")
  # -- make sure presets are ok
  presets <- match.arg(
    presets,
    choices = c("phytozome", "ncbi", "none"))
  # -- make sure the chrIdDictionary is ok
  if(!is.null(chrIdDictionary)){
    if(!is.vector(chrIdDictionary))
      stop("if specified, chrIdDictionary must be a vector")
    if(is.null(names(chrIdDictionary)))
      stop("if specified, chrIdDictionary must be a NAMED vector")
    if(any(duplicated(names(chrIdDictionary))))
      stop("if specified, chrIdDictionary must be a NAMED vector with unique names")
  }
  # -- set phytozome presets
  if(presets == "phytozome"){
    gffIdColumn <- "Name"
    headerEntryIndex <- 4
    headerSep <- " "
    gffStripText <- ""
    headerStripText <- "locus="
  }else{
    # -- set ncbi presets
    if(presets == "ncbi"){
      gffIdColumn <- "gene"
      headerEntryIndex <- 2
      headerSep <- " "
      gffStripText <- ""
      headerStripText <- "gene=|\\[|\\]"
      # -- get ncbi chromosome - seqid dictionary from the gff
      if(is.null(chrIdDictionary)){
        chrIDs2mask <- c("NA", "none", "unknown")
        chrIDs <- data.table(data.frame(rtracklayer::readGFF(
          filepath = path2gff,
          filter = list(type = "region"),
          tags = c("chromosome"))))

        # -- ensure that problematic chromosome names are not replaced
        wh <- which(is.na(chrIDs$chromosome))
        chrIDs$chromosome[wh] <- chrIDs$seqid[wh]
        wh <- which(tolower(chrIDs$chromosome) %in% chrIDs2mask)
        chrIDs$chromosome[wh] <- chrIDs$seqid[wh]

        seqid <- chromosome <- NULL
        chrIDs <- subset(chrIDs, !is.na(seqid) & !is.na(chromosome))
        chrIdDictionary <- as.character(chrIDs$chromosome)
        names(chrIdDictionary) <- as.character(chrIDs$seqid)
        chrIdDictionary <- chrIdDictionary[!duplicated(names(chrIdDictionary))]
      }
    }else{
      # -- ensure custom paramters are OK
      headerEntryIndex <- check_integer(
        headerEntryIndex, min = 1, max = Inf, default = 4, onlySingleValue = T)
      gffIdColumn <- as.character(gffIdColumn[1])
      headerSep <- as.character(headerSep[1])
      gffStripText <- as.character(gffStripText[1])
      headerStripText <- as.character(headerStripText[1])
    }
  }

  # -- read in the peptide annotations
  fa <- Biostrings::readAAStringSet(path2fasta)
  if(troubleShoot){
    cat("\n### first 6 fasta headers before parsing ... \n")
    cat(head(names(fa)), sep = "\n")
  }
  # -- ensure that they are peptides, not DNA
  tmp <- gsub("[^A-Za-z]", "", paste(fa[1:100], collapse = ""))
  dnaa <- paste(DNA_ALPHABET, collapse = "|")
  tmp <- gsub(gsub("[^A-Za-z]", "", dnaa), "", tmp)
  if(nchar(tmp) == 0)
    stop("Only DNA characters found ... is this a peptide sequence?\n")

  # -- parse the peptide headers
  names(fa) <- gsub(
    headerStripText, "", sapply(names(fa), function(y)
      strsplit(y, headerSep)[[1]][headerEntryIndex]))

  # -- remove all special characters
  names(fa) <- gsub("[^a-zA-Z0-9_.-]", convertSpecialCharacters, names(fa))

  if(troubleShoot){
    cat("\n### first 6 fasta headers after parsing ... \n")
    cat(head(names(fa)), sep = "\n")
  }
  # -- ensure no duplicated headers
  if(any(duplicated(names(fa)))){
    fa <- fa[order(-width(fa)),]
    fa <- fa[!duplicated(names(fa)),]
  }

  nfa <- length(fa)
  # -- read in the gff
  end <- start <- width <- id <- seqid <- chr <- NULL
  if(troubleShoot){
    cat("\n### first 6 gff lines before parsing ... \n")
    if(grepl(".gz$", path2gff)){
      system(sprintf("cat %s | gunzip | grep -v '#' | head", path2gff))
    }else{
      system(sprintf("cat %s | grep -v '#' | head", path2gff))
    }
  }
  gff <- data.table(rtracklayer::readGFF(path2gff, tags = gffIdColumn))
  setnames(gff, gffIdColumn, "id")
  gff[,id := gsub(gffStripText, "", id)]
  if(troubleShoot){
    cat("\n### first 6 gff lines after parsing ... \n")
    print(head(gff))
  }

  # -- remove special characters
  gff[,id := gsub("[^a-zA-Z0-9_.-]", convertSpecialCharacters, id)]

  gff <- subset(gff, id %in% names(fa))
  gff[,seqid := as.character(seqid)]

  # -- check if any duplicates
  if(any(duplicated(gff$id))){
    gff[,width := end - start]
    setorder(gff, id, -width)
    gff <- subset(gff, !duplicated(id))
  }
  # -- match positions
  fa <- fa[gff$id]

  # -- rename chromosomes
  if(!is.null(chrIdDictionary)){
    if(any(duplicated(names(chrIdDictionary))))
      stop("all names in chrIdDictionary must be unique\n")
    gff[,chr := chrIdDictionary[seqid]]
    wh <- which(!gff$seqid %in% names(chrIdDictionary))
    gff$chr[wh] <- gff$seqid[wh]
  }else{
    gff[,chr := seqid]
  }

  # -- convert to 0-index 1-open bed
  gff <- gff[,c("chr", "start", "end", "id")]
  gff[,start := start - 1]

  if(troubleShoot){
    cat("\n### first 6 bed lines after full parsing (and potential chr re-name)\n")
    print(head(gff))
  }

  # -- check for .'s or -'s in the fa
  nDot <- sum(vcountPattern(pattern = ".", subject = fa))
  nDash <- sum(vcountPattern(pattern = "-", subject = fa))

  if((nDot + nDash) > 0){
    if(!removeNonAAs)
      stop("some of the peptides have '.' or '-' in the sequence. Orthofinder can't handle this. Either remove them manually or set removeNonAAs to TRUE\n")

    if(nDash > 0)
      suppressWarnings(fa <- replaceAt(
        x = fa,
        at = vmatchPattern("-", fa, fixed = TRUE),
        value = ""))

    if(nDot > 0)
      suppressWarnings(fa <- replaceAt(
        x = fa,
        at = vmatchPattern(".", fa, fixed = TRUE),
        value = ""))
  }

  # -- drop short aas
  if(minPepLen > 0)
    fa <- fa[width(fa) >= minPepLen]

  # -- drop duplicates
  if(dropDuplicates)
    fa <- fa[!duplicated(as.character(fa))]

  # -- final match of ids
  gff <- subset(gff, id %in% names(fa))
  setkey(gff, chr, start, end, id)
  fa <- fa[gff$id]

  # -- write output
  wdPep <- file.path(wd, "peptide")
  if(!dir.exists(wdPep))
    dir.create(wdPep)
  wdBed <- file.path(wd, "bed")
  if(!dir.exists(wdBed))
    dir.create(wdBed)
  cat(sprintf("%s: n unique sequences = %s, n matched to gff = %s\n",
              genomeID, nfa, nrow(gff)))

  gffout <- file.path(wdBed, sprintf("%s.bed", genomeID))
  faout <- file.path(wdPep, sprintf("%s.fa", genomeID))
  fwrite(gff, file = gffout, col.names = FALSE, quote = FALSE, sep = "\t")
  writeXStringSet(fa, filepath = faout)
  return(list(bed = gffout, fasta = faout))
}

#' @title Parse ncbi-formatted annotations
#' @description
#' \code{parse_ncbi} a shortcut for parse_annotations(preset = "ncbi") to
#' maintain backwards compatibility with < v1.0.0.
#' @rdname parse_annotations
#' @export
parse_ncbi <- function(rawGenomeRepo,
                       genomeDirs,
                       genomeIDs = genomeDirs,
                       gffString = "gff$|gff3$|gff3\\.gz$|gff\\.gz",
                       faString = "fa$|fasta$|faa$|fa\\.gz$|fasta\\.gz|faa\\.gz",
                       genespaceWd,
                       troubleShoot = FALSE,
                       ...){

  outPaths <- parse_annotations(
    rawGenomeRepo = rawGenomeRepo,
    genomeDirs = genomeDirs,
    genomeIDs = genomeIDs,
    gffString = gffString,
    faString = faString,
    genespaceWd = genespaceWd,
    presets = "ncbi",
    chrIdDictionary = NULL,
    troubleShoot = troubleShoot,
    ...)

  return(outPaths)
}

#' @title Parse phytozome-formatted annotations
#' @description
#' \code{parse_phytozome} a shortcut for parse_annotations(preset = "phytozome")
#' to maintain backwards compatibility with < v1.0.0.
#' @rdname parse_annotations
#' @export
parse_phytozome <- function(rawGenomeRepo,
                            genomeDirs,
                            genomeIDs = genomeDirs,
                            gffString = "gff$|gff3$|gff3\\.gz$|gff\\.gz",
                            faString = "fa$|fasta$|faa$|fa\\.gz$|fasta\\.gz|faa\\.gz",
                            genespaceWd,
                            troubleShoot = FALSE,
                            ...){

  outPaths <- parse_annotations(
    rawGenomeRepo = rawGenomeRepo,
    genomeDirs = genomeDirs,
    genomeIDs = genomeIDs,
    gffString = gffString,
    faString = faString,
    genespaceWd = genespaceWd,
    presets = "phytozome",
    chrIdDictionary = NULL,
    troubleShoot = troubleShoot,
    ...)

  return(outPaths)
}

#' @title Deprecated
#' @description
#' \code{parse_faHeader} function to maintain backwards compatibility
#' @rdname parse_annotations
#' @export
parse_faHeader <- function(...){
  cat(strwrap(
    "Thanks for switching to GENESPACE v1! To avoid confusion due to the new
  specification, parse_faHeader is now deprecated. See ?parse_annotations
  to format your annotations."),
    sep = "\n")
}
