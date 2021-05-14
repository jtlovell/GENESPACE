#' @title Functions to convert from raw to genespace annotation formats
#' @description
#' \code{format_annotations} Format gff and fasta for genespace. The two parse
#' functions (parse_gff & parse_faHeader) take a single file and return a
#' formatted R object. format_annotations is a wrapper for these two that reads
#' and writes the full genespace-formatted files.
#'
#' @name format_annotations
#' @param genomeIDs character vector of the same length as paths. Names for
#' each genome.
#' @param path2rawGff3 character vector coercible to file.paths. This points to
#' the raw gff3-formatted annotation (.gz or gff3)
#' @param path2rawPeptide character vector coercible to a file.path. This points
#' to the raw fasta-formatted primary peptide annotation.
#' @param path2gff character vector coercible to file.paths. This points to
#' the location to save the parsed gff3-like file (.gz or gff3)
#' @param path2peptide character vector coercible to a file.path. This points
#' to location to save the re-named fasta-formatted primary peptide annotation.
#' @param gffEntryType character, specifying which attribute type should be
#' retained. This is the 3rd column of a gff3-formatted annotation. Can be a
#' vector of length 1, if all files should be parsed identically, or a vector
#' with the same length as file paths and genomeIDs, specifying different
#' parsing parameters for each genome.
#' @param gffIdColumn character, specifying the field name in the gff3
#' attributes column.Can be a vector of length 1, if all files should be parsed
#' identically, or a vector with the same length as file paths and genomeIDs,
#' specifying different parsing parameters for each genome.
#' @param headerEntryIndex integer specifying the field index in the fasta
#' header which contains the gene ID information to match with the gff. Can be a
#' vector of length 1, if all files should be parsed identically, or a vector
#' with the same length as file paths and genomeIDs, specifying different
#' parsing parameters for each genome.
#' @param headerSep character used as a field delimiter in the fasta header.
#' Can be a
#' vector of length 1, if all files should be parsed identically, or a vector
#' with the same length as file paths and genomeIDs, specifying different
#' parsing parameters for each genome.
#' @param headerStripText character or regex to remove from the fasta header.
#' Can be a vector of length 1, if all files should be parsed identically,
#' or a vector with the same length as file paths and genomeIDs, specifying
#' different parsing parameters for each genome.
#' @param troubleshoot logical, should the raw and parsed files be printed?
#' @param path2rawFasta character string coercible to a file path, pointing to
#' the location of the unparsed fasta file.
#' @param overwrite logical, should existing files be overwritten? If FALSE,
#' and all files are present, just returns a named vector of gff/pep file
#' locations
#' @param verbose logical, should updates be printed to the console?
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run

#' @note \code{format_annotations} is a generic name for the functions documented.
#' \cr
#' If called, \code{format_annotations} returns its own arguments.
#'

#' @title parse phytozome-formatted annotations
#' @description
#' \code{parse_phyotozome} parse phytozome-formatted annotations
#' @rdname format_annotations
#' @import data.table
#' @export
parse_phyotozome <- function(gsParam, overwrite = F, genomeIDs = NULL){

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomeIDs

  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  pytz_fun <- function(gffIn, gffOut, pepIn, pepOut, verbose){
    if(verbose)
      cat("\tReading gff ... ")
    # read in the gene gff
    gff <- data.table(data.frame(rtracklayer::readGFF(
      filepath = gffIn,
      filter = list(type = "gene"),
      tags = "Name")))
    if(verbose)
      cat(sprintf("found %s protein coding genes\n", nrow(gff)))

    # read in the peptide fa
    if(verbose)
      cat("\tReading peptide fasta ... ")
    fa <- readAAStringSet(pepIn)

    # rename peptide headers
    tmp <- sapply(names(fa), function(x) strsplit(x, " ")[[1]][4])
    tmp2 <- substr(tmp, 7, nchar(tmp))
    names(fa) <- tmp2
    nfa <- length(fa)

    # drop duplicates, keeping the longest
    fa <- fa[order(-width(fa))]
    fa <- fa[!duplicated(names(fa))]
    if(verbose)
      cat(sprintf("found %s / %s total and unique entries\n\tMerging fa and gff",
                  nfa, nrow(gff)))

    # join the two
    gff <- subset(gff, Name %in% names(fa))
    fa <- fa[gff$Name]
    gff <- subset(gff, Name %in% names(fa))

    if(verbose)
      cat(sprintf(" found %s matching entires\n", nrow(gff)))

    gff[,seqid := factor(seqid, levels = unique(seqid))]
    setkey(gff, seqid, start, end)

    # write output
    gff <- with(gff, data.table(
      chr = seqid,
      start = start,
      end = end,
      id = Name,
      strand = strand,
      ord = 1:length(Name)))

    fwrite(gff, file = gffOut, sep = "\t", quote = F)
    writeXStringSet(fa, filepath = pepOut, compress = "gzip")
  }

  rawGff <- gsParam$gffRaw
  rawPep <- gsParam$peptideRaw
  outGff <- file.path(gsParam$gff, sprintf("%s.gff.gz", genomeIDs))
  outPep <- file.path(gsParam$peptide, sprintf("%s.fa.gz", genomeIDs))
  verbose <- gsParam$verbose
  names(outPep) <- names(outGff) <- genomeIDs

  if(!all(file.exists(outGff)) | !all(file.exists(outPep)) | overwrite){
    for(i in genomeIDs){
      if(verbose)
        cat(sprintf("Parsing annotations: %s\n",i))
      pytz_fun(
        gffIn = rawGff[i], gffOut = outGff[i],
        pepIn = rawPep[i], pepOut = outPep[i], verbose = verbose)
    }
    if(verbose)
      cat("Done!")
  }
  return(list(gff = outGff, peptide = outPep))
}


#' @title parse NCBI-formatted annotations
#' @description
#' \code{parse_ncbi} parse NCBI-formatted annotations
#' @rdname format_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width readAAStringSet
#' @export
parse_ncbi <- function(gsParam, overwrite = F, genomeIDs = NULL){

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomeIDs

  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  ncbi_fun <- function(gffIn, gffOut, pepIn, pepOut, verbose){
    if(verbose)
      cat("\tReading gff ... ")
    # read in the gene gff
    gff <- data.table(data.frame(rtracklayer::readGFF(
      filepath = gffIn,
      filter = list(type = "gene"),
      tags = c("gene", "gene_biotype"))))
    gff <- subset(gff, gene_biotype == "protein_coding")
    gff <- subset(gff, !duplicated(gene))
    gff[,`:=`(seqid = as.character(seqid), gene = as.character(gene))]
    gchr <- gff[,list(nGene = .N), by = "seqid"]
    if(verbose)
      cat(sprintf("found %s protein coding genes\n", nrow(gff)))

    # read in the region gff
    chrIDs <- data.table(data.frame(rtracklayer::readGFF(
      filepath = gffIn,
      filter = list(type = "region"),
      tags = c("chromosome"))))
    chrIDs[,`:=`(seqid = as.character(seqid), chromosome = as.character(chromosome))]

    uchrs <- c(unique(chrIDs$chromosome), unique(chrIDs$seqid))
    uchrs <- uchrs[!duplicated(uchrs)]

    # choose largest regions / chromosome as rep
    gid <- chrIDs[,list(nbp = end - start), by = c("chromosome","seqid")]
    gchr <- merge(gchr, gid, by = "seqid", all.x = T)
    gchr[,isBest := nGene == max(nGene) & !is.na(chromosome) & chromosome != "Unknown", by = "chromosome"]
    gchr[,chr := ifelse(isBest, chromosome, seqid)]
    setorder(gchr, -nbp)
    gchr <- subset(gchr, !duplicated(seqid))

    # rename sequences in gene gff with regions
    gff <- merge(gff, gchr[,c("seqid", "chr")], by = "seqid")
    uchrs <- uchrs[uchrs %in% unique(gff$chr)]

    # order genes
    gff[,chr := factor(chr, levels = uchrs)]
    gff$chr[is.na(gff$chr)] <- gff$seqid[is.na(gff$chr)]
    setorder(gff, chr, start, end, na.last = T)

    # read in the peptide fa
    if(verbose)
      cat("\tReading peptide fasta ... ")
    fa <- readAAStringSet(pepIn)

    # rename peptide headers
    tmp <- sapply(names(fa), function(x) strsplit(x, " ")[[1]][2])
    tmp2 <- substr(tmp, 7, nchar(tmp)-1)
    names(fa) <- tmp2
    nfa <- length(fa)

    # drop duplicates, keeping the longest
    fa <- fa[order(-width(fa))]
    fa <- fa[!duplicated(names(fa))]
    if(verbose)
      cat(sprintf("found %s / %s total and unique entries\n\tMerging fa and gff",
                  nfa, nrow(gff)))
    # join the two
    gff <- subset(gff, gene %in% names(fa))
    fa <- fa[gff$gene]
    if(verbose)
      cat(sprintf(" found %s matching entires\n", nrow(gff)))

    # write output
    gff <- with(gff, data.table(
      chr = chr,
      start = start,
      end = end,
      id = gene,
      strand = strand,
      ord = 1:length(gene)))

    fwrite(gff, file = gffOut, sep = "\t", quote = F)
    writeXStringSet(fa, filepath = pepOut, compress = "gzip")
  }

  rawGff <- gsParam$gffRaw
  rawPep <- gsParam$peptideRaw
  outGff <- file.path(gsParam$gff, sprintf("%s.gff.gz", genomeIDs))
  outPep <- file.path(gsParam$peptide, sprintf("%s.fa.gz", genomeIDs))
  verbose <- gsParam$verbose
  names(outPep) <- names(outGff) <- genomeIDs

  if(!all(file.exists(outGff)) | !all(file.exists(outPep)) | overwrite){
    for(i in genomeIDs){
      if(verbose)
        cat(sprintf("Parsing annotations: %s\n",i))
      ncbi_fun(
        gffIn = rawGff[i], gffOut = outGff[i],
        pepIn = rawPep[i], pepOut = outPep[i], verbose = verbose)
    }
    if(verbose)
      cat("Done!")
  }
  return(list(gff = outGff, peptide = outPep))
}

#' @title build matching gene fasta and coordinate datasets
#' @description
#' \code{match_gffFasta} parse fasta headers and gff3 attributes
#' @rdname format_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width
#' @export
match_gffFasta <- function(
  genomeIDs,
  path2rawGff3,
  path2rawPeptide,
  path2gff = file.path(getwd(), "genome", "gff", paste0(genomeIDs, ".gff.gz")),
  path2peptide = file.path(getwd(), "genome", "peptide", paste0(genomeIDs, ".fa.gz")),
  gffEntryType = "gene",
  gffIdColumn = "Name",
  headerEntryIndex = 4,
  headerSep = " ",
  headerStripText = "locus=",
  overwrite = FALSE,
  troubleshoot = FALSE,
  verbose = TRUE,
  nCores = detectCores()/2){

  ord <- id <- start <- end <- nbp <- NULL

  on.exit(expr = setDTthreads(getDTthreads()))

  setDTthreads(threads = nCores)
  verbose <- check_logicalArg(verbose)
  troubleshoot <- check_logicalArg(troubleshoot)
  overwrite <- check_logicalArg(overwrite)

  gffdone <- all(file.exists(path2gff))
  pepdone <- all(file.exists(path2peptide))
  if(!overwrite & gffdone & pepdone){
    warning("Output files exist; to re-compute, rerun with overwrite = TRUE\n")
    names(path2gff) <- genomeIDs
    names(path2peptide) <- genomeIDs
  }else{
    ##############################################################################
    # check arguments
    if(length(genomeIDs) < 1 | !is.character(genomeIDs))
      stop("must specify genomeIDs as a character vector of length >= 1\n")
    if(length(path2rawGff3) != length(genomeIDs))
      stop("length of path2rawGff3 must match length of genomeIDs\n")
    names(path2rawGff3) <- genomeIDs
    if(length(path2rawPeptide) != length(genomeIDs))
      stop("length of path2rawPeptide must match length of genomeIDs\n")
    names(path2rawPeptide) <- genomeIDs
    if(length(path2peptide) != length(genomeIDs))
      stop("length of path2peptide must match length of genomeIDs\n")
    names(path2peptide) <- genomeIDs
    if(length(path2gff) != length(genomeIDs))
      stop("length of path2gff must match length of genomeIDs\n")
    names(path2gff) <- genomeIDs

    if(length(headerEntryIndex) == 1)
      headerEntryIndex <- rep(headerEntryIndex, length(genomeIDs))
    if(length(headerEntryIndex) != length(genomeIDs))
      stop("headerEntryIndex must be a vector of length 1, or match length of genomeIDs\n")
    names(headerEntryIndex) <- genomeIDs

    if(length(headerSep) == 1)
      headerSep <- rep(headerSep, length(genomeIDs))
    if(length(headerSep) != length(genomeIDs))
      stop("headerSep must be a vector of length 1, or match length of genomeIDs\n")
    names(headerSep) <- genomeIDs

    if(length(headerStripText) == 1)
      headerStripText <- rep(headerStripText, length(genomeIDs))
    if(length(headerStripText) != length(genomeIDs))
      stop("headerStripText must be a vector of length 1, or match length of genomeIDs\n")
    names(headerStripText) <- genomeIDs

    if(length(gffEntryType) == 1)
      gffEntryType <- rep(gffEntryType, length(genomeIDs))
    if(length(gffEntryType) != length(genomeIDs))
      stop("gffEntryType must be a vector of length 1, or match length of genomeIDs\n")
    names(gffEntryType) <- genomeIDs

    if(length(gffIdColumn) == 1)
      gffIdColumn <- rep(gffIdColumn, length(genomeIDs))
    if(length(gffIdColumn) != length(genomeIDs))
      stop("gffIdColumn must be a vector of length 1, or match length of genomeIDs\n")
    names(gffIdColumn) <- genomeIDs

    if(!dir.exists(dirname(path2gff[1])))
      dir.create(path2gff[1], recursive = T)
    if(!dir.exists(dirname(path2peptide[1])))
      dir.create(path2peptide[1], recursive = T)

    ##############################################################################
    # -- loop through, parsing each genome ID
    if(verbose)
      cat("Parsing annotation files ...\n")
    for(i in genomeIDs){
      if(verbose)
        cat("\t",i, " ... \n\t\t", sep = "")
      gff <- parse_gff(
        path2rawGff3 = path2rawGff3[i],
        gffEntryType = gffEntryType[i],
        gffIdColumn = gffIdColumn[i],
        verbose = verbose,
        troubleshoot = troubleshoot)

      if(any(duplicated(gff$id))){
        cat("\n**FOUND DUPLICATE GFF ENTRIES. THERE COULD BE A PROBLEM***",
            "\nsubsetting to longest model for each gene\n")
        gff[,nbp := abs(end - start)]
        setorder(gff, -nbp)
        gff <- subset(gff, !duplicated(id))
        setkey(gff, ord)
      }

      if(verbose)
        cat("\t\t")
      pep <- parse_faHeader(
        path2rawFasta = path2rawPeptide[i],
        headerEntryIndex = headerEntryIndex[i],
        headerSep = headerSep[i],
        headerStripText = headerStripText[i],
        verbose = verbose,
        troubleshoot = troubleshoot)

      if(any(duplicated(names(pep)))){
        cat("\n**FOUND DUPLICATE PEPTIDE FASTA HEADERS. IS THIS JUST A PRIMARY TRANSCRIPT FILE?***",
            "\nsubsetting to longest model for each gene\n")
        pepSize <- width(pep)
        pep <- pep[order(-pepSize)]
        pep <- pep[!duplicated(names(pep))]
      }
      writeXStringSet(
        pep,
        filepath = path2peptide[i],
        compress = "gzip")
      if(verbose)
        cat("\t\tPeptides written to", path2peptide[i], "\n\t\tChecking overlaps: ")
      gff <- subset(gff, id %in% names(pep))
      if(verbose)
        cat(nrow(gff),"gff matches; ")
      pep <- pep[gff$id]
      if(verbose)
        cat(length(pep),"peptide matches\n")

      fwrite(
        gff,
        file = path2gff[i],
        sep = "\t")
      if(verbose)
        cat("\t\tGff written to", path2gff[i], "\n")
    }
    cat("\tDone!\n")
  }

  return(list(
    gff = path2gff,
    peptide = path2peptide))
}

#' @title convert gff3-formatted annotation files for genespace
#' @description
#' \code{parse_gff} parse the gff3-formatted annotation into genespace-readable
#' data.table/csv file. This has a unique ID for each gene that matches the
#' fasta annotation entries exactly.
#' @rdname format_annotations
#' @import data.table
#' @export
parse_gff <- function(
  path2rawGff3,
  gffEntryType,
  gffIdColumn,
  verbose = TRUE,
  troubleshoot = FALSE){

  ord <- end <- start <- chr <- tmp2 <- tmp1 <- type <- NULL

  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  if(length(path2rawGff3) > 1){
    warning("can only parse one gff, 1st the first\n")
    path2rawGff3 <- path2rawGff3[1]
  }
  if(length(gffEntryType) > 1){
    warning("can only parse one gff, taking the 1st gffEntryType\n")
    gffEntryType <- gffEntryType[1]
  }
  if(length(gffIdColumn) > 1){
    warning("can only parse one gff, taking the 1st gffIdColumn\n")
    gffIdColumn <- gffIdColumn[1]
  }
  if(!file.exists(path2rawGff3))
    stop("cant find", path2rawGff3, "\n")
  verbose <- check_logicalArg(verbose)
  troubleshoot <- check_logicalArg(troubleshoot)

  # -- Import gff into R env
  if(verbose)
    cat("Importing gff ... ")
  gff <- data.table(data.frame(rtracklayer::readGFF(
    filepath = path2rawGff3,
    filter = list(type = gffEntryType),
    tags = c(gffIdColumn))))
  if(troubleshoot)
    print(head(gff, 10))
  if(verbose)
    cat("found", nrow(gff), "gff entires")

  # -- subset to the type of entry desired and reformat
  if(sum(gff$type == gffEntryType) < 1)
    stop("cant find", gffEntryType, "in the type (3rd) column\n")
  if(!gffIdColumn %in% colnames(gff))
    stop("cant find column named ", gffIdColumn,
         "\n\tAvailable column names are:", paste(colnames(gff), collapse = ", "),"\n")
  gffSub <- subset(gff, type == gffEntryType)
  gff <- with(gffSub, data.table(
    chr = seqid,
    start = start,
    end = end,
    id = get(gffIdColumn),
    strand = strand))

  # -- sort by chromosome number (strip text), then chr size, the position
  gff[,tmp1 := as.numeric(gsub("[^0-9]", "", chr))]
  gff[,tmp2 := .N, by = "chr"]
  setorder(gff, tmp1, -tmp2, start, end)
  gff[,ord := 1:.N]
  gff[,`:=`(tmp1 = NULL, tmp2 = NULL)]
  if (troubleshoot)
    print(head(gff, 10))
  if(verbose)
    cat(", and", nrow(gffSub), gffEntryType, "entries\n")
  return(gff)
}

#' @title convert fasta files for genespace
#' @description
#' \code{parse_faHeader} parse fasta headers to gene gff entries
#' @rdname format_annotations
#' @import data.table
#' @importFrom Biostrings readAAStringSet readDNAStringSet
#' @importFrom utils head
#' @export
parse_faHeader <- function(
  path2rawFasta,
  headerEntryIndex,
  headerSep,
  headerStripText,
  verbose = TRUE,
  troubleshoot = FALSE){

  # -- argument checking
  if(length(path2rawFasta) > 1){
    warning("can only parse one fasta, taking the 1st\n")
    path2rawFasta <- path2rawFasta[1]
  }
  if(length(headerEntryIndex) > 1){
    warning("can only parse one fasta, taking the 1st headerEntryIndex\n")
    headerEntryIndex <- headerEntryIndex[1]
  }
  if(length(headerSep) > 1) {
    warning("can only parse one fasta, taking the 1st headerSep\n")
    headerSep <- headerSep[1]
  }
  if(length(headerStripText) > 1){
    warning("can only parse one fasta, taking the 1st headerStripText\n")
    headerStripText <- headerStripText[1]
  }
  verbose <- check_logicalArg(verbose)
  troubleshoot <- check_logicalArg(troubleshoot)

  # -- Check if the fasta is a peptide or not (if necessary)
  nu <- check_isPeptideFasta(path2rawFasta)

  # -- Read fasta file in
  if (verbose)
    cat("Importing fasta ... ")
  fa <- readAAStringSet(path2rawFasta)
  if (troubleshoot)
    print(head(names(fa), 10))

  # parse the header, following use-specified rules
  if (verbose)
    cat("found", length(fa), "fasta entires\n")
  nparse <- tstrsplit(names(fa), headerSep)[[headerEntryIndex]]
  if (!is.null(headerStripText)) {
    nparse <- gsub(headerStripText, "", nparse)
  }
  names(fa) <- nparse
  if (troubleshoot)
    print(head(names(fa), 10))
  return(fa)
}

