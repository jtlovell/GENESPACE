#' @title Functions to convert from raw to genespace annotation formats
#' @description
#' \code{parse_annotations} Format gff and fasta for genespace. The two parse
#' functions (parse_gff & parse_faHeader) take a single file and return a
#' formatted R object. format_annotations is a wrapper for these two that reads
#' and writes the full genespace-formatted files.
#'
#' @name parse_annotations
#' @param gsParam a list containing all parameters for a GENESPACE run. See
#' init_genespace
#' @param genomeIDs character vector of the same length as paths. Names for
#' each genome.
#' @param path2rawGff3 character vector coercible to file.paths. This points to
#' the raw gff3-formatted annotation (.gz or gff3)
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
#' @param gffStripText like headerStripText but for the gff geneID entry
#' @param troubleshoot logical, should the raw and parsed files be printed?
#' @param path2rawFasta character string coercible to a file path, pointing to
#' the location of the unparsed fasta file.
#' @param overwrite logical, should existing files be overwritten? If FALSE,
#' and all files are present, just returns a named vector of gff/pep file
#' locations
#' @param verbose logical, should updates be printed to the console?

#' @note \code{parse_annotations} is a generic name for the functions documented.
#' \cr
#' If called, \code{parse_annotations} returns its own arguments.
#'

#' @title parse phytozome-formatted annotations
#' @description
#' \code{parse_phytozome} parse phytozome-formatted annotations
#' @rdname parse_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width readAAStringSet AAStringSet
#' @export
parse_phytozome <- function(gsParam, overwrite = F, genomeIDs = NULL){


  pytz_fun <- function(gffIn, gffOut, pepIn, pepOut, verbose, minPepLen){
    Name <- seqid <- start <- end <- NULL
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

    # cull out "."s very short sequences
    fa <- AAStringSet(gsub(".","",fa, fixed = T))
    fa <- fa[width(fa) > minPepLen]

    if(verbose)
      cat(sprintf("found %s / %s total and unique entries\n\tMerging fa and gff",
                  nfa, nrow(gff)))

    # join the two
    gff <- subset(gff, Name %in% names(fa))
    fa <- fa[gff$Name]
    gff <- subset(gff, Name %in% names(fa))

    if(nrow(gff) == 0){
      stop("PARSING FAILED - is this a phytozome-formatted annotation?\n")
    }else{
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
      writeXStringSet(fa, filepath = pepOut)
    }
  }

  on.exit(expr = setDTthreads(getDTthreads()))
  setDTthreads(threads = gsParam$params$nCores)

  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(!all(genomeIDs %in% gsParam$genomes$genomeIDs)){
    wh <- which(!genomeIDs %in% gsParam$genomes$genomeIDs)
    stop(sprintf("specified %s genomes not in gsParam\n", paste(genomeIDs[wh], collapse = ", ")))
  }

  for(i in genomeIDs){
    if(file.exists(gsParam$paths$gff[i]) && file.exists(gsParam$paths$peptide[i]) && !overwrite){
      cat(sprintf("parsed annotations for %s exist and !overwrite, skipping\n", i))
    }else{
      if(gsParam$params$verbose)
        cat(sprintf("Parsing annotations: %s\n",i))
      pytz_fun(
        gffIn = gsParam$paths$rawGff[i],
        gffOut = gsParam$paths$gff[i],
        pepIn = gsParam$paths$rawPeptide[i],
        pepOut = gsParam$paths$peptide[i],
        verbose = gsParam$params$verbose,
        minPepLen = gsParam$params$minPepLen)
    }
  }
}


#' @title parse NCBI-formatted annotations
#' @description
#' \code{parse_ncbi} parse NCBI-formatted annotations
#' @rdname parse_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width readAAStringSet AAStringSet
#' @export
parse_ncbi <- function(gsParam, overwrite = F, genomeIDs = NULL){

  gene_biotype <- gene <- seqid <- chromosome <- end <- start <- isBest <- NULL
  chr <- nbp <- nGene <- NULL

  ncbi_fun <- function(gffIn, gffOut, pepIn, pepOut, verbose, minPepLen){
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
    fa <- AAStringSet(gsub(".","",fa, fixed = T))
    fa <- fa[width(fa) > minPepLen]

    if(verbose)
      cat(sprintf("found %s / %s total and unique entries\n\tMerging fa and gff",
                  nfa, nrow(gff)))
    # join the two
    gff <- subset(gff, gene %in% names(fa))
    fa <- fa[gff$gene]

    if(nrow(gff) == 0){
      stop("PARSING FAILED - is this a ncbi-formatted annotation?\n")
    }else{
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
      writeXStringSet(fa, filepath = pepOut)
    }
  }

  on.exit(expr = setDTthreads(getDTthreads()))
  setDTthreads(threads = gsParam$params$nCores)
  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(!all(genomeIDs %in% gsParam$genomes$genomeIDs)){
    wh <- which(!genomeIDs %in% gsParam$genomes$genomeIDs)
    stop(sprintf("specified %s genomes not in gsParam\n", paste(genomeIDs[wh], collapse = ", ")))
  }

  for(i in genomeIDs){
    if(file.exists(gsParam$paths$gff[i]) && file.exists(gsParam$paths$peptide[i]) && !overwrite){
      cat(sprintf("parsed annotations for %s exist and !overwrite, skipping\n", i))
    }else{
      if(gsParam$params$verbose)
        cat(sprintf("Parsing annotations: %s\n",i))
      ncbi_fun(
        gffIn = gsParam$paths$rawGff[i],
        gffOut = gsParam$paths$gff[i],
        pepIn = gsParam$paths$rawPeptide[i],
        pepOut = gsParam$paths$peptide[i],
        verbose = gsParam$params$verbose,
        minPepLen = gsParam$params$minPepLen)
    }
  }
}

#' @title build matching gene fasta and coordinate datasets
#' @description
#' \code{parse_annotations} parse fasta headers and gff3 attributes
#' @rdname parse_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width readAAStringSet AAStringSet
#' @export
parse_annotations <- function(
  genomeIDs = NULL,
  gsParam,
  gffEntryType = "gene",
  gffIdColumn = "Name",
  gffStripText = "",
  headerEntryIndex = 4,
  headerSep = " ",
  headerStripText = "locus=",
  overwrite = FALSE,
  troubleshoot = FALSE){

  nbp <- end <- start <- id <- ord <- NULL

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(!all(genomeIDs %in% gsParam$genomes$genomeIDs)){
    wh <- which(!genomeIDs %in% gsParam$genomes$genomeIDs)
    stop(sprintf("specified %s genomes not in gsParam\n", paste(genomeIDs[wh], collapse = ", ")))
  }

  on.exit(expr = setDTthreads(getDTthreads()))
  setDTthreads(threads = gsParam$params$nCores)

  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to parse gff files, install rtracklayer from bioconductor\n")

  ##############################################################################
  # check arguments
  troubleshoot <- check_logicalArg(troubleshoot)
  overwrite <- check_logicalArg(overwrite)
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

  if(length(gffStripText) == 1)
    gffStripText <- rep(gffStripText, length(genomeIDs))
  if(length(gffStripText) != length(genomeIDs))
    stop("gffStripText must be a vector of length 1, or match length of genomeIDs\n")
  names(gffStripText) <- genomeIDs

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

  ##############################################################################
  # -- loop through, parsing each genome ID
  if(gsParam$params$verbose)
    cat("Parsing annotation files ...\n")

  for(i in genomeIDs){
    if(gsParam$params$verbose)
      cat("\t",i, " ... \n\t\t", sep = "")

    if(file.exists(gsParam$paths$gff[i]) && file.exists(gsParam$paths$peptide[i]) && !overwrite){
      cat(sprintf("parsed annotations for %s exist and !overwrite, skipping\n", i))
    }else{
      gff <- parse_gff(
        path2rawGff3 = gsParam$paths$rawGff[i],
        gffEntryType = gffEntryType[i],
        gffIdColumn = gffIdColumn[i],
        gffStripText = gffStripText[i],
        verbose = gsParam$params$verbose,
        troubleshoot = troubleshoot)

      if(any(duplicated(gff$id))){
        cat("\n**FOUND DUPLICATE GFF ENTRIES. THERE COULD BE A PROBLEM***",
            "\nsubsetting to longest model for each gene\n")
        gff[,nbp := abs(end - start)]
        setorder(gff, -nbp)
        gff <- subset(gff, !duplicated(id))
        setkey(gff, ord)
      }

      if(gsParam$params$verbose)
        cat("\t\t")
      pep <- parse_faHeader(
        path2rawFasta = gsParam$paths$rawPeptide[i],
        headerEntryIndex = headerEntryIndex[i],
        headerSep = headerSep[i],
        headerStripText = headerStripText[i],
        verbose = gsParam$params$verbose,
        troubleshoot = troubleshoot)

      if(any(duplicated(names(pep)))){
        cat("\n**FOUND DUPLICATE PEPTIDE FASTA HEADERS. IS THIS JUST A PRIMARY TRANSCRIPT FILE?***",
            "\nsubsetting to longest model for each gene\n")
        pepSize <- width(pep)
        pep <- pep[order(-pepSize)]
        pep <- pep[!duplicated(names(pep))]
      }
      pep <- AAStringSet(gsub(".","",pep, fixed = T))
      pep <- pep[width(pep) > gsParam$params$minPepLen]

      gff <- subset(gff, id %in% names(pep))
      pep <- pep[gff$id]
      if(gsParam$params$verbose)
        cat("\t\t",nrow(gff)," gff-peptide matches\n", sep = "")

      if(nrow(gff) == 0){
        if(troubleshoot)
          stop("PARSING FAILED - try different parameters\n")
        if(!troubleshoot)
          stop("PARSING FAILED - run with troubleshoot = TRUE to see where the issues might be\n")
      }else{
        writeXStringSet(
          pep,
          filepath = gsParam$paths$peptide[i])

        fwrite(
          gff,
          file = gsParam$paths$gff[i],
          sep = "\t")
      }
      if(gsParam$params$verbose)
        cat("\tDone!\n")
    }
  }
}

#' @title convert gff3-formatted annotation files for genespace
#' @description
#' \code{parse_gff} parse the gff3-formatted annotation into genespace-readable
#' data.table/csv file. This has a unique ID for each gene that matches the
#' fasta annotation entries exactly.
#' @rdname parse_annotations
#' @import data.table
#' @export
parse_gff <- function(
  path2rawGff3,
  gffEntryType,
  gffIdColumn,
  gffStripText,
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
  if(length(gffStripText) > 1){
    warning("can only parse one gff, taking the 1st gffStripText\n")
    gffStripText <- gffStripText[1]
  }
  if(!file.exists(path2rawGff3))
    stop("cant find", path2rawGff3, "\n")

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
    id = gsub(gffStripText, "", get(gffIdColumn)),
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
#' @rdname parse_annotations
#' @import data.table
#' @importFrom Biostrings writeXStringSet width readAAStringSet AAStringSet
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

