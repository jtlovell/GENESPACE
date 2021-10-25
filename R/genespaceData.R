#' @title Test data for GENESPACE examples
#'
#' @description
#' \code{genespaceData} Raw data parsing, viewing and exploration
#'
#' @name genespaceData
#'
#' @param writeDir file.path to move the raw files to. Must be an empty or non-
#' existent path within an existing directory.
#'
#'
#' @format A protien (i.e. translated cds) fasta file with headers matching
#' the locus entry in the gff, or a simplified gff3-formatted annotation text
#' file. The gff3 annotation columns are as follows:
#' \describe{
#'   \item{seqid}{name of the chromosome or scaffold}
#'   \item{source}{name of the program that generated this feature, or the data source (database or project name)}
#'   \item{type}{type of feature. Subset to 'gene' here}
#'   \item{start}{Start position of the feature, with sequence numbering starting at 1.}
#'   \item{end}{End position of the feature, with sequence numbering starting at 1.}
#'   \item{score}{set to NA, not used}
#'   \item{strand}{defined as + (forward) or - (reverse).}
#'   \item{phase}{set to NA, not used}
#'   \item{attributes}{parsed to just the gene name, appended with locus=}
#' }
#'
#' @details ...
#'
#' @return ...
#' @examples
#' \dontrun{
#' make_exampleData(writeDir = "inst/extdata")
#' }
#' @note \code{genespaceData} is a generic name for the functions documented.
#' \cr
#' If called, \code{genespaceData} returns its own arguments.
#'
#' @title make_exampleData
#' @description
#' \code{make_exampleData} make_exampleData
#' @rdname genespaceData
#' @import data.table
#' @importFrom Biostrings readAAStringSet writeXStringSet AAStringSet width
#' @importFrom rtracklayer readGFF
#' @export
make_exampleData <- function(writeDir){

  # -- make sure the writeDir is empty or non-existent in an existing parent dir
  if(dir.exists(writeDir))
    if(length(dir(writeDir)) > 0)
      stop(
        writeDir, "exists and is not empty\n")
  if(dir.exists(writeDir))
    if(!dir.exists(dirname(writeDir)))
      stop(
        dirname(writeDir), "doesnt exist\n")

  # -- get the raw data
  rawFileLoc <- download_rawData(
    writeDir = "/Users/jlovell/Desktop/GENESPACE/inst/extdata")

  gffs <- sapply(names(ftps), function(i)
    file.path(wd, sprintf("%s_%s_gene.gff.gz", i, names(ftps[[i]]))))
  peps <- sapply(names(ftps), function(i)
    file.path(wd, sprintf("%s_%s_pep.fa.gz", i, names(ftps[[i]]))))


  gffl <- lapply(names(rawFileLoc), function(i){
    cat(i)
    gffLoc <- rawFileLoc[[i]]$gff
    pepLoc <- rawFileLoc[[i]]$pep
    gff <- data.table(data.frame(readGFF(
      filepath = gffLoc,
      filter = list(type = "gene"),
      tags = c("gene", "gene_biotype"))))
    gff <- subset(gff, gene_biotype == "protein_coding")
    gff <- subset(gff, !duplicated(gene))
    gff[,`:=`(seqid = as.character(seqid), gene = as.character(gene))]
    gchr <- gff[,list(nGene = .N), by = "seqid"]

    # -- read in the region gff
    chrIDs <- data.table(data.frame(readGFF(
      filepath = gffLoc,
      filter = list(type = "region"),
      tags = c("chromosome"))))
    chrIDs[,`:=`(seqid = as.character(seqid), chromosome = as.character(chromosome))]

    # -- choose largest regions / chromosome as rep
    gid <- chrIDs[,list(nbp = end - start), by = c("chromosome","seqid")]
    gchr <- merge(gchr, gid, by = "seqid", all.x = T)
    gchr[,isBest := nGene == max(nGene) & !is.na(chromosome) & chromosome != "Unknown", by = "chromosome"]
    gchr[,chr := ifelse(isBest, chromosome, seqid)]
    setorder(gchr, -nbp)
    gchr <- subset(gchr, !duplicated(seqid))

    # -- rename sequences in gene gff with regions
    gff <- merge(gff, gchr[,c("seqid", "chr")], by = "seqid")
    uchrs <- c(1,2,"2a","2b",3)
    gff[,seqid := chr]
    gff <- subset(gff[,1:9], seqid %in% uchrs)
    cat("", nrow(gff), "gff entries ... ")

    # -- subset the peptide
    fa <- readAAStringSet(pepLoc)

    # rename peptide headers
    tmp <- sapply(names(fa), function(x) strsplit(x, " ")[[1]][2])
    tmp2 <- substr(tmp, 7, nchar(tmp)-1)
    names(fa) <- tmp2
    nfa <- length(fa)

    # drop duplicates, keeping the longest
    fa <- fa[order(-width(fa))]
    fa <- fa[!duplicated(names(fa))]
    fa <- AAStringSet(gsub(".","",fa, fixed = T))
    fa <- fa[width(fa) > 20]

    # -- subset the gff to genes in fa
    gff <- subset(gff, gene %in% names(fa))
    fa <- fa[gff$gene]
    cat(nrow(gff), "merged\n")
    # -- write the subsetted gff
    gff[[9]] <- paste0("locus=",gff[[9]])
    fwrite(
      gff,
      file = gffLoc,
      sep = "\t",
      quote = F,
      col.names = F)

    # -- write the peptides
    writeXStringSet(
      fa,
      filepath = pepLoc,
      compress = TRUE)
  })
  cat("Done!")
}

#' @title download_rawData
#' @description
#' \code{download_rawData} download_rawData
#' @rdname genespaceData
#' @importFrom utils download.file
#' @export
download_rawData <- function(writeDir){
  wd <- writeDir
  ftps <- find_rawData()
  # wd <- "~/Desktop/GENESPACE/inst/extdata"
  if(!dir.exists(wd))
    dir.create(wd, recursive = T)
  out <- sapply(names(ftps), USE.NAMES = T, simplify = F, function(i){
    cat(i,": peptides ... ")
    pepf <- file.path(wd, sprintf("%s_%s_pep.fa.gz", i, names(ftps[[i]])))
    gfff <- file.path(wd, sprintf("%s_%s_gene.gff.gz", i, names(ftps[[i]])))
    download.file(
      url = ftps[[i]][[1]][[1]],
      destfile = pepf,
      method = "wget",
      quiet = T)
    cat("Done! gff: ...")
    download.file(
      url = ftps[[i]][[1]][[2]],
      destfile = gfff,
      quiet = T,
      method = "wget")
    cat(" Done!\n")
    return(list(pep = pepf, gff = gfff))
  })
  return(out)
}

#' @title location_ofRawData
#' @description
#' \code{location_ofRawData} location_ofRawData
#' @rdname genespaceData
#' @export
find_rawData <- function(){
  ftps <- list(
    human = list(GRCh38.p13 = list(
      pep = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_translated_cds.faa.gz",
      gff = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz")),
    chimp = list(Clint_PTRv2 = list(
      pep = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Pan_troglodytes/latest_assembly_versions/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_translated_cds.faa.gz",
      gff = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Pan_troglodytes/latest_assembly_versions/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_genomic.gff.gz")),
    rhesus = list(Mmul_10 = list(
      pep = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_translated_cds.faa.gz",
      gff = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gff.gz")))
  return(ftps)
}

#' @title find_exampleData
#' @description
#' \code{find_exampleData} find_exampleData
#' @rdname genespaceData
#' @export
find_exampleData <- function(){
  paths <- list(
    human = list(
      gff = system.file(
        "extdata",
        "human_GRCh38.p13_gene.gff.gz",
        package = "GENESPACE",
        mustWork = TRUE),
      pep = system.file(
        "extdata",
        "human_GRCh38.p13_pep.fa.gz",
        package = "GENESPACE",
        mustWork = TRUE)),

    chimp = list(
      gff = system.file(
        "extdata",
        "chimp_Clint_PTRv2_gene.gff.gz",
        package = "GENESPACE",
        mustWork = TRUE),
      pep = system.file(
        "extdata",
        "chimp_Clint_PTRv2_pep.fa.gz",
        package = "GENESPACE",
        mustWork = TRUE)),

    rhesus = list(
      gff = system.file(
        "extdata",
        "rhesus_Mmul_10_gene.gff.gz",
        package = "GENESPACE",
        mustWork = TRUE),
      pep = system.file(
        "extdata",
        "rhesus_Mmul_10_pep.fa.gz",
        package = "GENESPACE",
        mustWork = TRUE)))
  return(paths)
}

#' @title make_exampleDataDir
#' @description
#' \code{find_exampleData} find_exampleData
#' @rdname genespaceData
#' @export
make_exampleDataDir <- function(writeDir){

  # -- find the data within the R package
  exDat <- find_exampleData()

  # -- check the write directory
  if(dir.exists(writeDir))
    if(length(dir(writeDir)) > 0)
      stop(writeDir, "exists and is not empty\n")
  if(dir.exists(writeDir))
    if(!dir.exists(dirname(writeDir)))
      stop(dirname(writeDir), "doesnt exist\n")
  if(!dir.exists(writeDir))
    dir.create(writeDir)

  # -- make a raw genome dir in write dir
  fps <- file.path(writeDir, "rawGenomes", names(exDat), names(exDat), "annotation")
  names(fps) <- names(exDat)
  nu <- sapply(fps, dir.create, recursive = T)

  # -- copy files into new loc
  for(i in names(exDat))
    nu <- file.copy(from = unlist(exDat[[i]]), to = fps[i])

  cat("Done! Run GENESPACE with rawGenomeRepo =", file.path(writeDir, "rawGenomes"), "\n")
}

#' Dataset human chr1-3 peptide fasta
#'
#' @name human_GRCh38.p13_pep.fa.gz
#' @section human_GRCh38.p13_pep.fa.gz
#' @source \url{https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_translated_cds.faa.gz}
#' This data is used in examples throughout genespace help files.


#' Dataset chimpanzee chr1-3 peptide fasta
#'
#' @name chimp_Clint_PTRv2_pep.fa.gz
#' @section chimp_Clint_PTRv2_pep.fa.gz
#' @source \url{https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Pan_troglodytes/latest_assembly_versions/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_translated_cds.faa.gz}
#' This data is used in examples throughout genespace help files.


#' Dataset rhesus chr1, 2A, 2B, 3 peptide fasta
#'
#' @name rhesus_Mmul_10_pep.fa.gz
#' @section rhesus_Mmul_10_pep.fa.gz
#' @source \url{https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_translated_cds.faa.gz}
#' This data is used in examples throughout genespace help files.

#' Dataset human chr1-3 gff3
#'
#' @name human_GRCh38.p13_gene.gff.gz
#' @section human_GRCh38.p13_gene.gff.gz
#' @source \url{https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz}
#' This data is used in examples throughout genespace help files.


#' Dataset chimpanzee chr1-3 gff3
#'
#' @name chimp_Clint_PTRv2_gene.gff.gz
#' @section chimp_Clint_PTRv2_gene.gff.gz
#' @source \url{https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Pan_troglodytes/latest_assembly_versions/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_genomic.gff.gz}
#' This data is used in examples throughout genespace help files.


#' Dataset rhesus chr1, 2A, 2B, 3 gff3
#'
#' @name rhesus_Mmul_10_gene.gff.gz
#' @section rhesus_Mmul_10_gene.gff.gz
#' @source \url{https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gff.gz}
#' This data is used in examples throughout genespace help files.


