#' @title Move and organize raw genome data
#'
#' @description
#' \code{convert_genomes} Take phytozome-formatted genomes, rename and reformat them
#' for the GENESPACE pipeline
#'
#' @param genomeIDs genome identifiers
#' @param directory file path to the base directory. Must contain two subdirectories:
#' raw_annotations and raw_assemblies. The directory structures of these are outlined
#' below.
#' @param cds_str character string to identify primary CDS fasta file
#' @param peptide_str character string to identify primary peptide fasta file
#' @param gff_str character string to identify gff annotation file
#' @param parse_fastaHeader.FUN The function to be used to parse the fasta headers.
#' @param verbose should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return The function does not return anything to the R console.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet
#' @export
convert_genomes <- function(genomeIDs,
                            directory,
                            peptide.only = F,
                            peptide_str = "protein_primaryTranscriptOnly",
                            cds_str = "cds_primaryTranscriptOnly",
                            gff_str = "gene.gff3",
                            parse_fastaHeader.FUN = function(y)
                              strsplit(gsub(".*locus=", "", y)," ")[[1]][1],
                            verbose = T,
                            ...){

  raw_assembly.dir <- file.path(directory, "raw_assemblies")
  raw_annot.dir <- file.path(directory, "raw_annotations")

  # 0. Check input
  if (!all(genomeIDs %in% dir(raw_annot.dir)))
    stop("all specified genomeIDs must be folder names in raw_annot.dir\n")

  if(peptide.only){
    strs <- c(peptide_str,
              gff_str)
    names(strs) <- c("peptide","gff")
  }else{
    strs <- c(peptide_str,
              cds_str,
              gff_str)
    names(strs) <- c("peptide","cds","gff")
  }


  file.lists <- sapply(names(strs), USE.NAMES = T, simplify = F, function(j){
    sapply(genomeIDs, function(i){
      list.files(file.path(raw_annot.dir, i),
                 pattern = strs[j], full.names = T)
    })
  })

  if (any(sapply(file.lists, function(x) any(is.null(x)))))
    stop("some annotation files are missing\n")

  if (verbose)
    cat("Making genome directories\n")

  # 1. make the directories
  input.dir <- file.path(directory, "genome")
  if (file.exists(input.dir))
    unlink(input.dir, recursive = T)
  dir.create(input.dir)

  if(peptide.only){
    ftypes <- c("peptide",
                "gff")
  }else{
    ftypes <- c("peptide",
                "cds",
                "gff",
                "assembly")
  }


  subdirs <- sapply(ftypes, USE.NAMES = T, simplify = F, function(x){
    fp <- file.path(input.dir, x)
    if (dir.exists(fp))
      nu <- unlink(fp, recursive = T)
    if (!dir.exists(fp))
      dir.create(fp)
    return(fp)
  })

  # 2. Move files into directories, uncompressing if necessary
  if (verbose)
    cat("Moving and unzipping annotation files\n")
  for (i in names(file.lists)) {
    for (j in names(file.lists[[i]])) {
      x <- file.lists[[i]][j]
      y <- subdirs[[i]]
      outname <- file.path(y, ifelse(grepl("gff", x),
                                     paste0(j, ".gff3"),
                                     paste0(j, ".fa")))
      if (grepl(".gz$", x)) {
        system(paste("gunzip -c", x, ">", outname))
      }else{
        system(paste("cp", x, outname))
      }
    }
  }

  # 3. Index genomes
  if(!peptide.only){
    if (verbose)
      cat("Moving and indexing assembly fastas\n")

    assem.dir <- subdirs[["assembly"]]
    ass.file <- file.path(raw_assembly.dir,paste0(genomeIDs,".fa"))
    nu <- file.copy(ass.file,
                    assem.dir)

    fais <- list.files(assem.dir,
                       pattern = ".fai$",
                       full.names = T)
    if (length(fais) > 0)
      nu <- file.remove(fais)

    assem.fas <- list.files(assem.dir,
                            pattern = ".fa$",
                            full.names = T)
    if (Sys.which("samtools") == "") {
      warning("samtools is not in the path, keep in mind that indexed assemblies are required for pseudogene inference.\n")
    }else{
      for (i in assem.fas)
        system(paste("samtools faidx", i))
    }
  }


  # 4. re-name annotation fastas with gene name (not model).
  if (verbose)
    cat("Renaming annotation fasta headers\n")
  if(!peptide.only){
    tmp <- parse_fastaHeader(fasta.dir = subdirs$cds,
                             is.peptide = F,
                             verbose = F,
                             parse_fastaHeader.FUN = parse_fastaHeader.FUN)
  }

  tmp <- parse_fastaHeader(fasta.dir = subdirs$peptide,
                           is.peptide = T,
                           verbose = F,
                           parse_fastaHeader.FUN = parse_fastaHeader.FUN)

  if (verbose)
    cat("Done!\n")
}
