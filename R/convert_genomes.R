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
#' @param transcript_str character string to identify primary transcript fasta file
#' @param cds_str character string to identify primary CDS fasta file
#' @param peptide_str character string to identify primary peptide fasta file
#' @param gff_str character string to identify gff annotation file
#' @param parse_fastaHeader.FUN The function to be used to parse the fasta headers.
#' @param verbose should updates be printed?
#' @param ... Not currently in use
#' @details Not required, but makes the formatting for the program a bit easier.
#' The raw_assembly subdirectory MUST contain a single fasta file for each
#' genomeID, named $genomeID.fa.
#' The raw_annotation subdirectory MUST contain a single subdirectory for each
#' genomeID, with a label matching the genomeID exactly.
#' Each raw_annotation/genomeID subdirectory MUST contain four files: cds.fasta
#' peptide.fasta, annotation.gff3 and transcript.fasta. These can be named whatever,
#' but must contain a string that can be identified in the specified parameters of
#' this function
#' @return nothing
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet
#' @export
convert_genomes = function(genomeIDs,
                           directory,
                           transcript_str = "transcript_primaryTranscriptOnly",
                           peptide_str = "protein_primaryTranscriptOnly",
                           cds_str = "cds_primaryTranscriptOnly",
                           gff_str = "gene.gff3",
                           parse_fastaHeader.FUN = function(y)
                             strsplit(gsub(".*locus=","",y)," ")[[1]][1],
                           verbose = T,
                           ...){

  parse_fastaHeader <- function(fasta.dir,
                                is.peptide = T,
                                pattern = "fa",
                                verbose = T,
                                parse_fastaHeader.FUN,
                                ...){

    files <- list.files(fasta.dir,
                        pattern = pattern,
                        full.names = T)

    if(verbose)
      cat("Renaming fasta headers ...\n")
    ss <- lapply(files, function(i){
      if(verbose)
        cat("...", i, "\n\t")
      if(is.peptide){
        x <- readAAStringSet(i)
      }else{
        x <- readDNAStringSet(i)
      }
      if(verbose)
        cat("original names (e.g.):",
            names(x)[1])
      names(x) <- sapply(names(x), parse_fastaHeader.FUN)
      if(verbose)
        cat("\n\tparsed names (e.g.):",
            names(x)[1],"\n")
      writeXStringSet(x, filepath = i)
    })
  }

  raw_assembly.dir = file.path(directory, "raw_assemblies")
  raw_annot.dir = file.path(directory, "raw_annotations")

  # 0. Check input
  if(!all(genomeIDs %in% dir(raw_annot.dir))){
    stop("all specified genomeIDs must be folder names in raw_annot.dir\n")
  }

  strs <- c(transcript_str,
            peptide_str,
            cds_str,
            gff_str)
  names(strs)<-c("transcript","peptide","cds","gff")
  file.lists <- sapply(names(strs), USE.NAMES = T, simplify = F, function(j){
    sapply(genomeIDs, function(i){
      list.files(file.path(raw_annot.dir, i), pattern = strs[j], full.names = T)
    })
  })

  if(any(sapply(file.lists, function(x) any(is.null(x))))){
    stop("some annotation files are missing\n")
  }

  if(verbose)
    cat("Making genome directories\n")

  # 1. make the directories
  input.dir = file.path(directory, "genome")
  if(!file.exists(input.dir)){
    system(paste("mkdir", input.dir))
  }else{
    system(paste("rm -r", input.dir))
    system(paste("mkdir", input.dir))
  }
  ftypes <- c("transcript",
              "peptide",
              "cds",
              "gff",
              "assembly")

  subdirs <- sapply(ftypes, USE.NAMES = T, simplify = F, function(x){
    fp <- file.path(input.dir, x)
    if(file.exists(fp)){
      system(paste("rm -r", fp))
    }
    system(paste("mkdir", fp))
    return(fp)
  })

  # 2. Move files into directories, uncompressing if necessary
  if(verbose)
    cat("Moving and unzipping annotation files\n")
  for(i in names(file.lists)){
    for(j in names(file.lists[[i]])){
      x <- file.lists[[i]][j]
      y <- subdirs[[i]]
      outname <- file.path(y, ifelse(grepl("gff", x),
                                     paste0(j, ".gff3"),
                                     paste0(j, ".fa")))
      if(grepl(".gz$", x)){
        system(paste("gunzip -c", x, ">", outname))
      }else{
        system(paste("cp", x, outname))
      }
    }
  }

  # 3. Index genomes
  if(verbose)
    cat("Moving and indexing assembly fastas\n")

  assem.dir <- subdirs[["assembly"]]
  system(paste("cp",
               file.path(raw_assembly.dir, "*"),
               assem.dir))

  if(any(grepl("*.fai", list.files(assem.dir)))){
    system(paste("rm", file.path(assem.dir, "*.fai")))
  }
  assem.fas <- list.files(assem.dir, full.names = T)
  for(i in assem.fas){
    system(paste("samtools faidx", i))
  }

  # 4. re-name annotation fastas with gene name (not model).
  if(verbose)
    cat("Renaming annotation fasta headers\n")
  tmp <- lapply(subdirs[c("transcript", "cds")], function(x){
    parse_fastaHeader(fasta.dir = x,
                      is.peptide = F,
                      verbose = F,
                      parse_fastaHeader.FUN = parse_fastaHeader.FUN)
  })
  tmp <- parse_fastaHeader(fasta.dir = subdirs$peptide,
                           is.peptide = T,
                           verbose = F,
                           parse_fastaHeader.FUN = parse_fastaHeader.FUN)
  if(verbose)
    cat("Done!\n")
}
