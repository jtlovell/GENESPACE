#' @title pipe_orphanBlast
#'
#' @description
#' \code{pipe_orphanBlast} pipe_syntenicBlocks
#'
#' @param genomeIDs Character vector of genome IDs to consider
#' @param dir.list Directory list, produced by `check_environment`.
#' @param orphan.metadata Directory containing peptide fasta files
#' @param n.cores
#' @param verbose Logical, should updates be printed?
#' @param ... Not in use yet.
#' @details More here
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @import Biostrings
#' @export
pipe_orphanBlast <- function(dir.list,
                             genomeIDs,
                             orphan.metadata,
                             n.cores = 6,
                             verbose = T){

  ########################################################
  ########################################################
  split_pepByOrphan <- function(orphan.md,
                                pep.fastas,
                                n.cores){
    spl <- split.data.table(orphan.md, "unique.block")

    pep.list <- mclapply(spl, mc.cores = n.cores, function(x){
      id <- x$gene2map
      return(pep.fastas[id])
    })
    return(pep.list)
  }

  ########################################################
  ########################################################
  load.annotations <- function(genomeIDs,
                               cds.dir,
                               peptide.dir,
                               assembly.dir){

    cds.fastas <- do.call(c, lapply(genomeIDs, function(i)
      readDNAStringSet(file.path(cds.dir,
                                 paste0(i, ".fa")))))

    pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
      readAAStringSet(file.path(peptide.dir,
                                paste0(i, ".fa")))))

    fais <- rbindlist(lapply(genomeIDs, function(i){
      tmp <- read.delim(file.path(assembly.dir,
                                  paste0(i, ".fa.fai")),
                        header = F,
                        stringsAsFactors = F,
                        col.names = c("chr",  "chr.length",
                                      "v1", "v2", "v3"))[, 1:2]
      tmp$genome <- i
      return(tmp)
    }))

    setkey(fais, genome, chr)
    return(list(peptide = pep.fastas,
                cds = cds.fastas,
                chrlen = fais))
  }
  ########################################################
  ########################################################
  split_bedByOrphan <- function(orphan.md){
    md <- orphan.md[,c("chr","start","end","unique.block","genome")]
    setnames(md, 4, "id")
    md <- md[!duplicated(md), ]
    out.bed <- split.data.table(md[ ,1:4, with = F], "id")
    out.genome <- split(md$genome, md$id)

    return(list(bed.list = out.bed,
                genome.list = out.genome))
  }
  ########################################################
  ########################################################
  run_getfasta <- function(bed, bed.file, ass.file, fa.file){
    write.table(bed,
                file = bed.file,
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = F)
    system(paste("bedtools getfasta -fi", ass.file,
                 "-bed", bed.file,
                 "-name -fo",fa.file))
  }
  ########################################################
  ########################################################
  prep_seq4blast <- function(block.dir,
                             pep.list,
                             genome.list,
                             assembly.dir,
                             bed.list,
                             verbose = T){
    block.dirs <- file.path(block.dir, names(pep.list))
    names(block.dirs) <- names(pep.list)
    n <- length(pep.list)
    nb <- ifelse(n < 100, 10,
                 ifelse(n < 1000, 100,
                        ifelse(n < 5000, 200, 1000)))
    blk.md <- rbindlist(lapply(names(pep.list), function(i){
      if (verbose)
        if (which(names(pep.list) == i) %% nb == 0)
          cat(which(names(pep.list) == i), "/",
              length(pep.list), "\n\t")

      if (dir.exists(block.dirs[[i]]))
        unlink(block.dirs[[i]], recursive = T)
      dir.create(block.dirs[[i]])

      bed.file <- file.path(block.dirs[[i]], "full.block.bed")
      fa.file <- file.path(block.dirs[[i]], "full.block.fa")
      pep.file <- file.path(block.dirs[[i]], "orphan.pep.fa")
      ass.file <- file.path(assembly.dir,
                            paste0(genome.list[[i]], ".fa"))

      writeXStringSet(pep.list[[i]], filepath = pep.file)
      run_getfasta(bed = bed.list[[i]],
                   bed.file = bed.file,
                   ass.file = ass.file,
                   fa.file = fa.file)

      blk.metadata <- data.frame(unique = i,
                                 genome = genome.list[[i]],
                                 path = block.dirs[[i]],
                                 fa.file = fa.file,
                                 pep.file = pep.file,
                                 bed.file = bed.file,
                                 assem.file = ass.file)
      return(blk.metadata)
    }))
    return(blk.md)
  }
  ########################################################
  ########################################################
  run_diamondBlastx <- function(db.file,
                                pep.fa,
                                fa.file,
                                blast.file,
                                max.target.seqs = 10000,
                                min.score = 20,
                                diamond.blastx.param = "--quiet"){
    system(paste("diamond makedb --quiet",
                 "--in", pep.fa,
                 "-d", db.file))
    system(paste("diamond blastx",
                 "--max-target-seqs", max.target.seqs,
                 "--min-score", min.score,
                 diamond.blastx.param,
                 "-d", db.file,
                 "-q", fa.file,
                 "-o", blast.file))
  }
  ########################################################
  ########################################################
  pipe_diamondBlastx <- function(blast.inputs,
                                 verbose = T){
    n <- nrow(blast.inputs)
    nb <- ifelse(n < 100,5,
                 ifelse(n < 1000,50,
                        ifelse(n < 5000, 100, 500)))

    out <- lapply(1:nrow(blast.inputs),  function(i){
      if (verbose)
        if (i %% nb == 0)
          cat(i,"/",nrow(blast.inputs),"\n\t")
      bl <- data.frame(blast.inputs[i,])
      db.file <- file.path(bl$path, "blast.db")
      blast.file <- file.path(bl$path, "out.blast")
      with(bl,
           run_diamondBlastx(db.file = db.file,
                             pep.fa = pep.file,
                             fa.file = fa.file,
                             blast.file = blast.file,
                             max.target.seqs = 10000,
                             min.score = 20,
                             diamond.blastx.param = "--quiet"))
      return(blast.file)
      names(out) <- blast.inputs$unique
    })
    return(out)
  }
  ########################################################
  ########################################################
  read_blastx <- function(blast.files,
                          n.cores = 1){
    out <- mclapply(blast.files, mc.cores = n.cores, function(x){
      rl <- length(readLines(x))
      if (rl < 1) {
        bl <- NA
      }else{
        bl <- fread(x)
      }
      return(bl)
    })
    return(out)
  }
  ########################################################
  ########################################################
  cds.dir <- dir.list$cds
  peptide.dir <- dir.list$peptide
  assembly.dir <- dir.list$assembly
  block.dir <- dir.list$block

  ########################################################
  if (verbose)
    cat("Importing fasta indices, and CDS / peptide sequences ... ")
  annot <- load.annotations(genomeIDs = genomeIDs,
                            cds.dir = cds.dir,
                            peptide.dir = peptide.dir,
                            assembly.dir = assembly.dir)
  cds.fastas <- annot$cds
  pep.fastas <- annot$peptide
  fais <- annot$chrlen
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Splitting orphan peptide sequences by block ... ")
  pep.list <- split_pepByOrphan(orphan.md = orphan.metadata,
                                pep.fastas = pep.fastas,
                                n.cores = n.cores)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Converting blocks to bed format ... ")
  tmp.list <- split_bedByOrphan(orphan.md = orphan.metadata)
  bed.list <- tmp.list$bed.list
  genome.list <- tmp.list$genome.list
  if (verbose)
    cat("Done!\n")


  ########################################################
  if (verbose)
    cat("Writing block directories, peptide fastas and block sequences ... \n\tCompleted\n\t")
  blast.inputs <- prep_seq4blast(block.dir = block.dir,
                                 pep.list = pep.list,
                                 genome.list = genome.list,
                                 assembly.dir = assembly.dir,
                                 bed.list = bed.list,
                                 verbose = T)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Running diamond blastx for", nrow(blast.inputs),
        "blocks ... Completed:\n\t")
  blast.loc <- pipe_diamondBlastx(blast.inputs)
  if (verbose)
    cat("Done!\n")

  ########################################################
  if (verbose)
    cat("Importing blast hits ... \n\t")
  blast <- read_blastx(blast.files = blast.loc, n.cores = 6)
  if (verbose)
    cat("Done!\n")

  return(list(blast.results = blast,
              blast.paths = blast.loc,
              blast.inputs = blast.inputs,
              orphan.metadata = orphan.metadata,
              bed.list = bed.list,
              genome.list = genome.list))
}
