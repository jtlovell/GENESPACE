#' @title pipe_syntenicBlocks
#'
#' @description
#' \code{pipe_syntenicBlocks} pipe_syntenicBlocks
#'
#' @param genomeIDs Character vector of genome IDs to consider
#' @param rerun.results Directory list, produced by `check_environment`.
#' @param peptide.dir Directory containing peptide fasta files
#' @param gff gff object from synteny analysis
#' @param tmp.dir Directory to write temporary files
#' @param cds.dir Directory containing cds fasta files
#' @param assembly.dir Directory containing assembly fasta files
#' @param n.cores
#' @param min.blast.score
#' @param min.exonerate.score
#' @param other.blastx.param
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
pipe_blast2exonerate <- function(genomeIDs,
                                 rerun.results,
                                 gff,
                                 cds.dir,
                                 peptide.dir,
                                 assembly.dir,
                                 tmp.dir,
                                 n.cores = 1,
                                 mode = "sensitive",
                                 min.blast.score = 50,
                                 min.exonerate.score = 50,
                                 other.blastx.param = NULL,
                                 verbose = T){

  if (verbose)
    cat("Importing fasta indices, and CDS / peptide sequences\n")
  cds.fastas <- do.call(c,lapply(genomeIDs, function(i)
    readDNAStringSet(file.path(cds.dir,
                               paste0(i,".fa")))))

  pep.fastas <- do.call(c,lapply(genomeIDs, function(i)
    readAAStringSet(file.path(peptide.dir,
                              paste0(i,".fa")))))

  fais <- rbindlist(lapply(genomeIDs, function(i){
    tmp <- read.delim(file.path(assembly.dir,
                                paste0(i,".fa.fai")),
                      header = F, stringsAsFactors = F,
                      col.names = c("chr",  "chr.length", "v1", "v2", "v3"))[, 1:2]
    tmp$genome <- i
    return(tmp)
  }))

  setkey(fais, genome, chr)


  rerun.results <- res.all$rerun.results
  uniq <- with(rerun.results$block,
               paste(genome1, genome2, block.id, sep = "_"))
  out <- lapply(1:length(uniq), function(i){

    ###################################
    # -- Find the genes in the blocks
    if (verbose)
      cat("Running block", i, "/",
          nrow(rerun.results$block),"... ")
    blk <- rerun.results$block[i,]
    map <- rerun.results$map
    map <- map[map$block.id == blk$block.id,]
    if (verbose)
      cat("N genes in OGs =",
          length(unique(map$id1)),
          paste0("(", blk$genome1, ")"),
          length(unique(map$id2)),
          paste0("(", blk$genome2, ")\n\t"))

    ###################################
    # -- Make the block-specific working directory
    wd <- file.path(tmp.dir, blk$block.id)
    if (!dir.exists(wd))
      dir.create(wd)

    ###################################
    # -- Find genes in annotation region, but without orthogrous: orphans
    gff1 <- gff[with(gff,
                     genome == blk$genome1 &
                       chr == blk$chr1 &
                       start <= blk$end1 &
                       end >= blk$start1),]

    gff2 <- gff[with(gff,
                     chr == blk$chr2 &
                       start <= blk$end2 &
                       end >= blk$start2),]

    orphan1 <- gff1$id[!gff1$id %in% map$id1]
    orphan2 <- gff2$id[!gff2$id %in% map$id2]
    cat(length(orphan1), "/",
        length(orphan2),"(orphan)")

    ###################################
    # - Pull CDS and peptides from orphans
    cds1 <- cds.fastas[orphan1]
    cds2 <- cds.fastas[orphan2]
    pep1 <- pep.fastas[orphan1]
    pep2 <- pep.fastas[orphan2]

    pep1.file <- file.path(wd, "pep1.fa")
    pep2.file <- file.path(wd, "pep2.fa")
    writeXStringSet(pep1, filepath = pep1.file)
    writeXStringSet(pep2, filepath = pep2.file)

    cds1.file <- file.path(wd, "cds1.fa")
    cds2.file <- file.path(wd, "cds2.fa")
    writeXStringSet(cds1, filepath = cds1.file)
    writeXStringSet(cds2, filepath = cds2.file)

    ###################################
    # - Make region bed files
    bed1 <- with(blk,
                 data.frame(chr1, start1, end1))
    bed1$id <- with(blk,
                    paste0(genome1,".", block.id,
                           ".", chr1, ".", start1))
    bed2 <- with(blk,
                 data.frame(chr2, start2, end2))
    bed2$id <- with(blk,
                    paste0(genome2,".",block.id,
                           ".", chr2, ".",start2))

    wd <- file.path(tmp.dir, blk$block.id)
    if (!dir.exists(wd))
      dir.create(wd)

    bed1.file <- file.path(wd, "bed1.bed")
    bed2.file <- file.path(wd, "bed2.bed")

    ###################################
    # - Extract assembly fasta for region
    ass1.file <- file.path(assembly.dir,
                           paste0(blk$genome1,".fa"))
    ass2.file <- file.path(assembly.dir,
                           paste0(blk$genome2,".fa"))

    fa1.file <- file.path(wd, "fa1.fa")
    fa2.file <- file.path(wd, "fa2.fa")

    run_getfasta(bed = bed1,
                 bed.file = bed1.file,
                 ass.file = ass1.file,
                 fa.file = fa1.file)
    run_getfasta(bed = bed2,
                 bed.file = bed2.file,
                 ass.file = ass2.file,
                 fa.file = fa2.file)

    ###################################
    # - Run blastX
    db1.file <- file.path(wd, "db1")
    db2.file <- file.path(wd, "db2")

    blast1.file <- file.path(wd, "1v2.blast")
    blast2.file <- file.path(wd, "2v1.blast")

    run_diamondBlastx(
      db.file = db1.file,
      pep.fa = pep1.file,
      fa.file = fa2.file,
      blast.file = blast1.file,
      mode = mode,
      n.cores = n.cores,
      min.score = min.blast.score,
      other.blastx.param = other.blastx.param)
    run_diamondBlastx(
      db.file = db2.file,
      pep.fa = pep2.file,
      fa.file = fa1.file,
      blast.file = blast2.file,
      mode = mode,
      n.cores = n.cores,
      min.score = min.blast.score,
      other.blastx.param = other.blastx.param)

    ###################################
    # - Parse blastx results
    rl1 = length(readLines(blast1.file))
    rl2 = length(readLines(blast2.file))
    if(rl1 > 0 | rl2 > 0){
      if(rl1 > 0 & rl2 > 0){
        blast1 <- fread(blast1.file)
        blast2 <- fread(blast2.file)
        bl = rbind(blast1,blast2)
      }else{
        if(rl2 > 0){
          bl = fread(blast2.file)
        }else{
          if(rl1 > 0){
            bl = fread(blast1.file)
          }
        }
      }

      blast <- parse_blastLoc(bl)
      hit.reg <- make_blastRegion(bl = blast)
      if(verbose)
        cat(", ", sum(hit.reg$genome == blk$genome2),
            "/", sum(hit.reg$genome == blk$genome1),
            "(blast) ")
      hits2exonerate <- add_buffer(bl = hit.reg,
                                   fais = fais,
                                   buffer = 1e3)

      ###################################
      # - Run exonerate
      exon.out <- pipe_exonerate(
        hit.reg = hits2exonerate,
        assembly.dir = assembly.dir,
        genomeIDs = genomeIDs,
        cds.fastas = cds.fastas,
        tmp.dir = wd,
        n.cores = n.cores,
        min.score = min.exonerate.score,
        verbose = F)

      exc <- exon.out$gff[complete.cases(exon.out$gff), ]
      exc <- sapply(exc$id, function(x)
        strsplit(x, "_", fixed = T)[[1]][2])

      if (verbose)
        cat(", ", sum(exc == blk$genome2),
            "/", sum(exc == blk$genome1),
            "(exonerate)\n")
    }else{
      exon.out = NULL
      hit.reg = NULL
      if(verbose)
        cat("... no blast hits found\n")
    }
    return(list(exonerate.output = exon.out,
                blast.output = hit.reg))
  })
  names(out) <- uniq
  if (verbose)
    cat("Done!\n")
  return(out)
}
