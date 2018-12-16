#' @title Rerun orthofinder utility functions
#' @name rerun_utilities
#' @aliases subset_blast
#' @aliases remake_blast
#' @aliases pull_gff
#' @aliases read_allBlast
#' @aliases make_mapFromOGs
#' @aliases make_blockMetadata
#' @aliases find_orphans
#' @aliases load.annotations
#'
#' @description
#' \code{rerun_utilities} Several utilities functions meant for internal calls in rerun_orthofinder
#' @name rerun_utilities
#' @param blast.file file path for blast output
#' @param genenum.list list of gene numbers / ids
#' @param blast.dir path for blast directory
#' @param cull.blast.dir path for culled blast directory
#' @param gff concatenated, simplified gff object
#' @param blk.line one line from a blk object
#' @param gene.index gene number / name data.table
#' @param gff.wNum gff with number and block.id added in
#' @param cds.md metadata of cds, including sequence lengths
#' @param blk blk object
#' @param gffog gff object with orthogroup info
#' @param spl.gff split gff object
#' @param blk.md2 simplified block metadata object
#' @param n.cores number of parallel processes to run
#' @param verbose should updates be printed?
#' @param ... not currently in use
#'

#' @title Subset blast results within blocks
#' @description
#' \code{subset_blast} Subset blast results.
#' @rdname blk_utilities
#' @import data.table
#' @export
subset_blast <- function(blast.file,
                         genenum.list,
                         n.cores,
                         verbose){
  if (verbose)
    cat("Subsetting:",basename(blast.file))
  f <- fread(blast.file,
             header = F,
             stringsAsFactors = F,
             check.names = F)

  if (verbose)
    cat(paste0(" (initial hits = ",nrow(f),") "))

  setkey(f, V1, V2)

  fl <- rbindlist(mclapply(genenum.list, mc.cores = n.cores, function(x){
    return(f[f$V1  %in% x & f$V2 %in% x,])
  }))
  fl <- fl[!duplicated(fl),]

  if (verbose)
    cat("to", nrow(fl),
        "hits in blocks\n")

  write.table(fl,
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F,
              file = blast.file)
}



#' @title Subset blasts in new directory
#' @description
#' \code{remake_blast} Subset blasts in new directory
#' @rdname blk_utilities
#' @import data.table
#' @export
remake_blast <- function(blast.dir,
                         cull.blast.dir,
                         genenum.list,
                         n.cores,
                         verbose = T){
  if (verbose)
    cat("Copying blast results to",cull.blast.dir,"... ")
  if (dir.exists(cull.blast.dir))
    unlink(cull.blast.dir)

  if (!dir.exists(cull.blast.dir))
    dir.create(cull.blast.dir)

  blast.dir <- cull.blast.dir
  tmp.dir <- blast.dir

  blast.loc <- dirname(list.files(tmp.dir,
                                  pattern = "SequenceIDs",
                                  recursive = T,
                                  full.names = T)[1])
  ortho.loc <- dirname(list.files(tmp.dir,
                                  pattern = "Orthogroups.txt",
                                  recursive = T,
                                  full.names = T)[1])


  blast.files <- list.files(blast.loc,
                            pattern = "Blast*",
                            full.names = T)
  fa.files <- list.files(blast.loc,
                         pattern = "Species*",
                         full.names = T)
  fa.files <- fa.files[grep(".fa$", fa.files)]

  dmnd.files <- list.files(blast.loc,
                           pattern = "diamondDBSpecies*",
                           full.names = T)
  og.files <- file.path(ortho.loc,
                        "Orthogroups.txt")

  sp.id.files <- file.path(blast.loc,"SpeciesIDs.txt")
  seq.id.files <- file.path(blast.loc,"SequenceIDs.txt")

  files <- c(blast.files,
             fa.files,
             dmnd.files,
             og.files,
             sp.id.files,
             seq.id.files)

  nu <- file.copy(files,
                  blast.dir)

  if (verbose)
    cat("Done!\n")

  nu <- file.remove(file.path(cull.blast.dir,
                              "Orthogroups.txt"))
  blast.files <- list.files(cull.blast.dir,
                           full.names = T,
                           pattern = "Blast")
  blast.files<-blast.files[grep(".txt$",blast.files)]

  for (i in blast.files){
    subset_blast(blast.file = i,
                 genenum.list = genenum.list,
                 n.cores = n.cores,
                 verbose = verbose)
  }
}

#' @title merge gff with gene numbers
#' @description
#' \code{pull_gff} merge gff with gene numbers
#' @rdname blk_utilities
#' @import data.table
#' @export
pull_gff <- function(gff,
                    blk.line,
                    gene.index){
  x <- blk.line
  setkey(gene.index, id)

  g1 <- gff[with(gff, genome == x$genome1 &
                  chr == x$chr1 &
                  start <= x$end1 &
                  end >= x$start1),]
  g1$block.id <- x$block.id
  setkey(g1, id)

  g2 <- gff[with(gff, genome == x$genome2 &
                  chr == x$chr2 &
                  start <= x$end2 &
                  end >= x$start2),]
  g2$block.id <- x$block.id
  setkey(g2, id)

  g1 <- merge(gene.index, g1)
  g2 <- merge(gene.index, g2)

  return(list(g1, g2))
}

#' @title read blast files
#' @description
#' \code{read_allBlast} read blast files
#' @rdname blk_utilities
#' @import data.table
#' @export
read_allBlast <- function(blast.dir,
                          verbose = T){

  if (verbose)
    cat("Importing orthofinder results:\n\t")
  gz <- list.files(blast.dir,
                   pattern = ".gz$")
  if (length(gz) > 0) {
    if (verbose)
      cat("Decompressing blast results\n")
    system(paste("gunzip -f",
                 file.path(blast.dir,
                           "*.gz")))
  }

  blast.files <- list.files(blast.dir,
                            full.names = T,
                            pattern = "^Blast")
  blast.files<-blast.files[grep(".txt$",blast.files)]
  print(blast.files)

  out <- rbindlist(lapply(blast.files, function(x)
    fread(x, header = F,
          stringsAsFactors = F,
          check.names = F,
          col.names = c("gene.num1", "gene.num2",
                        "perc.iden", "align.length",
                        "n.mismatch", "n.gapOpen",
                        "q.start", "q.end",
                        "s.start", "s.end",
                        "eval", "score"))))
  return(out)
}

#' @title buid gff networks from orthogroups
#' @description
#' \code{make_mapFromOGs} buid gff networks from orthogroups
#' @rdname blk_utilities
#' @import data.table
#' @export
make_mapFromOGs <- function(gff.wNum,
                            cull.blast.dir){

  og <- readLines(file.path(cull.blast.dir,
                            "Orthogroups.txt"))

  og <- lapply(og, function(x)
    strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x)
    x[1])
  names(og) <- ons

  og <- lapply(og, function(x)
    x[-1])
  og.name <- names(og)
  og.length <- sapply(og, length)

  od <- data.table(og.id = rep(og.name, og.length),
                  id = unlist(og),
                  stringsAsFactors = F)
  setkey(od, id)

  gffn <- rbindlist(unlist(gff.wNum,
                           recursive = F))
  setkey(gffn,id)
  gffn <- merge(od, gffn)

  gffn[, complete := all(duplicated(block.id) |
                           duplicated(block.id, fromLast = T)),
       by = og.id]
  setkey(gffn, id)
  return(gffn)
}



#' @title combine block information with gff of orthogroups
#' @description
#' \code{make_blockMetadata} combine block information with gff of orthogroups
#' @rdname blk_utilities
#' @import data.table
#' @export
make_blockMetadata <- function(cds.md,
                               gffog,
                               blk){
  gffog <- merge(cds.md, gffog)

  gffog.incomplete <- gffog[!gffog$complete,]
  gffog.incomplete$unique.block <- with(gffog.incomplete,
                                        paste0(genome,"_", block.id))

  gene.block <- gffog[,list(block = unique(block.id)),
                      by = id]
  gene.block <- split(gene.block$block, gene.block$id)

  blk.list <- split.data.table(blk, "block.id")
  spl.gff <- split.data.table(gffog.incomplete, "og.id")

  blk.md <- blk[,list(genome = c(genome1, genome2),
                      chr = c(chr1,chr2),
                      start = c(start1,start2),
                      end = c(end1,end2)),
                by = list(block.id)]

  blk.md$unique.block <- with(blk.md, paste0(genome,"_", block.id))
  setkey(blk.md, unique.block, block.id, genome)

  blk.md2 <- blk.md[,c("block.id","genome","unique.block")]
  setkey(blk.md2, unique.block)

  return(list(spl.incompleteGff = spl.gff,
              blk.metadata = blk.md,
              simple.blk.metadata = blk.md2))

}

#' @title find orphan genes in blocks
#' @description
#' \code{find_orphans} find orphan genes in blocks
#' @rdname blk_utilities
#' @import data.table
#' @export
find_orphans <- function(spl.gff,
                         blk.md2,
                         n.cores){
  y <- rbindlist(mclapply(spl.gff, mc.cores = n.cores, function(x){
    xmd <- blk.md2[blk.md2$block.id %in% unique(x$block.id), ]
    xmd <- xmd[!xmd$unique.block %in% x$unique.block, ]
    xmd$gene2map <- x$id[which.max(x$length)]
    xmd$og.id <- x$og.id[1]
    return(xmd)
  }))

  setkey(y, unique.block, block.id, genome)

  return(y)
}

#' @title load annotation to memory
#' @description
#' \code{load.annotations} load annotation to memory
#' @rdname blk_utilities
#' @import data.table
#' @export
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


convert_blast2blk <- function(){

}
