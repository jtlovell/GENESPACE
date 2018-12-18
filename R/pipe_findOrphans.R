  #' @title Find orphan genes
  #'
  #' @description
  #' \code{pipe_findOrphans} Search orthofinder output for genes withou
  #' hits in one of the blocks in the network.
  #'
  #' @param blk block results data.table
  #' @param gff concatenated gff data.table
  #' @param gene.index geneID gene number mapping
  #' @param dir.list directory list, created by check_environment
  #' @param n.cores numeric, the number of parallel processes to run
  #' @param verbose Should updates be printed?
  #' @param ... Not currently in use
  #' @details ...
  #' @return Nothing.
  #'
  #' @examples
  #' \dontrun{
  #' none yet
  #' }
  #' @import data.table
  #' @export
  pipe_findOrphans = function(blk,
                              gff,
                              gene.index,
                              dir.list,
                              n.cores = 1,
                              verbose = T){


    ########################################################
    ########################################################
    run_pullGff = function(gff, blk, gene.index, n.cores = 1){
      out <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
        pull_gff(gff = gff,
                 blk.line = blk[i,],
                 gene.index = gene.index))
      names(out)<-blk$block.id
      return(out)
    }
    ########################################################
    ########################################################
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
    ########################################################
    ########################################################
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
    ########################################################
    ########################################################
    add_cdsLength2gff = function(cds.dir, genomeIDs, gff){
      cds.fastas <- do.call(c,lapply(genomeIDs, function(i)
        readDNAStringSet(file.path(cds.dir,
                                   paste0(i,".fa")))))
      cds.md <- data.table(id = names(cds.fastas),
                           length = width(cds.fastas),
                           stringsAsFactors = F)
      setkey(cds.md, id)
      setkey(gff, id)
      out <- merge(gff, cds.md)
      return(out)
    }
    ########################################################
    ########################################################
    make_incompleteGfflist = function(gffog){
      tmp <- gffog[!gffog$complete,]
      tmp$unique.block <- with(tmp,
                               paste0(genome,"_", block.id))
      spl.gff <- split.data.table(tmp, "og.id")
      return(spl.gff)
    }
    ########################################################
    ########################################################
    make_blkMetadata = function(blk){
      tmp <- blk[,list(genome = c(genome1, genome2),
                       chr = c(chr1,chr2),
                       start = c(start1,start2),
                       end = c(end1,end2)),
                 by = list(block.id)]
      tmp$unique.block <- with(tmp, paste0(genome,"_", block.id))
      setkey(tmp, unique.block, block.id, genome)
      return(tmp)
    }
    ########################################################
    ########################################################
    find_orphan <- function(gff.incomplete,blk.md){
      ublk = unique(gff.incomplete$block.id)
      uublk = unique(gff.incomplete$unique.block)

      xmd <- blk.md[blk.md$block.id %in% ublk, ]
      xmd <- xmd[!xmd$unique.block %in% uublk, ]
      xmd$gene2map <- gff.incomplete$id[which.max(gff.incomplete$length)]
      xmd$og.id <- gff.incomplete$og.id[1]
      return(xmd)
    }
    ########################################################
    ########################################################
    run_findOrphans <- function(spl.gff,
                                blk.md,
                                n.cores = 1){
      y <- rbindlist(mclapply(spl.gff, mc.cores = n.cores, function(x){
        find_orphan(gff.incomplete = x, blk.md = blk.md)
      }))

      setkey(y, unique.block, block.id, genome)
      return(y)
    }


    ########################################################
    ########################################################

    ########################################################
    ########################################################

    ########################################################
    ########################################################

    if(verbose)
      cat("Pulling gff file for each block ... ")
    gff.wNum <- run_pullGff(gff = gff,
                            blk = blk,
                            gene.index = gene.index,
                            n.cores = n.cores)
    if(verbose)
      cat("Done!\n")

    ########################################################
    if(verbose)
      cat("Adding orthogroup info into gff ... ")
    gffog <- make_mapFromOGs(gff.wNum = gff.wNum,
                             cull.blast.dir = dir.list$cull.blast)
    if(verbose)
      cat("Done!\n")

    ########################################################
    if(verbose)
      cat("Adding CDS sequence length into gff ... ")
    gffog <- add_cdsLength2gff(cds.dir = dirs$cds,
                               genomeIDs = genomeIDs,
                               gff = gffog)
    if(verbose)
      cat("Done!\n")

    ########################################################
    if(verbose)
      cat("Finding orthogroups that do not include all possible blocks ... ")
    spl.gff <- make_incompleteGfflist(gffog = gffog)
    blk.md <- make_blkMetadata(blk)
    if(verbose)
      cat("Done!\n")

    ########################################################
    if(verbose)
      cat("Extracting orphan genes in orthogroups ... ")
    orphan.md <- run_findOrphans(spl.gff = spl.gff,
                                 blk.md = blk.md,
                                 n.cores = n.cores)
    if(verbose)
      cat("Done!\n")

    ########################################################
    return(list(orphan.metadata = orphan.md,
                gene.metadata = gffog))
  }
