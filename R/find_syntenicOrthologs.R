#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param blk The block data.frame or data.table
#' @param map The map data.frame or data.table
#' @param dir.list The directory list produced by check_environment
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param gene.index orthofinder geneID data.table or data.frame
#' giving a dictionary between the 'id' column in the gff object
#' and the 'gene.num' numeric geneIDs in the orthofinder-formatted
#' blast files.
#' @param species.index The 'species.mappings' data.table or data.frame
#' from form_syntenicBlocks, giving a dictionary between pairwise
#' blast genome IDs from orthofinder, the blast file locations and
#' the genomeIDs.
#' @param n.cores Number of parallel processes to run, when possible
#' @param min.block.size The minimum block size to retain.
#' @param verbose Logical, should updates be printed
#' @param ... Not currently in use
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
find_syntenicOrthogs <- function(map,
                                 blk,
                                 dir.list,
                                 gff,
                                 gene.index,
                                 species.index,
                                 n.cores = 1,
                                 min.block.size = 5,
                                 blk.bp.buffer = 1e5,
                                 verbose = T){

  #######################################################
  #######################################################
  add_ofnum2gff <- function(gff,
                            gene.index){
    setkey(gff, id)
    setkey(gene.index, id)
    m <- merge(gene.index, gff)
    return(m)
  }
  #######################################################
  #######################################################
  pull_idByBlk <- function(gff,
                           blk,
                           n.cores,
                           buffer){
    out <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
      pull_gff(gff = gff,
               blk.line = blk[i,],
               buffer = buffer))
    names(out) <- blk$block.id
    return(out)
  }
  #######################################################
  #######################################################
  pull_gff <- function(gff,
                       blk.line,
                       buffer){
    x <- blk.line

    g1 <- gff[with(gff, genome == x$genome1 &
                     chr == x$chr1 &
                     start <= x$end1+buffer &
                     end >= x$start1-buffer),]

    g2 <- gff[with(gff, genome == x$genome2 &
                     chr == x$chr2 &
                     start <= x$end2+buffer &
                     end >= x$start2-buffer),]

    return(rbind(g1,g2))
  }
  #######################################################
  #######################################################
  parse_blast <- function(species.index,
                          blk,
                          gff,
                          gene.index,
                          dir.list,
                          verbose,
                          genomeIDs,
                          n.cores,
                          buffer){

    cmb <- combn(genomeIDs, 2)
    for (i in genomeIDs)
      cmb <- cbind(cmb, rep(i,2))

    out <- apply(cmb, 2, function(gs){

      g1 <- gs[1]
      g2 <- gs[2]

      if (verbose)
        cat(g1,"-->",
            g2,"\n")

      gff.tmp <- gff[gff$genome %in% gs,]
      blk.tmp <- blk[(blk$genome1 == g1 &
                       blk$genome2 == g2) |
                       (blk$genome1 == g2 &
                          blk$genome2 == g1), ]

      gene.list <- pull_idByBlk(gff.tmp,
                                blk = blk.tmp,
                                n.cores = n.cores,
                                buffer = buffer)
      id.list <- lapply(gene.list, function(x) x$id)
      gene.list <- lapply(gene.list, function(x) x$gene.num)

      blast.files <- species.index$filename[with(species.index,
                                                 (genome1 == g1 &
                                                    genome2 == g2) |
                                                   (genome1 == g2 &
                                                      genome2 == g1))]
      blast.files <- file.path(dir.list$tmp, basename(blast.files))

      for (i in blast.files) {
        subset_blast(blast.file = i,
                     blk = blk.tmp,
                     verbose = T)
      }
      return(list(gene.list = gene.list,
                  id.list = id.list))
    })

    return(out)
  }

  #######################################################
  #######################################################
  subset_blast <- function(blast.file,
                           blk,
                           verbose){
    if (verbose)
      cat("\tSubsetting:",basename(blast.file))
    f <- fread(blast.file,
               header = F,
               stringsAsFactors = F,
               check.names = F)

    if (verbose)
      cat(paste0(" (initial hits = ", nrow(f),") "))


    f$ns = f$V12 * (-1)
    setkey(f, ns)
    fo <- f[!duplicated(f[,f[,1:2,with = F]]),]
    fo$ns <- NULL

    setkey(fo, V1, V2)

    blast = get_nearestDistance(map = cleaned2$map,
                                blast = res.all$init.results$orthogroup.blast)
    with(blast, table(distance.1 < 1e6 & distance.2 < 1e6))
    blast2keep= blast[with(blast, distance.1 < 1e6 & distance.2 < 1e6),]

    fl$ns = fl$V12 * (-1)
    setkey(fl, ns)
    fl <- fl[!duplicated(fl[,fl[,1:2,with = F]]),]
    fl$ns <- NULL
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
  ########################################################
  ########################################################
  convert_blast2blk = function(id.list,
                               blast,
                               n.cores){

    fl <- rbindlist(mclapply(names(id.list), mc.cores = n.cores, function(i){
      x = id.list[[i]]
      tmp = blast[blast$id1 %in% x & blast$id2 %in% x,]
      tmp$block.id = i
      return(tmp)
    }))
    return(make_blocks(fl))
  }
  ########################################################
  ########################################################

  #######################################################
  if (verbose)
    cat("Adding orthofinder gene IDs to gff data.table ... ")

  gff = add_ofnum2gff(gff,
                      gene.index = gene.index)
  gff = gff[gff$genome %in% unique(c(map$genome1, map$genome2)),]

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Cleaning out tmp and culled blast directories ... ")

  unlink(dir.list$tmp, recursive = T)
  dir.create(dir.list$tmp)

  unlink(dir.list$cull.blast, recursive = T)
  dir.create(dir.list$cull.blast)

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Copying blast results to tmp directory ... ")

  files <- list.files(dir.list$blast,
                      full.names = T)
  nu <- file.copy(files,
                  dir.list$tmp)
  nu <- file.remove(file.path(dir.list$tmp,
                              "Orthogroups.txt"))
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Parsing blast results to hits within blocks ...\n")

  gene.lists = parse_blast(blk = blk,
                           genomeIDs = genomeIDs,
                           dir.list = dir.list,
                           gene.index =  gene.index,
                           gff = gff,
                           n.cores = 6,
                           species.index = species.index,
                           verbose = T,
                           buffer = blk.bp.buffer)

  gene.lists = c(gene.lists, gene.lists)
  id.list = lapply(gene.lists, function(x) x$id.list)
  names(id.list)<-apply(combn(genomeIDs, 2),2, function(x) paste(x, collapse = "_"))
  gs.idlist = lapply(names(id.list), function(x) strsplit(x,"_")[[1]])

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Re-running orthofinder on culled blast hits ...\n")

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = dir.list$tmp,
    blast.dir = dir.list$cull.blast,
    og.threads = 6,
    og.silent = F,
    verbose = T)

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  gff$gene.num <- NULL
  of.blast <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = dir.list$cull.blast,
    verbose = T)
  #######################################################

  #######################################################
  all.blast <- import_ofBlast(
    species.mappings = of.blast$species.mappings,
    genomeIDs = genomeIDs,
    orthogroups = of.blast$orthogroups,
    gff = gff,
    gene.index = of.blast$gene.index,
    verbose = T)

  all.blast$unique = with(all.blast, paste0(genome1, "_", genome2))
  spl = split(all.blast, "unique")
  #######################################################

  #######################################################
  if (verbose)
    cat("Re-building map and block data.tables ... \n")
  block.list = lapply(names(spl), function(x){
    g = strsplit(x,"_")[[1]]
    g1 = g[1]
    g2 = g[2]
    if(g1 != g2){
      if(verbose)
        cat("\t",g1,"-->",g2,"... ")
      wh = which(sapply(gs.idlist, function(y) y[1] %in% g & y[2] %in% g))
      tmp = convert_blast2blk(id.list = id.list[[wh]],
                              blast = spl[[x]],
                              n.cores = 6)
      if (verbose)
        cat("Done!\n")

      return(tmp)
    }
  })
  #######################################################

  #######################################################
  if (verbose)
    cat("Combining block file across all pairwise comparisons ...")
  map.comb = rbindlist(lapply(block.list, function(x) x$map))
  blk.out = make_blocks(map.comb,
                        rename.blocks = T)
  if (verbose)
    cat("\n\tReturning a dataset of",
        nrow(blk.out$block),"blocks containing",
        nrow(blk.out$map),"hits\n")
  #######################################################

  #######################################################
  if(verbose)
    cat("Dropping small blocks ... ")
  syn.out = drop_smallBlocks(map = blk.out$map,
                             blk = blk.out$block,
                             min.block.size = min.block.size)
  if (verbose)
    cat("\n\tReturning a dataset of",
        nrow(syn.out$block),"blocks containing",
        nrow(syn.out$map),"hits\nDone!\n")
  #######################################################

  #######################################################
  return(list(blast = all.blast,
              of.results = of.blast,
              block = syn.out$block,
              map = syn.out$map,
              id.list = id.list))
}
