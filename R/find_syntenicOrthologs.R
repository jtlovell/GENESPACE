#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
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
#' @param species.mappings The 'species.mappings' data.table or data.frame
#' from form_syntenicBlocks, giving a dictionary between pairwise
#' blast genome IDs from orthofinder, the blast file locations and
#' the genomeIDs.
#' @param orthogroups The orthogroups list.
#' @param plotit Logical, should plots be made? Will not work with
#' n.core > 1.
#' @param rank.buffer The buffer, in gene rank order.
#' @param n.cores Number of parallel processes to run, when possible
#' @param min.block.size Numeric of length 1, specifying the minimum
#' block size in the cleaning procedure at the end.
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
#' @importFrom compiler cmpfun
#' @importFrom parallel mclapply
#' @export
find_syntenicOrthogs <- function(map,
                                 dir.list,
                                 gff,
                                 gene.index,
                                 species.mappings,
                                 n.cores = 1,
                                 orthogroups,
                                 plotit = n.cores == 1,
                                 rank.buffer = 250,
                                 verbose = T,
                                 min.block.size = 5,
                                 ...){

  #######################################################
  #######################################################
  add_ofnum2gff <- function(gff,
                            gene.index){
    setkey(gff, id)
    setkey(gene.index, id)
    m <- merge(gene.index, gff)
    return(m)
  }

  ########################################################
  ########################################################
  add_ofnum2gff <- cmpfun(add_ofnum2gff)
  ########################################################
  ########################################################

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
  all.blast <- import_ofBlast(species.mappings = species.mappings,
                              genomeIDs = genomeIDs,
                              orthogroups = orthogroups,
                              gene.index = gene.index,
                              gff = gff,
                              verbose = T,
                              only.orthologs = F)
  #######################################################
  if (verbose)
    cat("Adding orthofinder gene IDs to gff data.table ... ")
  gff.names <- colnames(gff)
  gff = add_ofnum2gff(gff,
                      gene.index = gene.index)

  if (verbose)
    cat("Done!\n")
  #######################################################
  if (verbose)
    cat("Culling blast results to regions near blocks ...\n")
  all.syn <- cull_syntenicBlast(gff = gff,
                                map = map,
                                blast = all.blast,
                                plotit = plotit,
                                verbose = verbose,
                                n.cores = n.cores,
                                rank.buffer = rank.buffer)
  if (verbose)
    cat("Done!\n")
  gi <- as.character(gene.index$gene.num)
  names(gi) <- gene.index$id
  ids2keep = with(all.syn,
                  data.table(gn1 = gi[c(id1,id2)],
                             gn2 = gi[c(id2,id1)],
                             stringsAsFactors = F))
  setkey(ids2keep, gn1, gn2)
  #######################################################

  #######################################################
  if (verbose)
    cat("Writing results to file ...\n")
  sm <- res.all$init.results$ortho.info$species.mappings
  sm$tmp.filename <- file.path(dir.list$tmp, basename(sm$filename))

  test <- lapply(1:nrow(sm), function(i){
    out.file = sm$tmp.filename[i]
    g1 = sm$genome1[i]
    g2 = sm$genome2[i]
    is1 <- g1 == sm$ref[i]
    blast.in <- all.blast[with(all.blast,
                               (genome1 == g1 &
                                  genome2 == g2) |
                                 (genome1 == g2 &
                                    genome2 == g1)),]
    if (verbose)
      cat("\t",g1,"-->",g2,"... total/syntenic hits =", nrow(blast.in),"/ ")
    blast.tmp <- blast.in[,c("gn1", "gn2", "perc.iden",
                             "align.length", "n.mismatch",
                             "n.gapOpen", "q.start", "q.end",
                             "s.start", "s.end",
                             "eval", "score")]
    if (!is1){
      setnames(blast.tmp,c("gn1","gn2"),c("gn2","gn1"))
      blast.tmp <- data.table(blast.tmp[,c("gn1", "gn2", "perc.iden",
                                           "align.length", "n.mismatch",
                                           "n.gapOpen", "q.start", "q.end",
                                           "s.start", "s.end",
                                           "eval", "score")])
    }
    setkey(blast.tmp, gn1, gn2)
    blast2keep = merge(ids2keep, blast.tmp)
    if(verbose)
      cat(nrow(blast2keep),"\n")
    write.table(blast2keep,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = out.file)
  })

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
  gff <- data.table(gff[,gff.names,with = F])
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

  syn.map <- clean_blocks(map = all.blast,
                          radius = rank.buffer,
                          n.mappings = min.block.size,
                          clean.columns = F)
  syn.map = with(syn.map,
                 merge_blocks(map = map,
                              blk = block,
                              buffer = min.block.size*4,
                              verbose = F,
                              clean.columns = F))
  #######################################################

  #######################################################
  return(list(map = syn.map$map,
              block = syn.map$block,
              blast = all.blast,
              of.results = of.blast))
}
