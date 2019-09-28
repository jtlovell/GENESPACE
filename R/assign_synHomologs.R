#' @title Assign syntenic homologs
#'
#' @description
#' \code{assign_synHomologs} Assign syntenic homologs as paralogs, orthologs,
#' or just homologs.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param syn.blast data.table, optional, containing the
#' syntenic blast hits to feed to orthofinder.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param rank.buffer numeric, the radius to search within for syntenic
#' mappings.
#' @param min.homolog.score numeric, minimum blast bit score to be
#' included
#' @param n.cores integer length 1, the number of parallel processes
#' to run.
#' @param min.perc.iden numeric, minimum blast proprotion of sequence
#' identity to be included
#' @param min.prop.of.best numeric, minimum proportion of top score for
#' each gene to be included.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use

#' @details ...
#'
#' @return A blast data.table, with an additional column indicating
#' whether each is is homologous, paralogous or orthologous
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom igraph clusters graph_from_data_frame
#' @export
assign_synHomologs <- function(gff,
                               genomeIDs,
                               map,
                               dir.list,
                               rank.buffer,
                               min.propBestScore4homolog = .5,
                               min.score4homolog = 50,
                               verbose = T,
                               quiet.orthofinder = F,
                               n.cores = 1){

  ####################################################################
  # - Get the blast hits to feed to orthofinder
  if(verbose)
    cat("Pulling syntenic blast hits (this might take a while) ... ")
  syn.blast <- extend_blocks(
    gff = gff,
    genomeIDs = genomeIDs,
    map = mirror_map(map),
    dir.list = dir.list,
    use.score.cull.blast = F,
    rank.buffer = rank.buffer,
    return.reg.only = T,
    verbose = F)
  syn.blast[, genome.num1 := NULL]
  syn.blast[, genome.num2 := NULL]


  ####################################################################
  # - Move over peptides, make gene/genomeIDs, etc.
  if(verbose)
    cat("Done!\nMaking new orthofinder database ... ")
  make_newOFdb(
    tmp.dir = dir.list$tmp,
    peptide.dir = dir.list$peptide,
    genomeIDs = genomeIDs,
    verbose = F,
    n.cores = 1,
    output.dir = dir.list$syn.blast)

  ####################################################################
  # - Read in the new ids, merge with blast and write to file
  if(verbose)
    cat("Done!\nReformatting blast results ... ")
  si <- read_speciesIDs(
    of.dir = dir.list$syn.blast,
    genomeIDs = genomeIDs)

  gi <- read_geneIDs(
    of.dir = dir.list$syn.blast,
    gff = gff,
    species.num.id = si)

  blast.in <- merge_gffWithBlast(
    gff = gi,
    blast = syn.blast,
    include.gene.num = T,
    include.genome.num = T,
    mirror = T)

  if(verbose)
    cat("Done!\nWriting blast results to the orthofinder directory ... ")
  write_ofBlast2file(
    blast = blast.in,
    output.dir = dir.list$syn.blast)

  ####################################################################
  # - Run orthofinder on syntenic blast hits
  if(verbose)
    cat("Done!\nRunning orthofinder (this might take a while) ... ")
  com <- paste("orthofinder",
               "-b", dir.list$syn.blast,
               "-a", n.cores)
  if(quiet.orthofinder)
    com <- paste(com, "1>/dev/null 2>&1")
  system(com)

  ####################################################################
  # - read in ids and genomes of genes in orthogroups
  if(verbose)
    cat("Done!\nParsing orthogroups ... ")
  orthogroup.dt <- read_ogs(
    of.dir = dir.list$syn.blast,
    gff = gff)
  orthogroup.dt[,is.ortholog := id1 == id2 & genome1 == genome2]

  if(verbose)
    cat("Done!\nExtracting orthologs ... ")
  ortholog.dt <- pull_orthologs(
    of.dir = dir.list$syn.blast)

  if(verbose)
    cat("Done!\nMerging blast and qualifying orthologs, paralogs and syntenic homologs ... ")
  ogs <- rbind(ortholog.dt,
               orthogroup.dt)
  ogs <- ogs[!duplicated(ogs[,c("genome1","genome2","id1","id2")]),]

  blast.out <- merge(
    blast.in,
    ogs,
    by = c("genome1","genome2","id1","id2"),
    all.x = T)

  blast.out[,best1 := max(score),
            by = list(genome1, genome2, id1)]
  blast.out[,best2 := max(score),
            by = list(genome2, genome1, id2)]

  blast.homo <- subset(
    blast.out,
    (score >= best1 * min.propBestScore4homolog &
       score >= best2 * min.propBestScore4homolog &
       score >= min.score4homolog) |
      !is.na(is.ortholog))
  blast.homo[,best1 := NULL]
  blast.homo[,best2 := NULL]

  blast.homo[,hit.type := ifelse(is.na(is.ortholog), "syntenic.homolog",
                                 ifelse(is.ortholog, "ortholog","paralog"))]
  if(verbose)
    with(blast.homo,
         cat("Done!\nReturning a blast data set with:\n\tn.orthologs =",
             sum(hit.type == "ortholog"),"\n\tn.paralogs =",
             sum(hit.type == "paralog"),"\n\tn.other syntenic homologs =",
             sum(hit.type == "syntenic.homolog"),"\nDone!\n"))
  return(blast.homo)
}
