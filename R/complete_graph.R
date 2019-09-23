#' @title Fill out culled blast hits
#'
#' @description
#' \code{complete_graph} Takes the network of known (potentially syntenic)
#' blast hits assigned to individual orthogroups and extends to all possible
#' combinations.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param ignore.self logical, should the subgraphs be built agnostic to
#' the existence of within-genome hits?
#' @param expand.all.combs logical, should all possible combinations be
#' reported, or just the unique combinations?
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Most useful when building out more complete within-genome
#' hits resulting from whole genome duplication. Typically does not affect closely related regions
#' but can dramatically improve resolution of small blocks, or those that
#' are very diverged, such as ancient WGDs.
#'
#' @return A data.table formatted following typical 'map' specifications.
#' The 'og.id' column will be re-formatted so that each subgraph is now s
#' specified instead of orthofinder orthogroups.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom igraph clusters graph_from_data_frame
#' @export
complete_graph <- function(map,
                           gff,
                           genomeIDs = NULL,
                           ignore.self,
                           expand.all.combs,
                           verbose = T){

  if (is.null(genomeIDs))
    genomeIDs <- unique(gff$genome)

  map <- subset(map,
                genome1 %in% genomeIDs &
                  genome2 %in% genomeIDs)
  gff <- subset(gff, genome %in% genomeIDs)

  gff1 <- data.table(gff)
  setnames(gff1,
           paste0(colnames(gff1), "1"))

  gff2 <- data.table(gff)
  setnames(gff2,
           paste0(colnames(gff2), "2"))

  if (verbose)
    cat("Building graph ... ")

  if (ignore.self) {
    gene2gene <- subset(map,
                        genome1 != genome2)[,c("id1","id2")]
  }else{
    gene2gene <- map[,c("id1","id2")]
  }

  clusters <- clusters(
    graph_from_data_frame(
      gene2gene,
      directed = F)
  )

  df <- with(clusters,
             data.table(
               id = names(membership),
               group = membership))
  gffi <- data.table(gff)
  gffi[,genome := factor(genome,
                         levels = genomeIDs)]
  setkey(gffi, genome, chr, start)
  gffi[,gene.num := 1:nrow(gffi)]

  df <- merge(gffi[,c("id","gene.num")],
              df,
              by = "id")
  df[,id := NULL]
  setnames(df, 1, "id")

  if (verbose)
    cat("Done!\nCompleting subgraphs ... ")

  mdf <- setDT(df)[order(id),
                   unique(id),
                   by = group]
  if (expand.all.combs) {
    out.comb <- mdf[mdf,
                    on = .(group),
                    .(group, x.V1, i.V1),
                    nomatch = 0L,
                    allow.cartesian = TRUE]
  }else{
    out.comb <- mdf[mdf,
                    on = .(group, V1 <= V1),
                    .(group, x.V1, i.V1),
                    nomatch = 0L,
                    allow.cartesian = TRUE]
  }

  if (verbose)
    cat("Done!\nMerging with gff ... ")
  out.comb[,id1 := gffi$id[match(x.V1, gffi$gene.num)]]
  out.comb[,id2 := gffi$id[match(i.V1, gffi$gene.num)]]
  setnames(out.comb, "group", "og.id")
  out.comb[,x.V1 := NULL]
  out.comb[,i.V1 := NULL]
  gffo <- merge(gff1,
                merge(gff2,
                      out.comb,
                      by = "id2"),
                by = "id1")

  if (verbose)
    cat("Done!\n")
  gffo[, block.id := 1]
  gffo[, og.id := paste0("grp_", og.id)]
  return(gffo)
}
