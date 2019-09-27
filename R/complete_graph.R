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
                           keep.all.self.hits = T,
                           verbose = T){

  if (is.null(genomeIDs))
    genomeIDs <- unique(gff$genome)
  map.in <- data.table(map)
  map <- subset(map,
                genome1 %in% genomeIDs &
                  genome2 %in% genomeIDs)
  gff[,gn := paste0("gene",1:nrow(gff))]
  map[,gn1 := gff$gn[match(paste(genome1, id1), paste(gff$genome, gff$id))]]
  map[,gn2 := gff$gn[match(paste(genome2, id2), paste(gff$genome, gff$id))]]

  if (verbose)
    cat("Building graph ... ")

  if (ignore.self) {
    gene2gene <- subset(map,
                        genome1 != genome2)
  }else{
    gene2gene <- data.table(map)
  }
  gene2gene <- gene2gene[,c("gn1","gn2")]

  clusters <- clusters(
    graph_from_data_frame(
      gene2gene,
      directed = F)
  )

  df <- with(clusters,
             data.table(
               id = names(membership),
               group = membership))
  if (verbose)
    cat("Done!\nCompleting subgraphs ... ")

  out.comb <- merge(df,
                    df,
                    by = "group",
                    allow.cartesian = T)
  setnames(out.comb, c("og.id", "gn1","gn2"))
  out.comb[,id1 := gff$id[match(gn1, gff$gn)]]
  out.comb[,id2 := gff$id[match(gn2, gff$gn)]]

  if (verbose)
    cat("Done!\nMerging with gff ... ")

  gff1 <- data.table(gff)
  setnames(gff1,
           paste0(colnames(gff1), "1"))

  gff2 <- data.table(gff)
  setnames(gff2,
           paste0(colnames(gff2), "2"))

  oc <- out.comb[,c("gn1","gn2")]
  if(keep.all.self.hits){
    ugn <- unique(unlist(oc))
    ugi <- unique(unlist(subset(oc, gn1 == gn2)))
    ugn <- ugn[!ugn %in% ugn]
    oc <- rbind(oc,
                data.table(gn1 = ugn,
                           gn2 = ugn))
  }

  gffo <- merge(gff1,
                merge(gff2,
                      oc,
                      by = "gn2"),
                by = "gn1")
  if (verbose)
    cat("Done!\n")
  gffo[, block.id := 1]

  return(gffo)
}
