#' @title Clean form_syntenicBlocks output
#'
#' @description
#' \code{clean_blocks} Clusters hits and drops low-confidence blocks.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param rerank logical, should the ranks be re-calculated prior to cleaning?
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#' @param complete.graphs logical, should complete_graph be run?
#' @param ignore.self logical, should the subgraphs be built agnostic to
#' the existence of within-genome hits? Passed on to complete_graph.
#' @param expand.all.combs logical, should all possible combinations be
#' reported, or just the unique combinations? Passed on to complete_graph.
#' @param n.mappings numeric, the number of mappings required within a
#' given radius. Length must match that of radius
#' @param radius numeric, the radius to search within for syntenic
#' mappings. Length must match that of n.mappings.
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A list of length 2, block and map, as output by make_blocks.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom igraph clusters graph_from_data_frame
#' @export
clean_blocks <- function(map,
                         gff,
                         genomeIDs,
                         rerank = TRUE,
                         radius = 100,
                         n.mappings = 10,
                         verbose = TRUE,
                         expand.all.combs = TRUE,
                         complete.graphs = T,
                         ignore.self = T){

  #######################################################
  if ((length(radius) != length(n.mappings)))
    stop("radius and n.mappings must be of same length\n")
  #######################################################

  #######################################################
  # -- Prep the map object
  map <- subset(map, genome1 %in% genomeIDs &
                  genome2 %in% genomeIDs)
  map[, unique := paste(genome1, genome2, chr1, chr2)]
  map[, unique.genome := paste(genome1, genome2)]

  setkey(map, chr1, chr2, start1, start2)
  #######################################################

  #######################################################
  # -- Iteratively (or not) run the cleaning
  for (i in 1:length(n.mappings)) {
    n.map <- n.mappings[i]
    rad <- radius[i]
    if (verbose) {
      if (length(n.mappings) == 1) {
        cat("Cleaning mappings to", n.map, "hits within",
            rad, "gene-rank radius\n")
      }else{
        cat(paste0("Step", i, ":"),
            "blocks must have", n.map, "hits within",
            rad, "gene-rank radius\n")
      }
    }

    cleaned <- clean_it(
      map = map,
      genomeIDs = genomeIDs,
      rerank = rerank,
      radius = rad,
      n.mappings = n.map,
      verbose = verbose)
    map <- data.table(cleaned$map)

    if (verbose)
      cat("Cleaned n blocks / mappings =",
          nrow(cleaned$block),
          "/",
          nrow(cleaned$map), "\n")
  }
  #######################################################

  if (complete.graphs) {
    if (is.null(gff))
      stop("Must specify gff in order to complete.graphs\n")
    comp <- complete_graph(map = map,
                          gff = gff,
                          genomeIDs = genomeIDs,
                          verbose = verbose,
                          expand.all.combs = expand.all.combs,
                          ignore.self = ignore.self)

    comb <- data.table(
      rbind(
        t(combn(genomeIDs, 2, simplify = T)),
        t(sapply(genomeIDs, rep, 2))))
    setnames(comb, c("genome1", "genome2"))

    comp <- merge(
      comp,
      comb,
      by = c("genome1", "genome2"))

    if (verbose)
      cat("Re-building blocks with radius =",
          radius[length(radius)],
          "and n.mappings =",
          n.mappings[length(n.mappings)],
          "...\n")

    cleaned <- clean_it(
      map = comp,
      genomeIDs = genomeIDs,
      rerank = TRUE,
      radius = radius[length(radius)],
      n.mappings = n.mappings[length(n.mappings)],
      verbose = verbose)
    map <- data.table(cleaned$map)
    if (verbose)
      cat("\tDone!\n")
  }

  map[,block.id := paste0("blk_",
                          as.numeric(
                            as.factor(
                              paste(genome1, genome2,
                                    chr1, chr2,
                                    block.id))))]
  #######################################################

  if (verbose)
    cat("Compiling block coordinates ... ")
  alo <- make_blocks(map,
                     clean.columns = F,
                     rerank = T,
                     rename.blocks = F)
  if (verbose)
    cat("Done!\n")
  #######################################################
  return(alo)
}
