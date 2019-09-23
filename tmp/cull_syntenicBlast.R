#' @title Cull blast to syntenic hist
#'
#' @description
#' \code{cull_syntenicBlast} Subset blast hits to syntenic regions.
#'
#' @param map The map data.frame or data.table
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param blast the blast dataset to screen for syntenic hits
#' @param rank.buffer The buffer, in gene rank order.
#' @param verbose Logical, should updates be printed
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
cull_syntenicBlast <- function(map,
                               blast,
                               gff,
                               rank.buffer,
                               verbose = T){
  if (verbose)
    cat("Dropping chromosome combinations in blast not found in map... ")
  blast.in <- cull_blast2MapChr(blast = blast,
                                map = map)
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Dropping blast hits in the map already ... ")
  t.blast <- cull_blast2NewIds(blast = blast.in,
                               map = map)
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Making new blast and map with gff-based ranks ... ")
  id.names <- c("genome1", "genome2",
                "id1", "id2",
                "chr1", "chr2",
                "start1", "start2",
                "end1", "end2")
  all.ids <- rbind(map[, id.names, with = F],
                   t.blast[, id.names, with = F])
  all.ids <- all.ids[!duplicated(all.ids), ]
  all.ids[ ,rank1 := frank(start1, ties.method = "dense"),
           by = list(genome1, chr1)]
  all.ids[ ,rank2 := frank(start2, ties.method = "dense"),
           by = list(genome2, chr2)]
  setkey(all.ids, id1, id2)
  setkey(map, id1, id2)
  setkey(blast, id1, id2)
  r.map <- merge(all.ids, map[ ,c("id1", "id2")])
  r.blast <- merge(all.ids, t.blast[ ,c("id1", "id2")])

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Finding neighbors to mappings in blast hits ... completed:\n")
  r.map[, unique.genome := paste(genome1, genome2)]
  r.blast[, unique.genome := paste(genome1, genome2)]
  r.map[, unique.chr := paste(unique.genome, chr1, chr2)]
  r.blast[, unique.chr := paste(unique.genome, chr1, chr2)]

  r.map <- merge(r.map, map[,c("id1","id2","block.id")], by = c("id1","id2"))

  ids2keep <- find_hitsInBuffer(map = r.map,
                                blast = r.blast,
                                rank.buffer = rank.buffer,
                                verbose = verbose)

  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Reformatting blast output ... ")

  setkey(blast, id1, id2)

  out.blast <- merge(ids2keep, blast)
  re.out <- rerank_fromIDs(map = out.blast,
                           gff = gff)

  if (verbose)
    #######################################################
  cat("Done!\n")
  return(re.out)
}
