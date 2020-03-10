#' @title finalize_blocks
#'
#' @description
#' \code{finalize_blocks} finalize_blocks
#'
#' @param assigned.map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param MCScanX.s.param numeric of length 1, that specifies the 's'
#' (block size) parameter for MCScanX.
#' @param MCScanX.m.param numeric of length 1, that specifies the 'm'
#' (n. gaps) parameter for MCScanX.
#' @param MCScanX.path file.path, specifying the location of the
#' MCScanX program. This directory must contain the executable
#' '/MCScanX'.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#' @param radii numeric length 3, the radius to search within for syntenic
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
finalize_blocks <- function(assigned.map,
                            radii,
                            gff,
                            genomeIDs,
                            dir.list,
                            MCScanX.s.param = 5,
                            MCScanX.path,
                            MCScanX.m.param = 50,
                            verbose = T){

  split_homologsBySynteny <- function(map){
    zo <- subset(map, hit.type %in% c("paralog", "ortholog"))
    zo[,rank1 := frank(order1, ties.method = "dense"),
       by = list(genome1, genome2, chr1)]
    zo[,rank2 := frank(order2, ties.method = "dense"),
       by = list(genome1, genome2, chr2)]
    spl.map <- split(zo, by = c("genome1","genome2","chr1","chr2"))
    out <- rbindlist(lapply(names(spl.map), function(i){
      z <- spl.map[[i]]
      if(sum(z$hit.type == "ortholog") > 5){
        wh <- find_whichInBuffer(x = z$rank1,
                                 y = z$rank2,
                                 which.in.blk = which(z$hit.type == "ortholog"),
                                 rank.buffer = 100)
        yo <- rep(TRUE, nrow(z))
        yo[-wh] <- FALSE
        z[,syn.with.ortho := yo]
      }else{
        z[,syn.with.ortho := FALSE]
      }
      return(z)
    }))
    intra.o <- subset(out, syn.with.ortho & genome1 == genome2)
    intra.p <- subset(out, !syn.with.ortho & genome1 == genome2)
    inter.b <- subset(map, genome1 != genome2 &
                        hit.type %in% c("paralog", "ortholog"))
    return(list(intra.ortholog = intra.o,
                intra.paralog = intra.p,
                inter.either = inter.b))
  }

  clean_3x <- function(map, genomeIDs, radii = c(100,10,50), verbose = verbose){
    if(verbose)
      cat("\tInitial cleaning - radius =", radii[1],"...",
          length(unique(map$block.id)),"blks")
    ib <- clean_it(map = map,
                   genomeIDs = genomeIDs,
                   rerank = T,
                   n.mappings = 1,
                   radius = radii[1],
                   verbose = F)$map
    if(verbose)
      cat(" -->",length(unique(ib$block.id)),"blks\n")
    ib[,rank1 := frank(start1, ties.method = "dense"),
       by = "block.id"]
    ib[,rank2 := frank(start2, ties.method = "dense"),
       by = "block.id"]
    ib[,genome1 := block.id]
    ib[,genome2 := block.id]

    if(verbose)
      cat("\tCleaning within blocks - radius =", radii[2],"...",
          length(unique(ib$block.id)),"blks")
    ib <- clean_it(map = ib,
                   genomeIDs = unique(ib$block.id),
                   rerank = F,
                   n.mappings = 1,
                   radius = radii[2],
                   verbose = F)$map
    cn <- colnames(map)
    ib <- merge(ib[,c("block.id","id1","id2")],
                map[,cn[cn!="block.id"], with = F],
                by = c("id1","id2"))
    if(verbose)
      cat(" -->",length(unique(ib$block.id)),"blks\n")
    ib[,rank1 := frank(start1, ties.method = "dense"),
       by = list(genome1, genome2, chr1)]
    ib[,rank2 := frank(start2, ties.method = "dense"),
       by = list(genome1, genome2, chr2)]
    ib[,genome1 := paste(block.id)]
    ib[,genome2 := paste(block.id)]
    if(verbose)
      cat("\tFinal splitting within blocks - radius =", radii[3],"...",
          length(unique(ib$block.id)),"blks")
    ib <- clean_it(map = ib,
                   genomeIDs = unique(ib$block.id),
                   rerank = F,
                   n.mappings = 1,
                   radius = radii[3],
                   verbose = F)$map
    ib <- merge(ib[,c("block.id","id1","id2")],
                map[,cn[cn!="block.id"], with = F],
                by = c("id1","id2"))
    if(verbose)
      cat(" -->",length(unique(ib$block.id)),"blks\n")
    return(ib)
  }

  if(verbose)
    cat("Pruning hits to collinear blocks with MCScanX ... ")
  sblast <- subset(assigned.map, genome1 == genome2 & id1 == id2)
  collin.blast <- run_MCScanX(
    blast = assigned.map,
    gff = gff,
    mcscan.dir = dir.list$mcscanx,
    overwrite.output.dir = T,
    genomeIDs = genomeIDs,
    MCScanX.path = MCScanX.path,
    MCScanX.s.param = MCScanX.s.param,
    MCScanX.m.param = MCScanX.m.param,
    verbose = F)
  collin.blast <- rbind(collin.blast,sblast, fill = T)
  collin.blast <- collin.blast[!duplicated(collin.blast[,c("id1","id2","genome1","genome2")]),]

  if(verbose)
    cat("Done!\nCompleting subgraphs ... ")
  compl.blast <- complete_graph(
    map = collin.blast,
    gff = gff,
    ignore.self = T,
    verbose = F)
  compl.blast <- rbind(compl.blast,sblast, fill = T)
  compl.blast <- compl.blast[!duplicated(compl.blast[,c("id1","id2","genome1","genome2")]),]

  compl.blast <- merge(
    compl.blast[,c("genome1","genome2","id1","id2")],
    assigned.map,
    by = c("genome1","genome2","id1","id2"))

  if(verbose)
    cat("Done!\nSplitting hits into three categories ... \n")
  homo.list <- split_homologsBySynteny(compl.blast)

  if(verbose)
    cat("\tn. intra-specific orthologs (or syntenic paralogs):",nrow(homo.list[[1]]),
        "\n\tn. intra-specific paralogs (distant from ortholog blocks):",nrow(homo.list[[2]]),
        "\n\tn. inter-specific orthologs or paralogs:",nrow(homo.list[[3]]),"\n")

  if(verbose)
    cat("Re-building blocks for each group: \n")
  each.blks <- rbindlist(lapply(names(homo.list), function(i){
    if(verbose)
      cat(paste0("\t",i),"... ")
    cl3 <- clean_3x(homo.list[[i]],
                    genomeIDs = genomeIDs,
                    radii = radii,
                    verbose = F)
    if(verbose)
      cat("Done!\n")
    cl3[,blk.type := i]
    return(cl3)
  }), fill = T)

  return(each.blks)
}

