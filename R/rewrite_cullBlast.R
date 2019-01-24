#' @title Cull blast hits by a SpatialPolygons object
#'
#' @description
#' \code{cull_blastInPolygon} Cull blast hits by a SpatialPolygons object
#'
#' @param map The map data.frame or data.table
#' @param blast Some blast data that needs to be culled.
#' Has the same format as the map data, but does not need to have
#' a 'block.id' column.
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param gene.index data.table with a dictionary conneccting
#' the gff annotation and orthofinder gene IDs
#' @param species.mappings data.table with metadata for the
#' orthofinder run. Should come directly from build_syntenicBlocks.
#' @param rerank.gff Logical, should the rank order of genes in the
#' gff be re-calculated?
#' @param rank.buffer The size of the buffer (in gene order ranks)
#' to extend the chull around each block.
#' @param plotit Should buffer_blkChull plots be made?
#' @param tmp.dir file path, pointing to directory that contains
#'  blast files that should be overwritten by culled hits
#' @param n.cores numeric, the number or parallel processes to run.
#' @param verbose logical, should
#' @param ... Plotting arguments passed on to buffer_blkChull and
#' then to plot.
#' @details None yet

#' @return Nothing to R, but re-rewrites the blast files.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
rewrite_cullBlast <- function(gff,
                              map,
                              species.mappings,
                              gene.index,
                              rank.buffer,
                              tmp.dir,
                              rerank.gff = T,
                              verbose = T,
                              n.cores = 1,
                              plotit = F){

  # -- Invert map ids so that both blast directions can be merged
  m.tmp1 <- data.table(map)
  m.tmp2 <- data.table(map)
  setnames(m.tmp2, gsub("1", "xxx", colnames(m.tmp2)))
  setnames(m.tmp2, gsub("2", "1", colnames(m.tmp2)))
  setnames(m.tmp2, gsub("xxx", "2", colnames(m.tmp2)))
  m.tmp <- rbind(m.tmp1, m.tmp2)
  map <- m.tmp[!duplicated(m.tmp),]

  # -- load the map data, split by genome
  map$unique.genome <- with(map,
                            paste(genome1, genome2))
  maps <- split(map, "unique.genome")

  # -- load and re-rank the gff file
  if(rerank.gff){
    gff[,rank := frank(start, ties.method = "dense"),
        by = list(genome, chr)]
  }

  # -- add gene index data into the gff
  setkey(gff, id)
  setkey(gene.index, id)
  gff <- merge(gff, gene.index)

  # -- make two gff objects to merge
  gff1 <- data.table(gff)
  gff2 <- data.table(gff)
  setnames(gff1, paste0(colnames(gff1), "1"))
  setnames(gff2, paste0(colnames(gff2), "2"))
  setkey(gff1, gene.num1)
  setkey(gff2, gene.num2)

  for(x in species.mappings$filename){
    cat("running", x, "...\n\t")

    # -- read in the geneIDs of the blast file
    tmp <- fread(x,
                 select = 1:2,
                 col.names = paste0("gene.num",1:2),
                 key = "gene.num2")

    # -- merge with the gff for both ids
    t1 <- merge(tmp, gff2)
    setkey(t1, gene.num1)
    blast.x <- merge(t1, gff1)
    j <- paste(blast.x$genome1[1],
               blast.x$genome2[1])
    map.x <- maps[[j]]

    # -- split the map and blast datasets
    map.x$unique <- with(map.x,
                        paste(genome1, genome2,
                              chr1, chr2))
    blast.x$unique <- with(blast.x,
                          paste(genome1, genome2,
                                chr1, chr2))

    if(verbose)
      cat(map.x$genome1[1], "-->", map.x$genome2[1],
          paste0("(", nrow(map.x), "/", nrow(blast.x), ") ... "))

    spl.map <- split(map.x, "unique")
    spl.blast <- split(blast.x, "unique")

    # -- For each unique genome / chromosome, cull blasts to chull+buffer regions
    out <- rbindlist(mclapply(names(spl.map), mc.cores = n.cores, function(i){
      blast.x <- spl.blast[[i]]
      map.x <- spl.map[[i]]
      cull.blast <- cull_blastInPolygon(map = map.x,
                                       blast = blast.x,
                                       gff = gff,
                                       rank.buffer = rank.buffer,
                                       plotit = plotit)
      return(cull.blast)
    }))

    # -- Reformat blast ids
    bl.out <- data.table(blast.x)
    setkey(bl.out, id1, id2)
    genes.out <- data.table(out)
    setkey(genes.out, id1, id2)
    out.bl <- merge(bl.out[,c("id1","id2","gene.num1","gene.num2")],
                    genes.out)[,3:4]
    setnames(out.bl, c("V1","V2"))

    # -- Read in full blast data and write culled blast
    dat.in <- merge(fread(x), out.bl)
    if(verbose)
      cat("culled to", nrow(dat.in), "\n")

    out.path <- file.path(tmp.dir,
                          basename(x))
    write.table(dat.in,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = out.path)
  }
}
