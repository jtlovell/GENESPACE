#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs2} Subset blast hits to syntenic regions and
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
#' @importFrom sp CRS SpatialPointsDataFrame over Polygon Polygons SpatialPolygons
#' @importFrom rgeos gBuffer
#' @importFrom grDevices chull
#' @export
find_syntenicOrthogs2 <- function(map,
                                 dir.list,
                                 gff,
                                 gene.index,
                                 species.mappings,
                                 n.cores = 1,
                                 plotit = F,
                                 rerank.gff = T,
                                 buffer,
                                 verbose = T){

  parse_blastRaw <- function(species.mappings,
                             gff,
                             rerank.gff = T,
                             gene.index,
                             map,
                             buffer,
                             verbose = T,
                             n.cores = 1,
                             tmp.dir,
                             plotit = F){
    m.tmp1 = data.table(map)
    m.tmp2 = data.table(map)
    setnames(m.tmp2, gsub("1","xxx", colnames(m.tmp2)))
    setnames(m.tmp2, gsub("2","1", colnames(m.tmp2)))
    setnames(m.tmp2, gsub("xxx","2", colnames(m.tmp2)))
    m.tmp = rbind(m.tmp1, m.tmp2)
    map = m.tmp[!duplicated(m.tmp),]
    ## load the map data, split by genome
    map$unique.genome = with(map,
                             paste(genome1, genome2))
    maps <- split(map, "unique.genome")

    ## load the species mapping data
    sm = species.mappings

    ## load and re-rank the gff file
    if(rerank.gff){
      gff[,rank := frank(start, ties.method = "dense"),
          by = list(genome, chr)]
    }

    ## add gene index data into the gff
    setkey(gff, id)
    setkey(gene.index, id)
    gff <- merge(gff, gene.index)

    ## make two gff objects to merge
    gff1 = data.table(gff)
    gff2 = data.table(gff)
    setnames(gff1, paste0(colnames(gff1),"1"))
    setnames(gff2, paste0(colnames(gff2),"2"))
    setkey(gff1, gene.num1)
    setkey(gff2, gene.num2)

    for(x in sm$filename){
      cat("running",x,"...\n\t")
      tmp = fread(x, select = 1:2, col.names = paste0("gene.num",1:2), key = "gene.num2")
      t1 = merge(tmp, gff2)
      setkey(t1, gene.num1)
      blast.x = merge(t1, gff1)
      j = paste(blast.x$genome1[1], blast.x$genome2[1])
      map.x = maps[[j]]
      map.x$unique = with(map.x,
                          paste(genome1, genome2,
                                chr1, chr2))
      blast.x$unique = with(blast.x,
                            paste(genome1, genome2,
                                  chr1, chr2))

      if(verbose)
        cat(map.x$genome1[1],"-->", map.x$genome2[1],
            paste0("(", nrow(map.x),"/",nrow(blast.x),") ... "))

      spl.map = split(map.x, "unique")
      spl.blast = split(blast.x, "unique")

      out <- rbindlist(mclapply(names(spl.map), mc.cores = n.cores, function(i){
        blast.x = spl.blast[[i]]
        map.x = spl.map[[i]]
        cull.blast = cull_blastbyMap(map = map.x,
                                     blast = blast.x,
                                     gff = gff,
                                     buffer = buffer,
                                     plotit = plotit)
        return(cull.blast)
      }))
      bl.out = data.table(blast.x)
      setkey(bl.out, id1, id2)
      genes.out = data.table(out)
      setkey(genes.out, id1, id2)
      out.bl = merge(bl.out[,c("id1","id2","gene.num1","gene.num2")], genes.out)[,3:4]
      setnames(out.bl,c("V1","V2"))
      dat.in = merge(fread(x),out.bl)
      if(verbose)
        cat("culled to", nrow(dat.in),"\n")

      out.path = file.path(tmp.dir, basename(x))
      write.table(dat.in,
                  sep = "\t",
                  row.names = F,
                  col.names = F,
                  quote = F,
                  file = out.path)
    }
  }


  if (verbose)
    cat("Cleaning out tmp and culled blast directories ... ")

  unlink(dir.list$tmp, recursive = T)
  dir.create(dir.list$tmp)

  unlink(dir.list$cull.blast, recursive = T)
  dir.create(dir.list$cull.blast)

  if (verbose)
    cat("Done!\n")

  if (verbose)
    cat("Copying orthofiner results to tmp directory ... ")

  files <- list.files(dir.list$blast,
                      full.names = T)
  files = files[-grep("Blast|Orthogroups.txt", files)]
  nu <- file.copy(files,
                  dir.list$tmp)
  if (verbose)
    cat("Done!\nParsing blast results and moving to tmp directory\n")

  nothing = parse_blastRaw(species.mappings = species.mappings,
                           gff = gff,
                           rerank.gff = rerank.gff,
                           gene.index = gene.index,
                           map = map,
                           buffer = buffer,
                           n.cores = n.cores,
                           tmp.dir = dir.list$tmp)
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
  gff$rank <- NULL
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

  return(list(blast = all.blast, of.info = of.blast))
}
