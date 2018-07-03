#' @title Make input metadata for pipe_Diamond2MCScanX
#'
#' @description
#' \code{make_inputFileMatrix} Utility function to build metadata
#'
#' @param map The map object (data.frame or data.table)
#' @param blk The block object (data.frame)
#' @param buffer Numeric, the overlapping distance between two blocks.
#' 0 indicates that blocks that overlap by >=0 should be merged.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details Primarily used in the run_MCScanX pipeline.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export

find_genesWithoutOrthos = function(blast.dir, gff.dir, only.interspecies = T){
  og = readLines(file.path(blast.dir,"Orthogroups.txt"))
  og = lapply(og, function(x) strsplit(x," ")[[1]])
  ons = sapply(og, function(x) x[1])
  names(og)<-ons
  og = lapply(og, function(x) x[-1])
  ng = sapply(og, length)

  if(only.interspecies){
    gff.files = list.files(gff.dir, full.names = T)
    names(gff.files) = gsub(".gff3$","",basename(gff.files))

    parse_gff = function(gff){
      g = suppressWarnings(
        data.table::fread(gff,showProgress = F, verbose = F))
      g = g[g$V3 == "gene",c(9,1,4,5,7)]
      g$V9 = sapply(g$V9, function(x) gsub("Name=","",strsplit(x,";")[[1]][2]))
      data.table::setnames(g, c("id","chr","start","end","strand"))
      return(g)
    }

    gff = rbindlist(lapply(names(gff.files), function(i){
      tmp = parse_gff(gff.files[[i]])
      tmp$genome = i
      tmp$order = frank(tmp[,c("chr","start")], ties.method = "random")
      return(tmp)
    }))

    gs = og[ng>1]
    gs = rbindlist(lapply(names(gs), function(x) data.frame(id = gs[[x]], blk = x,
                                                            stringsAsFactors = F)))
    gs = merge(gs, gff, by = "id")
    spl = split(gs, gs$blk)

    ngenomes.perblock = sapply(spl, simplify = T, USE.NAMES = T, function(x) length(unique(x$genome)))
    single.genome.blocks = names(ngenomes.perblock)[ngenomes.perblock == 1]
    out = unlist(og[ng == 1 | names(og) %in% single.genome.blocks])
  }else{
    out = unlist(og[ng==1])
  }

  return(out)

}
