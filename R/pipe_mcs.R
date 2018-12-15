#' @title Run MCScanX program
#'
#' @description
#' \code{run_MCScanX} Wrapper to run MCScanX program
#'
#' @param blast blast results data.table
#' @param gff concatenated gff data.table
#' @param mcscan.dir path to the directory where MCScanX will be run
#' @param mcscan.param parameters to supply MCScanX
#' @param silent.mcs logical, should MCScanX be run silently?
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
pipe_mcs <- function(blast,
                     gff,
                     mcscan.dir,
                     mcscan.param,
                     silent.mcs = T,
                     verbose = T){
  if(verbose)
    cat("Culling",nrow(blast),"BLAST hits by MCScanX\n\t")

  blast$unique <- with(blast, paste(genome1, genome2))
  if("block.id" %in% colnames(blast))
    blast$block.id<-NULL

  spl = split.data.table(blast, "unique")

  out <- rbindlist(lapply(spl, function(x){
    genomes = c(x$genome1[1],x$genome2[2])

    if(verbose)
      cat(genomes[1],"-->", genomes[2],
          paste0("(initial hits = ",nrow(x),")"))

    gff.x <- gff[gff$genome %in% genomes,]
    tmp = run_mcs(blast = x,
                  gff = gff.x,
                  mcscan.dir = mcscan.dir,
                  mcscan.param = mcscan.param,
                  silent.mcs = silent.mcs)

    tmp$block.id<-with(tmp, paste0(unique, block.id))

    if(verbose)
      cat(" culled hits = ",nrow(tmp),"\n\t")

    return(tmp)
  }))

  out$block.id<- as.numeric(as.factor(out$block.id))
  if(verbose)
    cat("Done!\n")

  return(make_blocks(out))
}


