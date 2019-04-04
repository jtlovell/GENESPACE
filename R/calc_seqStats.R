#' @title build pwBlocks
#'
#' @description
#' \code{calc_seqStats} build_pwBlocks
#'
#' @param geneIDs dir.list
#' @param orthonet genomeIDs
#' @param cds.fastas path to the directory where MCScanX will be run
#' @param pep.fastas parameters to supply MCScanX
#' @param tmp.dir temporary directory
#' @param make.tree logical, should trees be generated?
#' @param pal2nal.tool path to pal2nal
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
#' @importFrom Biostrings writeXStringSet
#' @export
calc_seqStats <- function(geneIDs = NULL,
                          orthonet = NULL,
                          cds.fastas,
                          pep.fastas,
                          tmp.dir,
                          make.tree,
                          pal2nal.tool = "/Users/jlovell/Documents/comparative_genomics/programs/pal2nal.v14/pal2nal.pl",
                          verbose = T){

  if (is.null(geneIDs) & is.null(orthonet))
    stop("either geneIDs or a orthonet object must be specified")

  if (is.null(geneIDs) & is.data.table(orthonet)) {
    geneIDs <- lapply(1:nrow(orthonet), function(i) as.character(unlist(orthonet[i,-1, with = F])))
    names(geneIDs) <- orthonet$og
  }else{
    if(is.list(orthonet) & is.data.table(orthonet[[1]])){
      geneIDs <- unlist(lapply(orthonet, function(x)
        lapply(1:nrow(x), function(i) as.character(unlist(x[i,-1, with = F])))),
        recursive = F)
      names(geneIDs) <- unlist(lapply(orthonet, function(x) x$og))
    }else{
      if(is.character(geneIDs)){
        geneIDs <- list(geneIDs)
        names(geneIDs) <- "NA"
      }else{
        stop("geneIDs must be a character vector.\n")
      }
    }
  }

  if(!all(unlist(geneIDs) %in% names(cds.fastas) &
          unlist(geneIDs) %in% names(pep.fastas)))
    stop("all geneIDs must be present in the names of both cds and peptide fastas.\n")

  if(verbose & is.null(orthonet))
    cat("Writing fasta files ... ")
  if(verbose)
    cat("Writing fasta files for",length(geneIDs),"orthogroups ... ")

  names(geneIDs) <- gsub("[^[:alnum:]]","",names(geneIDs))
  cds.tmp.files <- file.path(tmp.dir, paste0(names(geneIDs),".cds.tmp.fa"))
  names(cds.tmp.files) <- names(geneIDs)
  pep.tmp.files <- file.path(tmp.dir, paste0(names(geneIDs),".pep.tmp.fa"))
  names(pep.tmp.files) <- names(geneIDs)

  cds.list <- sapply(geneIDs, function(x)  cds.fastas[x])
  pep.list <- sapply(geneIDs, function(x)  pep.fastas[x])

  for(i in names(geneIDs)){
    writeXStringSet(cds.list[[i]], filepath = cds.tmp.files[i])
    writeXStringSet(pep.list[[i]], filepath = pep.tmp.files[i])
  }

  if(verbose)
    cat("Done\n")

  owd <- getwd()
  setwd(tmp.dir)

  out <- mclapply(names(geneIDs), mc.cores = 1, function(i){
    out <- calc_selectionStats(pep.file = pep.tmp.files[i],
                               cds.file = cds.tmp.files[i],
                               tmp.dir = tmp.dir,
                               codeml.msa.file = file.path(tmp.dir,paste0(i,".codon.aln")),
                               msa.clu = file.path(tmp.dir,paste0(i,".msa.clu")),
                               cds.msa.fa = file.path(tmp.dir,paste0(i,".msa.fa")),
                               pal2nal.tool = pal2nal.tool)

    if(make.tree){
      algn.fa <- out$files$cds.msa.fa
      tre <- system(paste("fasttree -nt -quiet -nopr", algn.fa), intern = T)
    }else{
      tre <- NULL
    }
    return(list(stats = out$stats,
                tree = tre))
  })

  stats <- rbindlist(lapply(out, function(x) x$stats))
  stats$og <- rep(names(geneIDs), sapply(lapply(out, function(x) x$stats), nrow))
  trees <- do.call(c, lapply(out, function(x) x$tree))
  names(trees) <- names(geneIDs)
  setwd(owd)
  return(list(stats = stats, trees = trees))
}
