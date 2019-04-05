#' @title calc_seqStats
#'
#' @description
#' \code{calc_seqStats} calc_seqStats
#'
#' @param geneIDs dir.list
#' @param orthonet genomeIDs
#' @param dir.list directory list
#' @param make.tree logical, should trees be generated?
#' @param pal2nal.tool path to pal2nal
#' @param n.cores number of parallel processes.
#' @param map optional map object
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
                          map,
                          dir.list,
                          make.tree = TRUE,
                          n.cores = 1,
                          clean = T,
                          pal2nal.tool = "/Users/jlovell/Documents/comparative_genomics/programs/pal2nal.v14/pal2nal.pl",
                          verbose = T){



  tmp.dir <- dir.list$tmp

  if(clean){
    unlink(dirs$tmp, recursive = T)
    dir.create(tmp.dir)

  }

  # -- Check an conver the gene IDs, either from a simple vector or an orthonet obj.
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


  names(geneIDs) <- gsub("[^[:alnum:]]","",names(geneIDs))

  cds.tmp.files <- file.path(tmp.dir, paste0(names(geneIDs),".cds.tmp.fa"))
  names(cds.tmp.files) <- names(geneIDs)
  pep.tmp.files <- file.path(tmp.dir, paste0(names(geneIDs),".pep.tmp.fa"))
  names(pep.tmp.files) <- names(geneIDs)

  if(verbose & is.null(orthonet))
    cat("Writing fasta files ... ")
  if(verbose)
    cat("Writing fasta files for",length(geneIDs),"orthogroups ... ")

  cds.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readDNAStringSet(file.path(dir.list$cds,
                               paste0(i, ".fa")))))

  pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readAAStringSet(file.path(dir.list$peptide,
                              paste0(i, ".fa")))))

  m <- map[map$id1 %in% unique(as.character(unlist(geneIDs))) &
             map$id1 %in% unique(as.character(unlist(geneIDs))),]
  mi <- data.table(id = c(m$id1, m$id2),
                   gn = c(m$gn1, m$gn2))
  mi <- mi[!duplicated(mi),]
  cds.fastas <- cds.fastas[mi$id]
  pep.fastas <- pep.fastas[mi$id]
  names(cds.fastas) <- mi$gn
  names(pep.fastas) <- mi$gn

  geneIDs <- lapply(geneIDs, function(x) mi$gn[mi$id %in% x])

  if(!all(unlist(geneIDs) %in% names(cds.fastas) &
          unlist(geneIDs) %in% names(pep.fastas)))
    stop("all geneIDs must be present in the names of both cds and peptide fastas.\n")

  cds.list <- sapply(geneIDs, function(x)  cds.fastas[x])
  pep.list <- sapply(geneIDs, function(x)  pep.fastas[x])

  for(i in names(geneIDs)){
    writeXStringSet(cds.list[[i]], filepath = cds.tmp.files[i])
    writeXStringSet(pep.list[[i]], filepath = pep.tmp.files[i])
  }

  if(verbose)
    cat("Done\nCalculating selection statistics using",n.cores,"parallel threads ... ")

  owd <- getwd()

  out <- mclapply(names(geneIDs), mc.cores = n.cores, function(i){

    tmp.dir2 <- file.path(tmp.dir, i)
    dir.create(tmp.dir2)
    setwd(tmp.dir2)
    out <- calc_selectionStats(pep.file = pep.tmp.files[i],
                               cds.file = cds.tmp.files[i],
                               tmp.dir = tmp.dir2,
                               pal2nal.tool = pal2nal.tool)
    if(clean){
      unlink(tmp.dir2, recursive = T)
    }
    return(out)
  })
  stats <- rbindlist(lapply(out, function(x) x$stats))
  setnames(stats, c("id1","id2"),c("gn1","gn2"))
  trees <- do.call(c, lapply(out, function(x) x$tree))
  names(trees) <- names(geneIDs)
  setwd(owd)


  mapo <- map[,c("genome1","genome2","id1","id2","gn1","gn2",
                 "chr1","chr2","start1","start2","end1","end2",
                 "strand1","strand2","score","og1","block.id")]
  map1 <- mapo[,c(2,1,4,3,6,5,8,7,10,9,12,11,14,13,15:17)]
  setnames(map1,colnames(mapo))
  mapi <- rbind(mapo,map1)
  mapi <- mapi[!duplicated(mapi),]
  stats <- merge(mapi, stats, by = c("gn1","gn2"))

  for(i in 1:length(trees)){
    for(j in 1:nrow(mi)){
      trees[[i]] <- gsub(mi[j,2],mi[j,1],trees[[i]] )
    }
  }


  if(verbose)
    cat("Done\n")
  return(list(stats = stats, trees = trees))
}
