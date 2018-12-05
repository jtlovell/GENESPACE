#' @title Find the coordinates of blocks, based on the original annotations
#'
#' @description
#' \code{make_genomeDB} Reduces the total number of block breakpoints, so that, when
#' blocks are concatenated, there are fewer small orphan blocks.
#'
#' @param blk The block object (data.frame/data.table)
#' @param map The map object (data.frame/data.table)
#' @param genomeIDs Character vector indicating the genome IDs to consider
#' @param assembly.dir The directory containing the assembly fasta sequences
#' @param gff.dir The directory containing gff3-formatted annotations
#' @param peptide.dir The directory containing peptide fasta sequences
#' @param cds.dir The directory containing cds fasta sequences
#' @param cds.dir The directory containing full-length fasta sequences
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details ...
#' @return A list of block datasets
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
make_genomeDB = function(blk.coords,
                         gff.dir,
                         assembly.dir,
                         peptide.dir,
                         cds.dir,
                         transcript.dir,
                         block.dir,
                         verbose = T){

  if(verbose)
    cat("Parsing gff files\n")
  gff.files <- list.files(gff.dir,
                          full.names = T)
  names(gff.files) <- gsub(".gff3$", "",
                           basename(gff.files))

  parse_gff <- function(gff){
    g <- suppressWarnings(
      data.table::fread(gff,
                        showProgress = F,
                        verbose = F))
    g <- g[g$V3 == "gene", c(9, 1, 4, 5, 7)]
    g$V9 <- sapply(g$V9, function(x) gsub("Name=", "",
                                          strsplit(x, ";")[[1]][2]))
    data.table::setnames(g, c("id", "chr", "start", "end", "strand"))
    return(g)
  }

  gff.list <- sapply(names(gff.files), simplify = F, USE.NAMES = T, function(i){
    tmp <- parse_gff(gff.files[[i]])
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "dense")
    return(tmp)
  })

  if(verbose)
    cat("Loading peptide fasta files\n")
  peptide.fastas = sapply(genomeIDs, simplify = F, USE.NAMES = T,function(i){
    Biostrings::readAAStringSet(file.path(peptide.dir,paste0(i,".fa")))
  })

  if(verbose)
    cat("Loading transcript fasta files\n")
  trs.fastas = sapply(genomeIDs, simplify = F, USE.NAMES = T,function(i){
    Biostrings::readDNAStringSet(file.path(transcript.dir,paste0(i,".fa")))
  })

  if(verbose)
    cat("Loading cds fasta files\n")
  cds.fastas = sapply(genomeIDs, simplify = F, USE.NAMES = T,function(i){
    Biostrings::readDNAStringSet(file.path(cds.dir,paste0(i,".fa")))
  })


  if(verbose)
    cat("Building block databases and storing them in", block.dir,"\n\tCompleted:")
  if(!file.exists(block.dir)){
    system(paste("mkdir", block.dir))
  }
  spl = split(blk.coords, blk.coords$block.id)
  out.data = rbindlist(lapply(1:length(spl), function(i){
    cat(paste0(i,", "))
    blk.id = paste0("block_",i)
    blk.path = file.path(block.dir,blk.id)
    blk.assem = file.path(blk.path,"assemblies")
    blk.pep = file.path(blk.path,"peptide")
    blk.cds = file.path(blk.path,"cds")
    blk.trs = file.path(blk.path,"transcript")
    blk.tmp = file.path(blk.path,"tmp")

    if(file.exists(blk.path)){
      system(paste("rm -rf", blk.path))
    }
    system(paste("mkdir", blk.path))
    system(paste("mkdir", blk.assem))
    system(paste("mkdir", blk.pep))
    system(paste("mkdir", blk.cds))
    system(paste("mkdir", blk.trs))
    system(paste("mkdir", blk.tmp))

    x = spl[[i]]
    x$unique = paste("block",x$block.id,x$genome, x$block.id,sep = "_")
    xo = x[!x$genome %in% genomeIDs[1:2],]
    x = x[x$genome %in% genomeIDs[1:2],]
    gff.l = rbindlist(lapply(1:nrow(x), function(j){
      x0 = x[j,]
      y = gff.list[[x0$genome]]
      y = y[with(y, chr == x0$chr & end >= x0$start & start <= x0$end),]
      ygenes = y$id
      ycds = cds.fastas[[x0$genome]][ygenes]
      ypep = peptide.fastas[[x0$genome]][ygenes]
      ytrs = trs.fastas[[x0$genome]][ygenes]
      x0$pep.path = file.path(blk.pep, paste0(x0$unique,".fa"))
      x0$trs.path = file.path(blk.cds, paste0(x0$unique,".fa"))
      x0$cds.path = file.path(blk.trs, paste0(x0$unique,".fa"))
      Biostrings::writeXStringSet(ypep, filepath =  x0$pep.path)
      Biostrings::writeXStringSet(ycds, filepath =  x0$cds.path)
      Biostrings::writeXStringSet(ytrs, filepath =  x0$trs.path)

      xbed = x0[,c("chr","start","end","unique")]
      bedf = file.path(blk.tmp,paste0(x0$unique,".bed"))
      faf = file.path(assembly.dir,paste0(x0$genome,".fa"))
      fafo = file.path(blk.assem,paste0(x0$unique,".fa"))
      x0$assem.path = fafo
      write.table(xbed, file=bedf,
                  quote=F, sep="\t", row.names=F, col.names=F)
      system(paste("bedtools getfasta -fi",
                   faf, "-bed", bedf,"-name -fo",fafo))
      return(x0)
    }))

    gff.o = rbindlist(lapply(unique(xo$genome), function(j){
      x0 = xo[xo$genome == j,]
      u = with(x0, paste("block",block.id[1],genome[1], sep = "_"))
      y = gff.list[[j]]
      y = rbindlist(lapply(1:nrow(x0), function(k){
        y[with(y, chr == x0$chr[k] & end >= x0$start[k] & start <= x0$end[k]),]
      }))
      ygenes = y$id
      ycds = cds.fastas[[j]][ygenes]
      ypep = peptide.fastas[[j]][ygenes]
      ytrs = trs.fastas[[j]][ygenes]
      x0$pep.path = file.path(blk.pep, paste0(u,".fa"))
      x0$trs.path = file.path(blk.cds, paste0(u,".fa"))
      x0$cds.path = file.path(blk.trs, paste0(u,".fa"))
      Biostrings::writeXStringSet(ypep, filepath =  x0$pep.path[1])
      Biostrings::writeXStringSet(ycds, filepath =  x0$cds.path[1])
      Biostrings::writeXStringSet(ytrs, filepath =  x0$trs.path[1])

      xbed = x0[,c("chr","start","end","unique")]
      bedf = file.path(blk.tmp,paste0(u,".bed"))
      faf = file.path(assembly.dir,paste0(j,".fa"))
      fafo = file.path(blk.assem,paste0(u,".fa"))
      x0$assem.path = fafo
      write.table(xbed, file=bedf,
                  quote=F, sep="\t", row.names=F, col.names=F)
      system(paste("bedtools getfasta -fi",
                   faf, "-bed", bedf,"-name -fo",fafo))
      return(x0)
    }))


    return(rbind(gff.l,gff.o))
  }))

  if(verbose)
    cat("\nDone!\n")
  return(out.data)
}



