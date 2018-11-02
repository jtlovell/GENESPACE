#' @title Run MCScanX program
#'
#' @description
#' \code{run_MCScanX} Wrapper to run MCScanX program
#'
#' @param blast.results R object contain the blast results
#' @param MCScanX.params Parameters to pass to MCScanX
#' @param mcscanx.input.dir Directory containing the MCScanX-formatted mapping files
#' @param abbrevs Genome abbreviations
#' @param verbose Logical, should updates be printed?
#' @param ... Not currently in use
#' @details See pipe_Diamond2MCScanX for more information.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
run_MCScanX = function(blast.results,
                       abbrevs,
                       mcscanx.input.dir,
                       verbose = T,
                       MCScanX.params = "-a -s 5 -m 25"){

  if(file.exists(mcscanx.input.dir)){
    system(paste("rm -r",mcscanx.input.dir))
  }
  system(paste("mkdir",mcscanx.input.dir))

  br = data.frame(blast.results, stringsAsFactors = F)
  bs = br[!duplicated(br[,c("genome1","genome2")]),c("genome1","genome2")]
  g1 = names(table(bs$genome1)[order(table(bs$genome1), decreasing = T)])
  glast = unique(bs$genome2[!bs$genome2 %in% g1])
  g = 1:length(c(g1,glast))
  names(g) = c(g1,glast)

  bs$index1 = NA
  for(i in names(g)) bs$index1[bs$genome1 == i]<-g[i]
  bs$index2 = NA
  for(i in names(g)) bs$index2[bs$genome2 == i]<-g[i]
  bs$abbrev1 = NA
  for(i in names(g)) bs$abbrev1[bs$genome1 == i]<-abbrevs[i]
  bs$abbrev2 = NA
  for(i in names(g)) bs$abbrev2[bs$genome2 == i]<-abbrevs[i]

  bs = bs[order(bs$index1, bs$index2),]

  spl = split(blast.results, paste(blast.results$genome1, blast.results$genome2))

  out = lapply(1:nrow(bs), function(i){
    id = paste(bs[i,1:2], collapse = " ")
    x = spl[[id]]
    gff1 = x[,c("chr1","id1","start1","end1")]
    chr1 =  as.numeric(as.factor(gff1$chr1))
    id1 = x$genome1[1]
    abbrev1 = abbrevs[id1]
    gff1$chr1 = paste0(abbrev1,chr1)

    gff2 = x[,c("chr2","id2","start2","end2")]
    chr2 =  as.numeric(as.factor(gff2$chr2))
    id2 = x$genome2[1]
    abbrev2 = abbrevs[id2]
    gff2$chr2 = paste0(abbrev2,chr2)
    names(gff2) = names(gff1)
    gff = data.table(rbind(gff1, gff2))
    gff = gff[order(gff$chr1, gff$start1),]

    blast = x[,c("id1","id2","perc.iden","align.length",
                 "n.mismatch", "n.gapOpen", "q.start", "q.end",
                 "s.start", "s.end", "eval", "score")]

    return(list(gff = gff, blast = blast))
  })

  gff.o = rbindlist(lapply(out, function(x) x$gff))
  blast.o = rbindlist(lapply(out, function(x) x$blast))

  write.table(gff.o, file = file.path(mcscanx.input.dir,"all.gff"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(blast.o, file = file.path(mcscanx.input.dir, "all.blast"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  mcscan.input = file.path(mcscanx.input.dir,"all")

  if(is.null(MCScanX.params)){
    if(verbose){
      com = paste("MCScanX",mcscan.input)
    }else{
      com = paste("MCScanX",mcscan.input,"&> /dev/null")
    }
  }else{
    if(verbose){
      com = paste("MCScanX",MCScanX.params,mcscan.input)
    }else{
      com = paste("MCScanX",MCScanX.params,mcscan.input,"&> /dev/null")
    }

  }
  system(com)
  mcscan.raw = read.delim(paste0(mcscan.input,".collinearity"),
                          sep = "\t", header  =F,
                          comment.char = "#", strip.white = T,
                          stringsAsFactors = F)

  fac = sapply(as.character(mcscan.raw$V1),function(x) strsplit(x,"-")[[1]][1])
  m = data.table(mcscan.raw[,2:3])
  setnames(m, c("id1","id2"))
  m2 = data.table(m)
  setnames(m, c("id2","id1"))
  m$block.id = as.numeric(as.factor(fac))
  m2$block.id = as.numeric(as.factor(fac))
  out1 = merge(m, blast.results, by = c("id1","id2"))
  out2 = merge(m2, blast.results, by = c("id1","id2"))
  out = rbind(out1, out2)
  out$mapping = paste(out$genome1, out$genome2)
  out$rank1 = frank(out[,c("mapping","chr1","start1")], ties.method = "dense")
  out$rank2 = frank(out[,c("mapping","chr2","start2")], ties.method = "dense")
  if(verbose)
    cat("Done\n")
  return(out)
}
