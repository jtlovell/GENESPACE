#' @title Make input metadata for pipe_Diamond2MCScanX
#'
#' @description
#' \code{make_inputFileMatrix} Utility function to build metadata
#'
#' @param ref.id Character, the reference genome identifier
#' @param MCScanX.path The path to the MCScanX program
#' @param mcs_mapping.dir Directory containing the MCScanX-formatted mapping files
#' @param mapping.obj The object returned from the mapping pipeline
#' @param MCScanX.params Additionally parameters to pass to MCScanX
#' @param buffer Buffer to find overlapping blocks (see merge_overlappingBlocks.R)
#' @param chr1toplot The chromosome ids in the reference genome to plot
#' @param chr2toplot The chromosome ids in the altGenome2plot to plot
#' @param altGenome2plot The identifier of the alternative genome to plot
#' @param plotit Logical, hould a plot be made
#' @param verbose Logical, should updates be printed?
#' @param ... Not currently in use
#' @details See pipe_Diamond2MCScanX for more information.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
run_MCScanX<-function(ref.id, MCScanX.path, mcs_mapping.dir, mapping.obj,
                      MCScanX.params = NULL, verbose = T, plotit = T,
                      chr1toplot = NULL, chr2toplot = NULL,
                      altGenome2plot = NULL, buffer = 4){

  if(is.null(altGenome2plot) & plotit)
    stop("If plotit = TRUE, must provide alternate genomeID\n")
  if(verbose)
    cat("concatenating gff files\n")
  ingff=file.path(mcs_mapping.dir,"*gff")
  outgff = file.path(mcs_mapping.dir,
                     paste0(ref.id,"_all.gff"))
  system(paste("cat", ingff,">", outgff))

  if(verbose)
    cat("concatenating blast results\n")

  inblast=file.path(mcs_mapping.dir,"*blast")
  outblast = file.path(mcs_mapping.dir,
                       paste0(ref.id,"_all.blast"))
  system(paste("cat", inblast,">", outblast))

  mcscan.input = file.path(mcscan.dir, paste0(ref.id,"_all"))
  if(is.null(MCScanX.params)){
    com = paste(file.path(MCScanX.path,"MCScanX"),mcscan.input)
  }else{
    com = paste(file.path(MCScanX.path,"MCScanX"),MCScanX.params,mcscan.input)
  }
  if(verbose)
    cat("running MCScanX via:\n\t", com,"\n")

  system(paste(com))
  system(paste("rm",outgff))
  system(paste("rm", outblast))

  if(verbose)
    cat("Merging output\n")
  init.mcscan2 = read.delim(paste0(mcscan.input,".collinearity"), sep = "\t", header  =F,
                     comment.char = "#", strip.white = T, stringsAsFactors = F)

  gff.list = sapply(unique(c(id.mat$gff1, id.mat$gff2)), simplify = F, USE.NAMES = T, parse_quickGff)
  names(gff.list)<-gsub(".gff3","",gsub(".*gff/","",names(gff.list)))

  gff.index = rbindlist(lapply(names(gff.list), function(x) data.table(genome = x, geneID = gff.list[[x]][,1])))
  gff.index1 = gff.index
  gff.index2 = gff.index


  setnames(init.mcscan2,c("info","geneID1","geneID2","score"))
  setnames(allmap2,c("geneID2","chr2","start2","end2","strand2"))
  mco = merge(allmap2, init.mcscan2, by = "geneID2")
  setnames(allmap1,c("geneID1","chr1","start1","end1","strand1"))
  mco = merge(allmap1,mco, by = "geneID1")
  setnames(gff.index2, c("genome2","geneID2"))
  mco = merge(gff.index2,mco, by = "geneID2")
  setnames(gff.index1, c("genome1","geneID1"))
  mco = merge(gff.index1,mco, by = "geneID1")

  setkey(mco, genome1, genome2, chr1, chr2)
  mco$block.id = gsub("-.*","",mco$info)

  out.map = mco
  out.map$rank1 = frank(out.map[,c("genome1","chr1","start1")], ties.method = "dense")
  out.map$rank2 = frank(out.map[,c("genome2","start2")], ties.method = "dense")

  out.blk = make_blocks(out.map)
  map = data.frame(out.blk[["map"]], stringsAsFactors = F)
  blk = data.frame(out.blk[["block"]], stringsAsFactors = F)

  if(plotit){

    if(verbose)
      cat("Genereating plots\n")

    tpb = blk[blk$genome1 == ref.id & blk$genome2 == altGenome2plot,]
    tp = map[map$genome1 == ref.id & map$genome2 == altGenome2plot,]
    if(!is.null(chr1toplot)){
      tpb = tpb[tpb$chr1 %in% chr1toplot,]
      tp = tp[tp$chr1 %in% chr1toplot,]
    }
    if(!is.null(chr2toplot)){
      tpb = tpb[tpb$chr2 %in% chr2toplot,]
      tp = tp[tp$chr2 %in% chr2toplot,]
    }
    sb1 = split(tpb, tpb$chr1)
    st1 = split(tp, tp$chr1)
    for(i in names(sb1)){
      sb2 = split(sb1[[i]], sb1[[i]]$chr2)
      st2 = split(st1[[i]], st1[[i]]$chr2)
      for(j in names(sb2)){
        t2 = st2[[j]]
        b2 = sb2[[j]]
        cols = rep_len(c("red3","salmon","darkorange3","gold",
                         "grey50","lightgreen","forestgreen","darkgreen",
                         "cyan","dodgerblue3","violet","purple"), length.out = nrow(blk))
        with(t2, plot(rank1, rank2,
                      col = cols[as.numeric(as.factor(block.id))],
                      pch = 16, cex = .5,
                      main = paste("Raw MCScanX output:",ref.id,
                                   "vs.",altGenome2plot,"... n.genes = ", nrow(t2)),
                      ylab = paste(altGenome2plot,j,"gene order"),
                      xlab = paste(ref.id,i,"gene order")))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 =rankend2,
                          col = "black", lwd = 1.5))
      }
    }
  }
  if(verbose)
    cat("MCScanX - Done!\n")

  merged = merge_overlappingBlocks(map = out.map, blk = out.blk,
                                   verbose = verbose, buffer = buffer)
  out.map = merged$map
  out.blk = merged$block

  blk = out.blk
  map = out.map
  if(plotit){

    if(verbose)
      cat("Genereating plots\n")

    tpb = blk[blk$genome1 == ref.id & blk$genome2 == altGenome2plot,]
    tp = map[map$genome1 == ref.id & map$genome2 == altGenome2plot,]
    if(!is.null(chr1toplot)){
      tpb = tpb[tpb$chr1 %in% chr1toplot,]
      tp = tp[tp$chr1 %in% chr1toplot,]
    }
    if(!is.null(chr2toplot)){
      tpb = tpb[tpb$chr2 %in% chr2toplot,]
      tp = tp[tp$chr2 %in% chr2toplot,]
    }
    sb1 = split(tpb, tpb$chr1)
    st1 = split(tp, tp$chr1)
    for(i in names(sb1)){
      sb2 = split(sb1[[i]], sb1[[i]]$chr2)
      st2 = split(st1[[i]], st1[[i]]$chr2)
      for(j in names(sb2)){
        t2 = st2[[j]]
        b2 = sb2[[j]]
        cols = rep_len(c("red3","salmon","darkorange3","gold",
                         "grey50","lightgreen","forestgreen","darkgreen",
                         "cyan","dodgerblue3","violet","purple"), length.out = nrow(b2))
        with(t2, plot(rank1, rank2,
                      col = cols[as.numeric(as.factor(block.id))],
                      pch = 16, cex = .5,
                      main = paste("Merged MCScanX output:",ref.id,
                                   "vs.",altGenome2plot,"... n.genes = ", nrow(t2)),
                      ylab = paste(altGenome2plot,j,"gene order"),
                      xlab = paste(ref.id,i,"gene order")))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 =rankend2,
                          col = "black", lwd = 1.5))
      }
    }
  }

  return(list(block = out.blk, map = out.map))
}
