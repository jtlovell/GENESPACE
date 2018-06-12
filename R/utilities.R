#' @title Utility functions
#' @description
#' \code{utilities} Five utilities functions meant for internal calls in compareGeneSpace
#' @name utilities
#' @param map A `map` object (data.table or data.frame) with the
#' folowing columns: block.id, genome1, genome2, chr1, chr2, start1,
#' start2, end1, end2, rankstart1, rankstart2, rankend1, rankend2.
#' @param gff The path to a gff3-formatted annotation file.
#' @param radii The vector of radii to search within via dbscan
#' @param mappings The number of mappings with a give radius
#' @param id1 First genome ID
#' @param id2 Second genome ID
#' @param col The color to which to add alpha
#' @param alpha see col
#' @param x Numeric or ranks, x values
#' @param y Numeric or ranks, y values
#' @param rerank Logical, if x/y are not ranks, rank them.
#' @param max.jump Numeric, length 1. The maximum euclidean distance between
#' two nodes in the shortest path within a block. If two points are farther
#' apart then this, the block is broken. If a point is isolated by this value,
#' it is dropped.
#' @param tsp.method The method to pass to TSP.solve
#' @param plotit Logical, should plots be made?
#' @param unique A unique identifier for the TSP solver
#' @param fasta.dir Directory containing fasta files to parse
#' @param is.peptide Logical, are fasta files peptide?
#' @param pattern The string identifying the peptide fastas.
#' @param ... Other arguments passed on. Just to solve.TSP for now.
#'
#'
#' @note \code{utilities} is a generic name for the functions documented.
#' \cr
#' If called, \code{utilities} returns its own arguments.
#'
#' @title make_blocks make blocks from mappings
#' @description
#' \code{make_blocks} parses mapping file by block identifier
#' @rdname utilities
#' @import data.table
#' @export
make_blocks<-function(map){
  out.blk = with(map,
                 data.table(
                   block.id = tapply(block.id, block.id, function(x) x[1]),
                   genome1 = tapply(genome1, block.id, function(x) x[1]),
                   genome2 = tapply(genome2, block.id, function(x) x[1]),
                   chr1 = tapply(chr1, block.id, function(x) x[1]),
                   chr2 = tapply(chr2, block.id, function(x) x[1]),
                   start1 = tapply(start1, block.id, min),
                   start2 = tapply(start2, block.id, min),
                   end1 = tapply(end1, block.id, max),
                   end2 = tapply(end2, block.id, max),
                   rankstart1 = tapply(rank1, block.id, min),
                   rankstart2 = tapply(rank2, block.id, min),
                   rankend1 = tapply(rank1, block.id, max),
                   rankend2 = tapply(rank2, block.id, max),
                   stringsAsFactors = F))
  orient = sapply(split(map, map$block.id), function(x) cor(x$start1, x$start2))
  orient = data.frame(block.id = names(orient), orient = ifelse(orient>0,"+","-"))
  out.blk = merge(out.blk, orient, by = "block.id")
  out.blk = out.blk[order(out.blk$genome1, out.blk$chr1, out.blk$start1),]
  map = data.frame(map[order(map$genome1, map$chr1, map$start1),], stringsAsFactors = F)
  blk = data.frame(out.blk, stringsAsFactors = F)
  return(list(block = blk, map = map))
}
#' @title Simple parsing of a gff file
#' @description
#' \code{parse_quickGff} Returns a simplified gff with five columns
#' @rdname utilities
#' @import data.table
#' @export
parse_quickGff = function(gff){
  g = suppressWarnings(
    data.table::fread(gff,showProgress = F, verbose = F))
  g = g[g$V3 == "gene",c(9,1,4,5,7)]
  g$V9 = sapply(g$V9, function(x) gsub("Name=","",strsplit(x,";")[[1]][2]))
  setnames(g, c("id","chr","start","end","strand"))
  return(g)
}
#' @title Subset mappings by density
#' @description
#' \code{loop_dbs} Do dbscan for a varity of radii
#' @rdname utilities
#' @import data.table
#' @import dbscan
#' @export
loop_dbs = function(map, radii, mappings, plotit = T,id1, id2){
  run_dbs = function(map, eps_radius, mappings){
    x=data.frame(map)
    x$x_a = frank(x[,paste0(c("chr_","start_"),id1)],
                  ties.method = "dense")
    x$x_b = frank(x[,paste0(c("chr_","start_"),id2)],
                  ties.method = "dense")

    nn = dbscan::frNN(x[,c("x_a","x_b")], eps = eps_radius)
    dbs = dbscan::dbscan(nn, minPts = mappings)
    return(data.frame(rank1 = x$x_a,
                      rank2 = x$x_b,
                      cluster = dbs$cluster,
                      stringsAsFactors = F))
  }
  if(length(radii)!=length(mappings))
    stop("radii and mappings must be same length\n")
  for(i in 1:length(radii)){
    dclus = run_dbs(map, eps = radii[i], mappings = mappings[i])
    mo = cbind(map, data.table(dclus))
    if(plotit){
      with(mo, plot(rank1, rank2,
                    col = ifelse(cluster == 0,rgb(0,0,0,.2),"darkred"),
                    pch = ".",
                    xlab = paste(id1, "rank"), ylab = paste(id2, "rank"),
                    main = paste("radius =", radii[i], "min mapping =", mappings[i])))
    }
    map = map[mo$cluster !=0,]
  }
  return(map)
}
#' @title Add transparency to a color
#' @description
#' \code{add.alpha} given an alpha value make a color transparent
#' @rdname utilities
#' @export
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
#' @title For a set of xy data, run a TSP solver
#' @description
#' \code{run_TSP} Traveling salesman problem solver pipeline
#' @rdname utilities
#' @import TSP
#' @import data.table
#' @export
run_TSP = function(x,y,
                   unique,
                   rerank = T,
                   tsp.method = "Concorde",
                   max.jump = 5, ...){
  o = data.frame(x,y)
  if(rerank){
    o[,1]<-frank(o[,1], ties.method = "dense")
    o[,2]<-frank(o[,2], ties.method = "dense")
  }
  rownames(o)<-unique
  etsp = TSP::TSP(as.ETSP(o), labels = unique)
  etsp = TSP::insert_dummy(etsp, n = 1, label = "cut")
  tour <- TSP::solve_TSP(etsp, method = tsp.method, control = list(verbose = F,clo = "-V"))
  tour <- TSP::cut_tour(tour, "cut")


  o = o[labels(tour)[order(tour)],]
  o$x1 = c(o$x[-1], o$x[length(o$x)])
  o$y1 = c(o$y[-1], o$y[length(o$y)])
  o$x2 = c(o$x[1], o$x[-length(o$x)])
  o$y2 = c(o$y[1], o$y[-length(o$y)])
  o$euc1 = with(o, sqrt(abs(x-x1)^2 + abs(y-y1)^2))
  o$euc2 = with(o, sqrt(abs(x-x2)^2 + abs(y-y2)^2))

  o$subclus = cumsum(o$euc2>max.jump)+1
  o$subclus[abs(o$euc1)>max.jump & abs(o$euc2)>max.jump]<-NA

  o = o[unique,]
  return(o$subclus)
}
#' @title Parse complex fasta header
#' @description
#' \code{parse_fastaHeader} Subset space-delimited fasta header, retaining the "locus" entry
#' @rdname utilities
#' @import Biostrings
#' @export
parse_fastaHeader = function(fasta.dir, is.peptide = T,
                             pattern = "pep.fa", verbose = T){

  in.pep = list.files(peptide.dir,
                      pattern = "pep.fa",
                      full.names = T)

  if(verbose)
    cat("Renaming fasta headers ...\n")
  ss = lapply(in.trs, function(i){
    if(verbose)
      cat("...",i,"\n\t")
    if(is.peptide){
      x = Biostrings::readAAStringSet(i)
    }else{
      x = Biostrings::readDNAStringSet(i)
    }
    if(verbose)
      cat("original names (e.g.):", names(x)[1])
    names(x)<-sapply(gsub(".*locus=","",names(x)),
                     function(y) strsplit(y," ")[[1]][1])
    if(verbose)
      cat("\n\tparsed names (e.g.):", names(x)[1],"\n")
    writeXStringSet(x, filepath = i)
  })
}


#' @title Split blocks based on density, iteratively
#' @description
#' \code{split_byDensity} takes a mapping object and splits it based on rules.
#' @rdname utilities
#' @import dbscan
#' @export
split_byDensity<-function(map, max.jump = 5, min.dens = max.jump-1,
                          zScore2split = 5){

  split_it = function(map, radius,max.jump,min.dens,zScore2split){
    spl = split(map, map$block.id)
    out = do.call(rbind, lapply(names(spl), function(i){
      x = spl[[i]]
      d1 = diff(x$x_a)
      d2 = diff(x$x_b)
      z1 = abs((d1-mean(d1))/(sd(d1)+.000001))
      z2 = abs((d2-mean(d2))/(sd(d2)+.000001))
      x$was_split = F
      if(nrow(x)>10 & max(c(d1,d2))>max.jump & max(c(z1,z2))>zScore2split){
        nn = dbscan::frNN(x[,c("x_a","x_b")], eps = radius)
        dbs = dbscan::dbscan(nn, minPts = min.dens)

        x$newclust = dbs$cluster
        x = x[x$newclust!=0,]
        x$block.id = paste0(x$block.id, "_",x$newclust)
        x$was_split = length(unique(dbs$cluster))>1
        x$newclust = NULL
      }
      return(x)
    }))
  }

  map$x_a = frank(map[,c("chr1","start1")],
                  ties.method = "dense")
  map$x_b = frank(map[,c("chr2","start2")],
                  ties.method = "dense")

  ntospl = 1
  map$was_split = TRUE
  while(any(map$was_split)){
    if(verbose)
      cat("checked", sum(map$was_split[!duplicated(map$block.id)]),"blocks\n")
    map = split_it(map,
                   max.jump= max.jump,
                   min.dens = min.dens,
                   radius = sqrt((max.jump^2)*2),
                   zScore2split = zScore2split)
  }
  out = make_blocks(map)
  if(verbose)
    cat("split into",nrow(out$block),"... Done")
  return(list(map = out$map, block = out$block))
}




#' @title Plot mapping and blocks
#' @description
#' \code{plot_blocksAndMapping} takes a mapping and block object and plots selected genomes / chromosomes.
#' @rdname utilities
#' @export
plot_blocksAndMapping = function(map,
                                 blk,
                                 ref.id,
                                 altGenome2plot,
                                 chr1toplot,
                                 chr2toplot,
                                 main = NULL){
  tpb = blk[blk$genome1 == ref.id & blk$genome2 == altGenome2plot,]
  tp = map[map$genome1 == ref.id & map$genome2 == altGenome2plot,]
  tpb$s1 = with(tpb, ifelse(orient == "+", rankstart1, rankend1))
  tpb$e1 = with(tpb, ifelse(orient == "+", rankend1, rankstart1))
  tpb$rankstart1 = tpb$s1
  tpb$rankend1 = tpb$e1
  tpb$s1 = NULL
  tpb$e1 = NULL
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
                    ylab = paste(altGenome2plot,j,"gene order"),
                    xlab = paste(ref.id,i,"gene order"),
                    main = main))
      with(b2, segments(x0 = rankstart1, x1 = rankend1,
                        y0 = rankstart2, y1 =rankend2,
                        col = "black", lwd = 1.5))
      with(b2, text(x = rowMeans(b2[,c("rankstart1","rankend1")]),
                    y = rowMeans(b2[,c("rankstart2","rankend2")]),
                    labels = block.id,
                    col = "black", cex = .5, adj = c(1,-1)))
    }
  }
}



