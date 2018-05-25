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
  map = data.frame(mapo, stringsAsFactors = F)
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
  tour <- TSP::solve_TSP(etsp, method = tsp.method)
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


