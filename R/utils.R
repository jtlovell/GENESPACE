#' @title Generic internal functions used by genespace
#' @description
#' \code{utils} Convience functions for genespace, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param gsParam list of genespace parameters
#' @param path file path character string
#' @param pattern regular expression to search for
#' @param recursive logical, should the search be recursive
#' @param which character specifying which method to use
#' @param x vector of observations
#' @param col a vector of colors
#' @param id1 character vector of ids with length matching id2
#' @param id2 character vector of ids with length matching id1
#' @param min numeric, length 1 specifying the minumum value in the scaled data
#' @param max numeric, length 1 specifying the maximum value in the scaled data
#' @param dt data.table with the first two columns as id1, id2.
#' @param raw logical, should a raw vector of peptide widths be returned?
#' @param scale1toMean logical, if single value, should it be the mean?
#' @param alpha numeric 0-1, specifying transparency
#' @param width integer, the number of characters to return in string.
#' @param x x position of the scale bar
#' @param y y position of the scale bar
#' @param xspan amount of span on the x axis
#' @param yspan amount of span on the y axis
#' @param label scale bar label
#' @param cex scale bar label character expansion
#' @param lwd line thickness for scale bar
#' @param xleft numeric, specifying the coordinate of left x position
#' @param ybottom numeric, specifying the coordinate of lower y position
#' @param xright numeric, specifying the coordinate of right x position
#' @param ytop numeric, specifying the coordinate of upper y position
#' @param start1 numeric, specifying the coordinate of blk start in genome1
#' @param end1 numeric, specifying the coordinate of  blk end in genome1
#' @param start2 numeric, specifying the coordinate of  blk start in genome2
#' @param end2 numeric, specifying the coordinate of blk end in genome2
#' @param y1 numeric, specifying the coordinate of y position in genome1
#' @param y2 numeric, specifying the coordinate of  y position in genome1
#' @param to what should the value be rounded to?
#' @param byGrpCol column which serves as the by factor
#' @param windowSize integer specifying the window size
#' @param yCol character specifying the column name containing the y values
#' @param fun function to pass to sliding window
#' @param gffFiles vector of file paths pointing to the gff files
#' @param refGenome character string specifying which genome is the reference
#' @param nCores integer, the number of parallel processes to run
#' @param blFile file path to the blast-like text file
#' @param ofID1 orthofinder ID of the first gene
#' @param ofID2 orthofinder ID of the second gene
#' @param onlyIDScore logical, should only the geneIDs and score be returned?
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'
#' @title Order files by date of last modification
#' @description
#' \code{order_filesByMtime} Order files by date of last modification
#' @rdname utils
#' @export
order_filesByMtime <- function(path = getwd(),
                               pattern = "*",
                               recursive = F){
  if (length(path) == 1) {
    allFiles <- list.files(
      path = path,
      full.names = T,
      pattern = pattern,
      recursive = recursive)
  }else{
    allFiles <- path
  }

  details <- file.info(allFiles, extra_cols = F)
  details <- details[rev(with(details, order(as.POSIXct(mtime)))), ]
  return(rownames(details))
}

#' @title Check if orthofinder is installed
#' @description
#' \code{check_orthofinderInstall} Check if orthofinder is installed
#' @rdname utils
#' @export
check_orthofinderInstall <- function(path){
  grepl("OrthoFinder",
        system(paste(path, "-h"),
               intern = T)[2])
}

#' @title Check logical argument
#' @description
#' \code{check_logicalArg} Ensure a logical (length 1) argument is coerced
#' correctly. If cannot be coerced, return error.
#' @rdname utils
#' @importFrom parallel detectCores
#' @export
check_logicalArg <- function(x){
  x <- as.logical(x)
  if(is.null(x)) stop(x, "must be logical\n")
  if(is.na(x)) stop(x, "must be logical\n")
  if(length(x) == 0) stop(x, "must be logical\n")
  if(length(x) > 1) x <- x[1]
  return(x)
}

#' @title convert vector to RLE
#' @description
#' \code{add_rle} Convert vector into run length equivalents
#' @rdname utils
#' @export
add_rle <- function(x, which = "n"){
  if (which == "n") {
    rep(rle(x)$lengths, rle(x)$lengths)
  }else{
    rep(1:length(rle(x)$lengths), rle(x)$lengths)
  }
}

#' @title clus_igraph
#' @description
#' \code{clus_igraph} Clus_igraph
#' @rdname utils
#' @import data.table
#' @importFrom igraph graph_from_data_frame clusters
#' @export
clus_igraph <- function(id1, id2){
  if(length(unique(id1)) == 1 & length(unique(id2)) == 2){
    return(rep(1, length(id1)))
  }else{
    return(clusters(graph_from_data_frame(
      data.frame(id1, id2),
      directed = F))$membership[id1])
  }
}

#' @title count the number of amino acids by gene
#' @description
#' \code{get_nAA} count the number of amino acids by gene
#' @rdname utils
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
get_nAA <- function(path, raw = FALSE){
  if(!raw){
    pepF <- list.files(path, pattern = "^Species", full.names = T)
    pepF <- pepF[grep(".fa$", pepF)]
    peps <- rbindlist(lapply(pepF, function(x){
      y <- readAAStringSet(x)
      return(data.table(ofID = names(y),
                        nAA = width(y)))
    }))
  }else{
    y <- readAAStringSet(path)
    o <- width(y)
    names(o) <- names(y)
    return(o)
  }
}

#' @title calculate the mode
#' @description
#' \code{find_modalValue} find the most commmon value in a series
#' @rdname utils
#' @export
find_modalValue <- function(x){
  tab <- table(x)
  return(names(tab[order(-tab)])[1])
}

#' @title scale a vector between a range
#' @description
#' \code{scale_between} scale a vector between a range
#' @rdname utils
#' @export
scale_between <- function(x, min, max){
  (x - min(x)) / (max(x) - min(x)) * (max - min) + min
}

#' @title complete subgraphs
#' @description
#' \code{complete_subgraph} given a hits data.table, expand unique orthofinder
#' IDs among all unique elements in a subgraph
#' @rdname utils
#' @export
complete_subgraph <- function(dt){
  clusters <- clusters(graph_from_data_frame(
    dt,
    directed = F))
  clusters <- with(clusters, data.table(
    id = names(membership),
    group = membership))
  out <- data.table(merge(
    clusters,
    clusters,
    by = "group",
    allow.cartesian = T))
  setnames(out, c("cluster", "x", "y"))
  return(out)
}

#' @title flatten a list
#' @description
#' \code{flatten_list} convert a list into a vector while keeping names
#' @rdname utils
#' @export
flatten_list <- function(x){
  y <- unlist(x)
  names(y) <- rep(names(x), sapply(x, length))
  return(y)
}


#' @title checks parameter for minimum length of a peptide
#' @description
#' \code{check_minPepLen} checks for init_genespace
#' @rdname utils
#' @export
check_minPepLen <- function(x){
  x <- as.integer(x)[1]
  if(is.na(x) || is.null(x))
    stop("minPepLen must be an integer\n")
  if(x < 0)
    x <- 0
  return(x)
}

#' @title checks parameter for interleave size between blocks
#' @description
#' \code{check_dropInterSize} checks for init_genespace
#' @rdname utils
#' @export
check_dropInterSize <- function(x){
  x <- as.integer(x)[1]
  if(is.na(x) || is.null(x))
    stop("dropInterleavesSmallerThan must be an integer\n")
  if(x < 1)
    x <- 1
  return(x)
}

#' @title check if a vector is coercible to R colors
#' @description
#' \code{are_colors} check if a vector is coercible to R colors
#' @rdname utils
#' @importFrom grDevices col2rgb
#' @export
are_colors <- function(col) {
  sapply(col, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}

#' @title add transparency to a color
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col,
                      alpha = 1){

  if (missing(col))
    stop("Please provide a vector of colors.")
  if (!all(are_colors(col)))
    stop("Please provide a vector of colors.")

  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1],
              x[2],
              x[3],
              alpha = alpha))
}

#' @title
#' @description
#' \code{scale_between} ...
#' @rdname utils
#' @export
scale_between <- function(x, min, max, scale1toMean = TRUE){
  if(length(unique(x)) > 1){
    return((x - min(x)) / (max(x) - min(x)) * (max - min) + min)
  }else{
    if(scale1toMean){
      return(mean(c(min, max)))
    }else{
      return(max)
    }
  }
}

#' @title pull_strWidth
#' @description
#' \code{pull_strWidth} pull_strWidth
#' @rdname utils
#' @export
pull_strWidth <- function(x, width){
  y <- substr(x, 1, width)
  if(nchar(y) < width)
    y <- sprintf("%s%s", y, paste(rep(" ", width-nchar(y)), collapse = ""))
  return(y)
}

#' @title round_toInteger
#' @description
#' \code{round_toInteger} round_toInteger
#' @rdname utils
#' @export
round_toInteger <- function(x, to){
  round(x / to, 0) * to
}

#' @title calc_partialSw
#' @description
#' \code{calc_partialSw} calc_partialSw
#' @rdname utils
#' @export
calc_partialSw <- function(dt,
                           byGrpCol,
                           windowSize,
                           yCol,
                           fun = function(x) mean(x, na.rm = T)){
  yDataFakeCol <- byGrpFakeCol <- swLeft <- swRight <- NULL
  dt <- data.table(dt)
  # rename variables
  if("yDataFakeCol" %in% colnames(dt) & yCol != "yDataFakeCol")
    stop("yDataFakeCol in colnames but not ycol, rename please\n")
  if("byGrpFakeCol" %in% colnames(dt) & byGrpCol != "byGrpFakeCol")
    stop("byGrpFakeCol in colnames but not byGrpCol, rename please\n")

  if(yCol != "yDataFakeCol")
    dt[,yDataFakeCol := dt[[yCol]]]
  if(byGrpCol != "byGrpFakeCol")
    dt[,byGrpFakeCol := dt[[byGrpCol]]]

  # apply sliding windows
  dt[,`:=`(swLeft = frollapply(yDataFakeCol, n = windowSize, align = "left", FUN = fun),
           swRight = frollapply(yDataFakeCol, n = windowSize, align = "right", FUN = fun)),
     by = "byGrpFakeCol"]

  # replace NAs with the alternative means
  dt$swLeft[is.na(dt$swLeft)] <- dt$swRight[is.na(dt$swLeft)]
  dt$swRight[is.na(dt$swRight)] <- dt$swLeft[is.na(dt$swRight)]

  # get mean of left and right aligned windows
  dt[,`:=`(sw = (swLeft + swRight)/2,
           swLeft = NULL, swRight = NULL)]

  # drop renamed variables
  if(yCol != "yDataFakeCol")
    dt[,yDataFakeCol := NULL]
  if(byGrpCol != "byGrpFakeCol")
    dt[,byGrpFakeCol := NULL]
  return(dt)
}

#' @title load gff annotations
#' @description
#' \code{read_gff} reads a genespace-formatted gff-like annotation file into
#' memory
#' @rdname utils
#' @export
read_gff <- function(gffFiles){
  genome <- NULL
  gff <- rbindlist(lapply(names(gffFiles), function(i){
    x <- fread(
      gffFiles[[i]],
      key = c("chr", "start","end","strand"))
    x[,genome := i]
    return(x)
  }))
  gffCols <- c("ord","chr","genome","start", "end", "strand","id")
  if(!all(gffCols %in% colnames(gff)))
    stop(paste(gffCols, collapse = ", "),
         " must all be column names in gff\n")
  return(gff)
}

#' @title Read orthofinder species IDs
#' @description
#' \code{read_orthofinderSpeciesIDs} Parses the SpeciesIDs.txt file into a
#' data.table and returns to R.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSpeciesIDs <- function(path){
  genome <- NULL
  si <- fread(
    file.path(path, "SpeciesIDs.txt"),
    sep = ":",
    header = F,
    col.names = c("genomeNum", "genome"),
    colClasses = c("numeric", "character"))
  si[,genome := gsub(".fa", "", genome, fixed = T)]
  sio <- si$genomeNum; names(sio) <- si$genome
  return(sio)
}

#' @title Read orthofinder sequence IDs
#' @description
#' \code{read_orthofinderSequenceIDs} Reads the sequence
#' IDs:gene name dictionary into memory.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSequenceIDs <- function(path){
  ofID <- NULL
  gi <- fread(
    file.path(path, "SequenceIDs.txt"),
    header = F,
    sep = ":",
    col.names = c("ofID","id"),
    colClasses = c("character","character"))
  gi[,c("genomeNum","geneNum") := tstrsplit(ofID, "_", type.convert = T)]
  return(gi)
}

#' @title choose_mostRecentOF
#' @description
#' \code{choose_mostRecentOF} choose_mostRecentOF
#' @rdname utils
#' @import data.table
#' @export
choose_mostRecentOF <- function(path){
  p <- file.path(path, "Orthofinder")
  ps <- order_filesByMtime(p)
  pschk <- lapply(ps, function(x){
    chks <- c(dir.exists(file.path(x, "Orthologues")),
              dir.exists(file.path(x, "Gene_Duplication_Events")),
              file.exists(file.path(x, "Gene_Duplication_Events","Duplications.tsv")))
    if(all(chks))
      return(x)
  })
  wh <- min(which(!sapply(pschk, is.null)))
  if(length(wh) > 0){
    return(pschk[[wh]])
  }else{
    stop("Cannot find orthogues and duplications, has the full orthofinder pipe been run?\n")
  }
}

#' @title parse orthogroups file into a data.table
#' @description
#' \code{parse_ogs} wide to long format conversion for orthogroups.tsv
#' @rdname utils
#' @export
parse_ogs <- function(gsParam){
  id <- genome <- Orthogroup <- NULL
  ogtsv <- file.path(gsParam$paths$orthogroupsDir, "Orthogroups.tsv")
  og <- fread(ogtsv, showProgress = F, verbose = F)
  og <- melt(
    og, id.vars = "Orthogroup", variable.name = "genome", value.name = "id")
  og <- og[,list(id = strsplit(id, ",")[[1]]), by = c("Orthogroup", "genome")]
  og[,`:=`(genome = trimws(genome), id = trimws(id), Orthogroup = trimws(Orthogroup))]
  setnames(og, 1, "ogID")
  return(og)
}

#' @title parse_orthologues
#' @description
#' \code{parse_orthologues} parse_orthologues
#' @rdname utils
#' @import data.table
#' @importFrom parallel mclapply
#' @export
parse_orthologues <- function(gsParam, refGenome, nCores = 1){
  id1 <- id2 <- orthID <- NULL
  od <- file.path(gsParam$paths$orthologuesDir,
                  sprintf("Orthologues_%s", refGenome))
  odf <- list.files(od, full.names = T, pattern = "__v__")
  ogo <- rbindlist(mclapply(odf, mc.cores = nCores, function(i){
    x <- fread(i, showProgress = F)
    refID <- colnames(x)[2]; altID <- colnames(x)[3]
    setnames(x, c("og", "id1", "id2"))
    x1 <- subset(x, !grepl(",", paste(id1, id2)))
    x2 <- subset(x, grepl(",", paste(id1, id2)))
    x2[,orthID := 1:.N]
    x2r <- x2[,list(id1 = unique(strsplit(id1, ",")[[1]])),
              by = "orthID"]
    x2a <- x2[,list(id2 = unique(strsplit(id2, ",")[[1]])),
              by = "orthID"]
    x2 <- merge(x2r, x2a, by = "orthID", all = T, allow.cartesian = T)
    x1[,orthID := (1:.N)+max(x2$orthID)]
    x <- rbind(x1[,colnames(x2), with = F], x2)
    x[,`:=`(gen1 = refID, gen2 = altID,
            id1 = gsub(" ", "", id1), id2 = gsub(" ", "", id2))]
    return(x)
  }))
  return(ogo)
}

#' @title Read orthofinder blast file
#' @description
#' \code{read_blast} Reads in a single pairwise orthofinder-formaatted blast
#' file
#' @rdname utils
#' @import data.table
#' @export
read_blast <- function(blFile = NULL,
                       ofID1 = NULL,
                       ofID2 = NULL,
                       path = NULL,
                       onlyIDScore = TRUE){
  V12 <- score <- NULL
  if(is.null(blFile)){
    blFile <- file.path(path, paste0("Blast", ofID1, "_", ofID2,".txt.gz"))
  }
  if(!file.exists(blFile))
    stop("cannot find ", blFile, "\n")

  if(!onlyIDScore){
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE)
    g1 <- strsplit(bl$V1[1], "_")[[1]][1]
    g2 <- strsplit(bl$V2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[,c(2,1,3:6,8,7,10,9,11,12)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl),with = F]
      bl <- rbind(bl, tmp)
      setorder(bl, -V12)
      bl <- subset(bl, !duplicated(bl[,c(1:2)]))
    }
  }else{
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE,
      select = c(1,2,12),
      col.names = c("ofID1","ofID2","score"))
    g1 <- strsplit(bl$ofID1[1], "_")[[1]][1]
    g2 <- strsplit(bl$ofID2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[,c(2,1,3)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl),with = F]
      bl <- rbind(bl, tmp)
      setorder(bl, -score)
      bl <- subset(bl, !duplicated(bl[,c(1:2)]))
    }
  }

  return(bl)
}

#' @title draw_scaleBar
#' @description
#' \code{draw_scaleBar} draw_scaleBar
#' @rdname utils
#' @export
draw_scaleBar <- function(x, y, yspan, xspan, label, lwd, cex){
  xstart <- x - (xspan / 2)
  xend <- x + (xspan / 2)
  ytop <- y + (yspan / 2)
  ybottom <- y - (yspan / 2)
  segments(x0 = xstart, x1 = xend, y0 = y, y1 = y, lwd = lwd)
  segments(x0 = xstart, x1 = xstart, y0 = ytop, y1 = ybottom, lwd = lwd)
  segments(x0 = xend, x1 = xend, y0 = ytop, y1 = ybottom, lwd = lwd)
  text(x = xstart + (xspan / 2), y = ytop + (yspan / 2), labels = label, adj = c(.5,0), cex = cex)
}


#' @title cosine curve source data
#' @description
#' \code{cosine_points} vector of points for polygons based on cosine curves
#' @rdname utils
#' @export
cosine_points <- function(){
  npts = 1e4 # initial number of points
  keepat = round(npts / 20) # grid to keep always
  grid <- seq(from = 0, to = pi, length.out = npts) # grid
  x <- (1 - cos(grid)) / max((1 - cos(grid))) # scaled cosine
  y <- grid / max(grid) # scaled grid
  # calculate slope for each point
  x1 <- x[-1];  y1 <- y[-1]
  x2 <- x[-length(x)];  y2 <- y[-length(y)]
  s <-  (y1 - y2) / (x1 - x2)
  # choose points that capture changes in slope
  ds <- cumsum(abs(diff(s)))*5
  wh <- c(1,which(!duplicated(round(ds))), length(x))
  wh2 <- c(wh, seq(from = 0, to = length(x), by = round(keepat)))
  wh <- c(wh, wh2)[!duplicated(c(wh, wh2))]
  wh <- wh[order(wh)]
  return(cbind(x[wh], y[wh]))
}

#' @title convert cosine points to polygon
#' @description
#' \code{calc_curvePolygon} from 2d coordinates, make a curve
#' @rdname utils
#' @export
calc_curvePolygon <- function(start1,
                              end1 = NULL,
                              start2,
                              end2 = NULL,
                              y1,
                              y2){
  scaledCurve <- cosine_points()
  if (!is.null(end1) | !is.null(end2)) {
    tp <- rbind(
      start1 = data.table(
        x = start1, y = y1),
      poly1 = data.table(
        x = scale_between(x = scaledCurve[,1], min = start1, max = start2),
        y = scale_between(x = scaledCurve[,2], min = y1, max = y2)),
      start2 = data.table(x = start2, y = y2),
      end2 = data.table(
        x = end2, y = y2),
      poly2 = data.table(
        x = scale_between(x = scaledCurve[,1], min = end2, max = end1),
        y = scale_between(x = scaledCurve[,2], min = y2, max = y1)),
      end1 = data.table(
        x = end1, y = y1))
  }else{
    tp <- data.table(
      x = scale_between(x = scaledCurve[,1], min = start1, max = start2),
      y = scale_between(x = scaledCurve[,2], min = y1, max = y2))
  }
  return(tp)
}

#' @title calculate coordinates for rounded rectange polygons
#' @description
#' \code{round_rect} from x-y coordinates, make a rounded rectangle
#' @rdname utils
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
round_rect <- function(xleft, ybottom, xright, ytop){

  if (ytop <= ybottom)
    stop("ytop must be > ybottom")
  if (xleft >= xright)
    stop("xleft must be < xright")

  # measure graphics device
  asp <- diff(par("usr")[3:4]) / diff(par("usr")[1:2])
  dev <- dev.size()[1] / dev.size()[2]

  # make a curve and split into left and right
  radius <- (ytop - ybottom) / 2
  centerY <- ytop - radius
  centerX <- mean(c(xleft, xright))
  theta <- seq(0, 2 * pi, length = 200)
  circX <- cos(theta)
  circY <- sin(theta)
  leftC <- which(circX <= 0)
  rightC <- which(circX >= 0)

  xR <- circX[rightC]
  yR <- circY[rightC]
  ordYR <- rev(order(yR))
  xR <- xR[ordYR]
  yR <- yR[ordYR]

  xL <- circX[leftC]
  yL <- circY[leftC]
  ordYL <- order(yL)
  xL <- xL[ordYL]
  yL <- yL[ordYL]

  # project onto graphics device and scale
  xRightS <- xright - (radius / asp / dev)
  xLeftS <- xleft + (radius / asp / dev)
  if (centerX < xLeftS)
    xLeftS <- centerX
  if (centerX > xRightS)
    xRightS <- centerX
  xLS <- scale_between(xL, xleft, xLeftS)
  xRS <- scale_between(xR, xRightS, xright)
  yLS <- scale_between(yL, ybottom, ytop)
  yRS <- scale_between(yR, ybottom, ytop)
  return(data.table(x = c(xRS,xLS), y = c(yRS, yLS)))
}

