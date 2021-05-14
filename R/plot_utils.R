#' @title Generic internal functions used by genespace
#' @description
#' \code{plot_utils} Convience functions for genespace, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name plot_utils
#'
#' @param col vector of length 1, coercible to an R color
#' @param alpha numeric  0<=x<=1 of length 1, specifying transparency
#' @note \code{plot_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{plot_utils} returns its own arguments.
#'
#'
#' @title check if a vector is coercible to R colors
#' @description
#' \code{areColors} check if a vector is coercible to R colors
#' @rdname plot_utils
#' @export
areColors <- function(col) {
  sapply(col, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}

#' @title add transparency to a color
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname plot_utils
#' @export
add_alpha <- function(col,
                      alpha = 1){

  if (missing(col))
    stop("Please provide a vector of colors.")
  if (!all(areColors(col)))
    stop("Please provide a vector of colors.")

  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1],
              x[2],
              x[3],
              alpha = alpha))
}

#' #' @title calc_linearGenePos
#' #' @description
#' #' \code{calc_linearGenePos} calc_linearGenePos
#' #' @rdname plot_utils
#' #' @export
#' calc_linearGenePos <- function(gsAnnot,
#'                                chrList,
#'                                genomeIDs,
#'                                gapProp = .01,
#'                                useOrder = TRUE){
#'   gff <- add_ofID2gff(read_gff(gsAnnot$gff[genomeIDs]), gsParam$blast)
#'   spl <- split(gff, by = "genome")
#'   xv <- lapply(spl, function(x){
#'     x[,chrn := as.numeric(factor(chr, levels = chrl[[x$genome[1]]]))]
#'     x <- subset(x, !is.na(chrn))
#'     setkey(x, chrn, start, end)
#'     if(useOrder){
#'       gapSize <- ceiling(nrow(x) * gapProp)
#'       x[,offset := gapSize * (chrn-1)]
#'       x[,xpos := (1:.N)+offset]
#'     }else{
#'       chroff <- x[,list(offset = max(end)), by = "chr"]
#'       chroff[,offset := cumsum(offset)]
#'       gapSize <- ceiling(max(chroff$offset) * gapProp)
#'       gapSize <- c(0, cumsum(rep(gapSize, nrow(chroff)-1)))
#'       chroff[,offset := c(0, offset[1:(.N-1)]) + gapSize]
#'       cov <- chroff$offset; names(cov) <- chroff$chr
#'       x[,offset := cov[chr]]
#'       x[,xpos := ((start + end)/2) + offset]
#'     }
#'     out <- x$xpos; names(out) <- x$ofID
#'     return(out)
#'   })
#'   names(xv) <- NULL
#'   xv <- do.call(c, xv)
#'   return(xv)
#' }

