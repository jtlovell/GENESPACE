#' @title Fast split of data.table
#'  @description
#' \code{split.data.table} Much faster than base split.
#' @param x a data.table
#' @param f a vector, not needed if by is specified
#' @param drop drop
#' @param flatten non-recursive unlisting
#' @param ... additional arguments passed to data.table
#' @details Nothing yet
#' @return list of split data.tables
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
split.data.table <- function(x,
                             f,
                             drop = FALSE,
                             by,
                             flatten = FALSE,
                             ...){
  if (missing(by) &&
     !missing(f))
    by = f

  stopifnot(!missing(by),
            is.character(by),
            is.logical(drop),
            is.logical(flatten),
            !".ll" %in% names(x),
            by %in% names(x),
            !"nm" %in% by)

  if (!flatten) {
    .by = by[1L]
    tmp = x[, list(.ll = list(.SD)),
            by = .by,
            .SDcols = if (drop) setdiff(names(x), .by) else names(x)]
    setattr(ll <- tmp$.ll,
            "names",
            tmp[[.by]])
    if (length(by) > 1L) {
      return(lapply(ll,
                    split.data.table,
                    drop = drop,
                    by = by[-1L]))
    } else {
      return(ll)
    }
  } else {
    tmp = x[, list(.ll = list(.SD)),
            by = by,
            .SDcols = if (drop) setdiff(names(x), by) else names(x)]
    setattr(ll <- tmp$.ll,
            'names',
            tmp[, .(nm = paste(.SD, collapse = ".")),
                by = by,
                .SDcols = by]$nm)
    return(ll)
  }
}
