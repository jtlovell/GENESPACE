#' @title build pwBlocks
#'
#' @description
#' \code{build_pwBlocks} build_pwBlocks
#'
#' @param pw.of pw.of
#' @param comb comb
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
#' @export
merge_ofGff <- function(comb,
                        pw.of,
                        verbose = T){
  pw.of2 <- lapply(1:length(comb), function(i){
    x = pw.of[[i]]
    if(verbose)
      cat(paste0("\t",comb[[i]][1]),"<-->", comb[[i]][2],"... ")
    x[,n := length(unique(genome)), by = og]
    x <- subset(x, n > 1)

    x$og.id <- paste0(comb[[i]][1],"_",comb[[i]][2],"_", as.numeric(as.factor(x$og)))

    x$n <- NULL; x$order <- NULL; x$strand <- NULL; x$og <- NULL

    x1 <- data.table(x)
    x2 <- data.table(x)
    setnames(x1, paste0(colnames(x1), "1"))
    setnames(x2, paste0(colnames(x2), "2"))
    setkey(x1, id1)
    setkey(x2, id2)
    g1.ids <- x$id[x$genome == comb[[i]][1]]
    g2.ids <- x$id[x$genome == comb[[i]][2]]

    y <- read_allBlasts(gff = subset(gff, genome %in% comb[[i]]),
                        genomeIDs = comb[[i]],
                        blast.dir = dirs$blast)
    y2 <- data.table(y[, c(2, 1, 3:6, 9:10, 7:8, 11:12)])
    setnames(y2, colnames(y))
    y0 <- rbind(y, y2)
    y <- subset(y, (id1 %in% g1.ids &
                      id2 %in% g2.ids) |
                  (id1 %in% g1.ids & id2 %in% g1.ids) |
                  (id1 %in% g2.ids & id2 %in% g2.ids))

    y$neg.score <- y$score * (-1)
    setkey(y, neg.score)
    y <- y[!duplicated(y[,c("id1","id2")]),]
    y$neg.score <- NULL

    if(verbose)
      cat("initial hits =", nrow(y),"... ")
    setkey(y, id2)
    out <- merge(x2, y)
    setkey(out, id1)
    out <- merge(x1, out)
    out <- subset(out, og.id1 == og.id2)

    if(verbose)
      cat("found", nrow(out),"in orthogroups\n ")
    return(out)
  })
  map <- rbindlist(pw.of2)
  map$og.id2 <- NULL
  setnames(map, "og.id1","og.id")
  return(map)
}
