#' @title finalize_brkpts
#'
#' @description
#' \code{finalize_brkpts} finalize_brkpts
#'
#' @param clmap data.table, containing the merged gff and blast results,
#' must be built by finalize_blocks.
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param verbose logical, should updates be printed to the console?
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A list of length 2, block and map, as output by make_blocks.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
finalize_brkpts <- function(clmap,
                            gff,
                            verbose = T){
  find_ovlpBlks <- function(map, buffer = 0, calc.n.dup = T){
    map[,ovlrank1 := frankv(map,
                            cols = c("genome1","chr1","start1","end1"),
                            ties.method = "random")]
    map[,ovlrank2 := frankv(map,
                            cols = c("genome2","chr2","start2","end2"),
                            ties.method = "random")]
    map[,ovlrank1 := as.numeric(as.factor(paste(genome1, chr1)))*1000 + ovlrank1]
    map[,ovlrank2 := as.numeric(as.factor(paste(genome2, chr2)))*1000 + ovlrank2]
    b1 <- map[,list(start = min(ovlrank1) - buffer,
                    end = max(ovlrank1) + buffer),
              by = "block.id"]
    b2 <- map[,list(start = min(ovlrank2) - buffer,
                    end = max(ovlrank2) + buffer),
              by = "block.id"]
    setkey(b1, start, end)
    setkey(b2, start, end)
    ovl1 <- foverlaps(b1, b1, type = "any", which = TRUE)
    ovl2 <- foverlaps(b2, b2, type = "any", which = TRUE)
    ovl1 <- subset(ovl1, xid != yid)
    ovl2 <- subset(ovl2, xid != yid)

    ovl.out1 <- data.table(block.id = b1$block.id[ovl1$xid],
                           block.ovl = b1$block.id[ovl1$yid])
    ovl.out2 <- data.table(block.id = b2$block.id[ovl2$xid],
                           block.ovl = b2$block.id[ovl2$yid])
    if (calc.n.dup) {
      spl.map <- split(map, by = "block.id")
      spl.blk1 <- split(b1, by = "block.id")
      spl.blk2 <- split(b2, by = "block.id")
      ovl.out1[,n.dup1 := apply(ovl.out1, 1, function(x){
        sum(duplicated(c(unique(spl.map[[x[1]]]$id1), unique(spl.map[[x[2]]]$id1))))
      })]
      ovl.out1[,prop.ovl1 := apply(ovl.out1, 1, function(x){
        x1 <- spl.blk1[[x[1]]]
        x2 <- spl.blk1[[x[2]]]
        u1 <- x1$start:x1$end
        u2 <- x2$start:x2$end
        return(length(intersect(u1,u2))/length(unique(c(u1,u2))))
      })]
      ovl.out2[,n.dup2 := apply(ovl.out2, 1, function(x){
        sum(duplicated(c(unique(spl.map[[x[1]]]$id2), unique(spl.map[[x[2]]]$id2))))
      })]
      ovl.out2[,prop.ovl2 := apply(ovl.out2, 1, function(x){
        x1 <- spl.blk2[[x[1]]]
        x2 <- spl.blk2[[x[2]]]
        u1 <- x1$start:x1$end
        u2 <- x2$start:x2$end
        return(length(intersect(u1,u2))/length(unique(c(u1,u2))))
      })]
    }
    ovl.out <- merge(ovl.out1, ovl.out2)
    return(list(ovl1 = ovl.out1,
                ovl2 = ovl.out2,
                ovlboth = ovl.out))
  }


  solve_blkBreakpoints <- function(map, block.index, solve = "both"){
    map.in <- data.table(map)
    bo <- map.in[,list(start1.ord = min(order1),
                       end1.ord = max(order1),
                       start2.ord = min(order2),
                       end2.ord = max(order2)),
                 by = list(genome1, genome2, chr1, chr2, block.id)]
    if (nrow(block.index) > 0) {
      map <- subset(map, block.id %in% unique(unlist(block.index)))

      loo <- map[,list(start1 = min(order1),
                       end1 = max(order1),
                       start2 = min(order2),
                       end2 = max(order2)),
                 by = list(genome1, genome2, chr1, chr2, block.id)]

      lo <- rbindlist(apply(block.index, 1, function(i){

        x <- subset(map, block.id %in% i)
        splx <- split(x, by = "block.id")
        y <- do.call(rbind, splx[i])
        y1 <- y[!duplicated(id1), ]

        if (solve %in% c("1","both")) {
          setkey(y1, order1)
          for (j in 1:5) {
            y1[,rlid1 := rleid(block.id)]
            y1[,rllen1 := .N, by = rlid1]
            y1 <- subset(y1, rllen1 > j)
          }
        }

        setkey(y1, block.id)
        if (solve %in% c("2","both")) {
          y1 <- y1[!duplicated(id2), ]
          setkey(y1, order2)
          for (j in 1:5) {
            y1[,rlid2 := rleid(block.id)]
            y1[,rllen2 := .N, by = rlid2]
            y1 <- subset(y1, rllen2 > j)
          }
        }

        out <- y1[,list(start1 = min(order1),
                        end1 = max(order1),
                        start2 = min(order2),
                        end2 = max(order2)),
                  by = list(genome1, genome2, chr1, chr2, block.id)]
        return(out)
      }))
      loo <- lo[,list(start1.ord = max(start1),
                      end1.ord = min(end1),
                      start2.ord = max(start2),
                      end2.ord = min(end2)),
                by = list(genome1, genome2, chr1, chr2, block.id)]
      bo <- subset(bo, !block.id %in% loo$block.id)
      bo <- rbind(bo, loo)
    }
    return(bo)
  }


  if (verbose)
    cat("Finding overlapping blocks ... ")

  g <- data.table(gff)

  clmap[,order1 := g$order[match(id1, g$id)]]
  clmap[,order2 := g$order[match(id2, g$id)]]
  spl.map <- split(clmap, by = "blk.type")
  ovl.list <- lapply(spl.map, function(x){
    find_ovlpBlks(map = x)
  })
  if (verbose)
    cat("Done!\nBuilding optimal breakpoints ... ")

  brl <- lapply(1:length(ovl.list), function(i){
    m <- spl.map[[i]]
    o <- ovl.list[[i]]
    br <- solve_blkBreakpoints(map = m,
                               block.index = o$ovlboth)
    splm <- split(m, by = "block.id")
    splb <- split(br, by = "block.id")
    m <- rbindlist(lapply(names(splm), function(j){
      y <- splb[[j]]
      return(subset(splm[[j]], order1 >= y$start1.ord &
                      order1 <= y$end1.ord &
                      order2 >= y$start2.ord &
                      order2 <= y$end2.ord))
    }))
    o1 <- subset(o$ovl1, n.dup1 < 5 & prop.ovl1 < .1)
    br <- solve_blkBreakpoints(map = m,
                               block.index = o1,
                               solve = "1")
    splm <- split(m, by = "block.id")
    splb <- split(br, by = "block.id")
    m <- rbindlist(lapply(names(splm), function(j){
      y <- splb[[j]]
      return(subset(splm[[j]], order1 >= y$start1.ord &
                      order1 <= y$end1.ord &
                      order2 >= y$start2.ord &
                      order2 <= y$end2.ord))
    }))
    o2 <- subset(o$ovl2, n.dup2 < 5 & prop.ovl2 < .1)
    br <- solve_blkBreakpoints(map = m,
                               block.index = o2,
                               solve = "2")
    br[,block.type := i]
    return(br)
  })

  br <- rbindlist(brl)

  if (verbose)
    cat("Done!\nPropagating gff onto block breakpoints ... ")
  br[,block.id := paste0("blk_",as.numeric(as.factor(paste(genome1, genome2, chr1, chr2, block.id, block.type))))]

  g1 <- data.table(g)
  setnames(g1, paste0(colnames(g1),"1"))
  g2 <- data.table(g)
  setnames(g2, paste0(colnames(g2),"2"))

  gindex1 <- rbindlist(lapply(1:nrow(br), function(i)
    data.table(block.id = br$block.id[i],
               genome2 = br$genome2[i],
               order1 = br$start1.ord[i]:br$end1.ord[i])))
  gindex1 <- merge(gindex1, g1, by = "order1")
  gi1.spl <- split(gindex1, by = "block.id")

  gindex2 <- rbindlist(lapply(1:nrow(br), function(i)
    data.table(block.id = br$block.id[i],
               genome1 = br$genome1[i],
               order2 = br$start2.ord[i]:br$end2.ord[i])))
  gindex2 <- merge(gindex2, g2, by = "order2")
  gi2.spl <- split(gindex2, by = "block.id")

  bl.ids <- assblks[,c("id1","id2")]
  if (verbose)
    cat("Done!\nBuilding new hit database for", length(names(gi1.spl)),"blocks\n")
  idl <- lapply(names(gi1.spl), function(i){
    if (which(names(gi1.spl) == i) %% 100  == 0 | i == length(names(gi1.spl)))
      cat("\tCompleted", which(names(gi1.spl) == i),"/",length(names(gi1.spl)),"blocks\n")
    bli <- subset(bl.ids, id1 %in% gi1.spl[[i]]$id1 & id2 %in% gi2.spl[[i]]$id2)
    out <- merge(gi1.spl[[i]], bli, by = "id1", all = T)
    out <- merge(gi2.spl[[i]], out, by = c("genome1","genome2","id2","block.id"), all = T)
    return(out)
  })
  if (verbose)
    cat("\tDone!\n")
  ido <- rbindlist(idl)
  return(ido)
}
