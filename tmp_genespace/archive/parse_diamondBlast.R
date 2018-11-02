#' @title Process mapping files from reciprocal blasts
#'
#' @description
#' \code{parse_diamondBlast} Combine two mapping files from reciprocal blasts
#' and subset to high-confidence regions with dense hits.
#'
#' @param blast1 The file path of the blast result for genome 1 is the querry
#' @param gff1 The file path of the gff3-formatted annotation for genome 1
#' @param id1 Single word identifier for genome 1. Used for file and column names.
#' @param ploidy1 Ploidy of genome 1. This informs the number of blast hits to
#' retain. For example, if the query is 4x, we will retain `nmapsPerHaplotype`
#' mappings for each haploid genome.
#' @param abbrev1 Abrreviation of genome 1. Appended to seqnames for MCScanX
#' input files.
#' @param blast2 The file path of the blast result for genome 2 is the querry
#' @param gff2 The file path of the gff3-formatted annotation for genome 2
#' @param id2 Single word identifier for genome 2. Used for file and column names.
#' @param ploidy2 Ploidy of genome 2. This informs the number of blast hits to
#' retain. For example, if the query is 4x, we will retain `nmapsPerHaplotype`
#' mappings for each haploid genome.
#' @param abbrev2 Abrreviation of genome 2. Appended to seqnames for MCScanX
#' input files.
#' @param nmapsPerHaplotype see ploidy
#' @param dbs_radii numeric vector of length equal to `dbs_mappingsInRadius`.
#' @param dbs_mappingsInRadius numeric vector of length equal to `dbs_radius`.
#' @param mcscan.dir Path to output directory for mcscan input files
#' @param topPerc The percent threshold for the maximum score to retain a mapping
#' @param plotit Logical, should a plot be drawn?
#' @param verbose Logical, should updates be printed?
#' @param ... Not currently in use
#' @details See `pipe_Diamond2MCScanX` documentation for details
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @import dbscan
#' @export
parse_diamondBlast = function(blast1, blast2,
                              gff1, gff2,
                              id2, id1,
                              ploidy1 = 2, ploidy2 = 2,
                              abbrev1 = NULL, abbrev2 = NULL,
                              nmapsPerHaplotype = 1,
                              plotit = T,
                              dbs_radii = c(100,50,20,40),
                              dbs_mappingsInRadius = c(10,8,5,10),
                              mcscan.dir,
                              topPerc = 50,
                              verbose = T){

  topProp = topPerc/100

  if(verbose)
    cat("\n4. Reading in mapping files ... \n\t",id1)
  d1 = suppressWarnings(
    data.table::fread(blast1,showProgress = F, verbose = F)
    )
  setnames(d1, c("v1", "v2","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))

  d1 <- d1[order(d1$v2, d1$v1, -d1$bitscore), ]
  d1 <- d1[!duplicated(d1[,1:2]),]

  d1[, rank := frank(bitscore, ties.method = "dense"),
     by = list(v1)]
  d1 <- d1[d1$rank <= ploidy1*nmapsPerHaplotype,]

  d1[, prop := bitscore/max(bitscore),
      by = list(v1)]
  d1 <- d1[d1$prop >= topProp,]

  if(verbose)
    cat(" found", nrow(d1),"with score within the top", ploidy1*nmapsPerHaplotype,"\n\t",id2)


  d2 = suppressWarnings(
    data.table::fread(blast2,showProgress = F, verbose = F)
  )
  d2 = d2[,c(2,1,3:6,9:10,7:8,11:12)]
  setnames(d2, c("v1", "v2","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))

  d2 <- d2[order(d2$v2, d2$v1, -d2$bitscore), ]
  d2 <- d2[!duplicated(d2[,1:2]),]

  d2[, rank := frank(bitscore, ties.method = "dense"),
     by = list(v1)]
  d2 <- d2[d2$rank <= ploidy2*nmapsPerHaplotype,]

  d2[, prop := bitscore/max(bitscore),
     by = list(v2)]

  d2 <- d2[d2$prop >= topProp,]

  if(verbose)
    cat(" found", nrow(d2),"with score within the top", ploidy2*nmapsPerHaplotype,"\n")

  if(verbose)
    cat("5. Culling by score ...")
  dat = rbind(d1,d2)
  setkey(dat, "bitscore")
  dat <- dat[order(dat$v2, dat$v1, -dat$bitscore), ]
  dat <- dat[!duplicated(dat[,1:2]),]

  if(verbose)
    cat("found", nrow(dat), "unique reciporal mappings\n")



  setnames(dat,1:2, c(id2,id1))

  g1 = parse_quickGff(gff1)
  g2 = parse_quickGff(gff2)
  setnames(g1,c(id1, paste0(colnames(g1)[-1],"_",id1)))
  setnames(g2,c(id2,  paste0(colnames(g2)[-1],"_",id2)))

  if(verbose)
    cat("6. Culling by density ...")
  m = merge(g1,
            merge(g2,dat,
                  by = id2, all.y = T),
            by = id1, all.y = T)

  m2 = loop_dbs(m,id1 = id1, id2= id2,
                radii = dbs_radii,
                mappings = dbs_mappingsInRadius,
                plotit = plotit)

  if(verbose)
    cat("\tretained", nrow(m2), "mappings in total\n")

  outdir = file.path(mcscan.dir)
  out.gffish = file.path(outdir, paste0(id1,"_",id2,".gff"))
  out.mapish = file.path(outdir, paste0(id1,"_",id2,".blast"))
  if(verbose)
    cat("Writing MCScanX formated output to:\n\t",outdir)

  allmap = data.frame(m2, stringsAsFactors = F)

  if(plotit){
    c1 = table(allmap[,2])
    c2 = table(allmap[,7])
    tp = allmap[allmap[,2] %in% names(c1)[c1>100] &
                  allmap[,7] %in% names(c2)[c2>100],]
    tp$x_a = frank(tp[,paste0(c("chr_","start_"),id1)],
                  ties.method = "dense")
    tp$x_b = frank(tp[,paste0(c("chr_","start_"),id2)],
                  ties.method = "dense")
    lab1 = tapply(tp$x_a,tp[,paste0("chr_",id1)], mean)
    end1 = tapply(tp$x_a,tp[,paste0("chr_",id1)], max)[-length(lab1)]

    lab2 = tapply(tp$x_b,tp[,paste0("chr_",id2)], mean)
    end2 = tapply(tp$x_b,tp[,paste0("chr_",id2)], max)[-length(lab2)]

    plot(tp$x_a, tp$x_b, axes = F, type= "n",
         xlab = paste0(id1, " genome order"),
         ylab = paste0(id2, " genome order"), asp = 1)
    axis(1, at = lab1, labels = names(lab1))
    axis(2, at = lab2, labels = names(lab2))
    segments(x0 = end1, x1 = end1, y0 = min(tp$x_b), y1 = max(tp$x_b),
             col = "lightgrey", lty = 2)
    segments(y0 = end2, y1 = end2, x0 = min(tp$x_a), x1 = max(tp$x_a),
             col = "lightgrey", lty = 2)
    pal  = rep(c("firebrick3","darksalmon",
                 "darkorange4","gold",
                 "dodgerblue4","skyblue",
                 "purple","violet")
               ,100)
    count = 0
    for(i in names(lab1)){
      count = count+1
      with(tp[tp[,2] == i,], points(x_a, x_b,
                                    col = add.alpha(pal[count],.25),
                                    pch = 16, cex = .25))
      title("Culled Diamond Protein blast mappings")
    }
  }

  g1 = allmap[,c(paste0("chr_",id1),id1, paste0(c("start_","end_"),id1))]
  g2 = allmap[,c(paste0("chr_",id2),id2, paste0(c("start_","end_"),id2))]

  g1[,1]<-paste0(abbrev1,as.numeric(as.factor(g1[,1])))
  g2[,1]<-paste0(abbrev2,as.numeric(as.factor(g2[,1])))

  colnames(g1) = colnames(g2)
  go = rbind(g1, g2)

  mo = allmap[,c(id1,id2,"pident","length","mismatch","gapopen",
                 "qstart","qend","sstart","send","evalue", "bitscore")]
  write.table(go, file = out.gffish, sep = "\t",
              row.names = F, col.names = F, quote = F)
  write.table(mo, file = out.mapish, sep = "\t",
              row.names = F, col.names = F, quote = F)


  return(m2)
}
