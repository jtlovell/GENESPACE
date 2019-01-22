split_fastaByBlock <- function(blk,
                               assembly.dir,
                               write.dir,
                               tmp.dir,
                               run.getfasta = T){
  ########################################################
  ########################################################
  get_fasta <- function(bed,
                        bed.file,
                        ass.file,
                        fa.file){

    write.table(bed,
                file = bed.file,
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = F)

    system(paste("bedtools getfasta",
                 "-fi", ass.file,
                 "-bed", bed.file,
                 "-name+ -fo", fa.file))
  }
  ########################################################
  ########################################################
  convert_blk2bed <- function(blk){
    bed.dt1 <- data.table(blk[,c("chr1","start1","end1","genome1")])
    bed.dt2 <- data.table(blk[,c("chr2","start2","end2","genome2")])
    bed.dt1$start1 <- bed.dt1$start1 - 1
    bed.dt1$end1 <- bed.dt1$end1 - 1
    bed.dt2$start2 <- bed.dt2$start2 - 1
    bed.dt2$end2 <- bed.dt2$end2 - 1
    bed.dt1$name <- paste0(blk$genome1, "xxxx", blk$block.id)
    bed.dt2$name <- paste0(blk$genome2, "xxxx", blk$block.id)
    setnames(bed.dt1,c("chr","start","end","genome","name"))
    setnames(bed.dt2,c("chr","start","end","genome","name"))
    bed.dt <- rbind(bed.dt1, bed.dt2)
    return(bed.dt)
  }
  ########################################################
  ########################################################
  if(verbose)
    cat("Preparing the block data.table ...")
  bed.dt<-convert_blk2bed(blk)

  nb <- ifelse(nrow(bed.dt) < 1000, 100,
               ifelse(nrow(bed.dt) < 5000, 500,
                      ifelse(nrow(bed.dt) < 10000, 1000, 5000)))

  bed.dt$unique <- 1:nrow(bed.dt)
  spl = split(bed.dt, "unique")
  if(verbose)
    cat("Done!\n")
  if(verbose)
    cat("Running bedtools getfasta ... Completed:\n\t")
  fa.locs <- unlist(lapply(1:length(spl), function(i){
    x <- spl[[i]]
    if (verbose)
      if (i %% nb == 0)
        cat(i, "/",
            nrow(bed.dt), "\n\t")

    bed <- x[,c("chr","start","end","genome")]
    bed.id <- x$name

    ass.file <- file.path(assembly.dir,
                          paste0(bed$genome, ".fa"))
    fa.file <- file.path(write.dir,
                         paste0(bed.id, ".fa"))
    bed.file <- file.path(tmp.dir,
                          paste0(bed.id, ".bed"))

    if(run.getfasta){
      get_fasta(bed = bed,
                ass.file = ass.file,
                fa.file = fa.file,
                bed.file = bed.file)

    }
    return(fa.file)
  }))
  bed.dt$unique<-NULL
  bed.dt$fa.file <- fa.locs
  if(verbose)
    cat("Done!\n")
  return(bed.dt)
}
