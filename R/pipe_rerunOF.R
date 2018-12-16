pipe_rerunOF = function(dir.list,
                        genomeIDs,
                        gene.index,
                        verbose = T,
                        blk,
                        n.cores = 1){
  ########################################################
  ########################################################
  convert_blast2blk = function(id.list,
                               blast,
                               n.cores = 1){

    fl <- rbindlist(mclapply(names(id.list), mc.cores = n.cores, function(i){
      x = id.list[[i]]
      tmp = blast[blast$id1  %in% x & blast$id2 %in% x,]
      tmp$block.id = i
      return(tmp)
    }))
    return(make_blocks(fl))
  }
  ########################################################
  ########################################################
  subset_blast <- function(blast.file,
                           genenum.list,
                           n.cores,
                           verbose){
    if (verbose)
      cat("\tSubsetting:",basename(blast.file))
    f <- fread(blast.file,
               header = F,
               stringsAsFactors = F,
               check.names = F)

    if (verbose)
      cat(paste0(" (initial hits = ",nrow(f),") "))


    f$ns = f$V12 * (-1)
    setkey(f, ns)
    fo <- f[!duplicated(f[,f[,1:2,with = F]]),]
    fo$ns <- NULL

    setkey(fo, V1, V2)

    fl <- rbindlist(mclapply(genenum.list, mc.cores = n.cores, function(x){
      return(fo[fo$V1  %in% x & fo$V2 %in% x,])
    }))


    if (verbose)
      cat("to", nrow(fl),
          "hits in blocks\n")

    write.table(fl,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = blast.file)
  }
  ########################################################
  ########################################################
  add_ofnum2gff = function(gff, gene.index){
    setkey(gff, id)
    setkey(gene.index, id)
    m = merge(gene.index, gff)
    return(m)
  }
  ########################################################
  ########################################################
  pull_idByBlk = function(gff, blk, what = "gene.num", n.cores = 1){
    gff.list <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
      pull_gff(gff = gff,
               blk.line = blk[i,]))
    gene.ids = lapply(gff.list, function(x) unique(x$id))
    gene.nums = lapply(gff.list, function(x) unique(x$gene.num))
    names(gene.ids) <- blk$block.id
    names(gene.nums) <- blk$block.id
    if(what == "gene.num"){
      return(gene.nums)
    }else{
      return(gene.ids)
    }
  }
  ########################################################
  ########################################################
  pull_gff <- function(gff,
                       blk.line){
    x <- blk.line

    g1 <- gff[with(gff, genome == x$genome1 &
                     chr == x$chr1 &
                     start <= x$end1 &
                     end >= x$start1),]

    g2 <- gff[with(gff, genome == x$genome2 &
                     chr == x$chr2 &
                     start <= x$end2 &
                     end >= x$start2),]

    return(rbind(g1,g2))
  }
  ########################################################
  ########################################################
  parse_gffByBlock = function(gff, blk){
    gff <- init.results$gff
    gff.wNum <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
      pull_gff(gff = gff,
               blk.line = blk[i,],
               gene.index = init.results$ortho.info$gene.index))

    gffn <- gff.wNum
    genenum.list <- lapply(gff.wNum, function(x)
      unique(c(x[[1]]$gene.num, x[[2]]$gene.num)))
    all.blkgff <- rbindlist(unlist(gff.wNum, recursive = F))
  }
  ########################################################
  ########################################################
  remake_blast <- function(blast.dir,
                           cull.blast.dir,
                           genenum.list,
                           n.cores,
                           verbose = T){
    if (verbose)
      cat("Copying blast results to", cull.blast.dir,"... ")
    if (dir.exists(cull.blast.dir))
      unlink(cull.blast.dir, recursive = T)

    if (!dir.exists(cull.blast.dir))
      dir.create(cull.blast.dir)

    files = list.files(blast.dir, full.names = T)

    nu <- file.copy(files,
                    cull.blast.dir)

    if (verbose)
      cat("Done!\n")

    nu <- file.remove(file.path(cull.blast.dir,
                                "Orthogroups.txt"))
    blast.files <- list.files(cull.blast.dir,
                              full.names = T,
                              pattern = "Blast")
    blast.files<-blast.files[grep(".txt$",blast.files)]

    if(verbose)
      cat("Subsetting blast files to hits within blocks\n")
    for (i in blast.files){
      subset_blast(blast.file = i,
                   genenum.list = genenum.list,
                   n.cores = n.cores,
                   verbose = verbose)
    }
  }
  ########################################################
  ########################################################

  ########################################################
  ########################################################

  ########################################################
  ########################################################

  ########################################################
  ########################################################

  ########################################################
  ########################################################

  ########################################################
  ########################################################

  gff <- import_gff(
    gff.dir = dir.list$gff,
    genomeIDs = genomeIDs)

  ########################################################
  if(verbose)
    cat("Adding orthofinder ids to gffs ...")

  gff = add_ofnum2gff(gff,
                      gene.index = gene.index)

  if(verbose)
    cat(" Done!\n")

  ########################################################
  if(verbose)
    cat("Generating list of gene numbers within each block ...")

  gene.list = pull_idByBlk(gff,
                           blk = blk,
                           n.cores = n.cores)

  if(verbose)
    cat(" Done!\n")

  ########################################################
  if(verbose)
    cat("Generating list of gene ids within each block ...")

  id.list = pull_idByBlk(gff,
                         blk = blk,
                         n.cores = n.cores, what = "id")

  if(verbose)
    cat(" Done!\n")

  ########################################################
  if(verbose)
    cat("Cleaning out the tmp and cull.blast directories ...")

  unlink(dir.list$tmp, recursive = T)
  dir.create(dir.list$tmp)

  unlink(dir.list$cull.blast, recursive = T)
  dir.create(dir.list$cull.blast)

  if(verbose)
    cat(" Done!\n")

  ########################################################
  remake_blast(blast.dir = dir.list$blast,
               cull.blast.dir = dir.list$tmp,
               genenum.list = gene.list,
               n.cores = n.cores,verbose = T)

  if(verbose)
    cat("\tDone!\n")

  ########################################################
  if(verbose)
    cat("Running orthofinder on culled blast data ...")

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = dir.list$tmp,
    blast.dir = dir.list$cull.blast,
    og.threads = n.cores,
    og.silent = F,
    verbose = T)

  if(verbose)
    cat("\tDone!\n")

  ########################################################

  gff$gene.num <- NULL
  of.blast <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = dir.list$tmp,
    verbose = T)

  print(of.blast$species.mappings)
  ########################################################

  all.blast <- import_blast(
    species.mappings = of.blast$species.mappings,
    genomeIDs = genomeIDs,
    orthogroups = of.blast$orthogroups,
    gff = gff,
    gene.index = of.blast$gene.index,
    verbose = T)

  ########################################################
  if(verbose)
    cat("Converting blast hits to block / map objects ...")

  tmp = convert_blast2blk(id.list = id.list,
                          blast = all.blast)

  if(verbose)
    cat(" Done!\n")

  ########################################################
  return(list(blast = all.blast,
              block = tmp$block,
              map = tmp$map,
              id.list = id.list))
}
