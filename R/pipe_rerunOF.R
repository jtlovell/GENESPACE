pipe_rerunOF = function(dir.list,
                        genomeIDs,
                        gene.index,
                        verbose = T,
                        blk,
                        map,
                        n.cores = 1,
                        min.block.size = 3){
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
    if(what == "gene.num"){
      out <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
        pull_gff(gff = gff,
                 blk.line = blk[i,])$gene.num)
    }else{
      out <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
        pull_gff(gff = gff,
                 blk.line = blk[i,])$id)
    }
    names(out) <- blk$block.id
    return(out)
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

  check_blockSize = function(orig.blk, new.blk){
    ob = with(orig.blk,
              data.table(block.id = block.id,
                         n.mapping1 = n.mapping))
    nb = with(new.blk,
              data.table(block.id = block.id,
                         n.mapping2 = n.mapping))
    setkey(ob, block.id)
    setkey(ob, block.id)
    m = merge(ob, nb)
    m$n.lost = with(m, n.mapping1 - n.mapping2)
    m$n.lost[m$n.lost<0]<-0
    m$n.gained = with(m, n.mapping2 - n.mapping1)
    m$n.gained[m$n.gained<0]<-0
    m$prop.lost = with(m, n.lost/n.mapping1)
    m$prop.gained = with(m, n.gained/n.mapping1)

    return(m)
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
    blast.dir = dir.list$cull.blast,
    verbose = T)

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

  if(verbose)
    cat("Checking blocks\n\t")

  new.block.meta = check_blockSize(orig.blk = blk,
                                   new.blk = tmp$block)
  drop.these = new.block.meta$block.id[new.block.meta$n.mapping2< min.block.size]
  if(verbose)
    cat("Found", length(drop.these),
        "blocks that are now smaller than", min.block.size)

  nbm = new.block.meta
  nbm = nbm[nbm$block.id %in% drop.these,]
  map.out = tmp$map[!tmp$map$block.id %in% drop.these,]
  blk.out = tmp$block[!tmp$block$block.id %in% drop.these,]
  if(verbose)
    cat("\n\tRemaining blocks contain", nrow(map.out), "hits\n\tThat is",
        round(((nrow(map.out)-nrow(map))/nrow(map))*100),
        "% more genes in orthogroups\n\tDone!\n")
  ########################################################
  return(list(blast = all.blast,
              block = blk.out,
              map = map.out,
              id.list = id.list,
              all.block.metadata = new.block.meta))
}
