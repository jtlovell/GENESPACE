#
# #######################################################
# #######################################################
# # Subset blast results within blocks
#
# subset_blast <- function(blast.file,
#                         genenum.list,
#                         n.cores,
#                         verbose){
#   if(verbose)
#     cat("Subsetting:",blast.file)
#   f <- fread(blast.file,
#              header = F,
#              stringsAsFactors = F,
#              check.names = F)
#
#   if(verbose)
#     cat(paste0(" (initial hits = ",nrow(f),") ")
#
#   setkey(f, V1, V2)
#
#   fl <- rbindlist(mclapply(genenum.list, mc.cores = n.cores, function(x){
#     return(f[f$V1  %in% x & f$V2 %in% x,])
#   }))
#
#   if(verbose)
#     cat("to",nrow(f),"hits in blocks\n")
#
#   write.table(fl,
#               sep = "\t",
#               row.names = F,
#               col.names = F,
#               quote = F,
#               file = blast.file)
# }
#
#
# remake_blast <- function(blast.dir,
#                          cull.blast.dir,
#                          genenum.list,
#                          n.cores,
#                          verbose){
#   if(verbose)
#     cat("Copying blast results to",cull.blast.dir,"... ")
#   if (dir.exists(cull.blast.dir)) {
#     nu <- file.remove(list.files(cull.blast.dir,
#                                  full.names = T))
#   }
#   dir.create("mkdir", cull.blast.dir)
#
#   for(i in list.files(blast.dir,
#                       full.names = T)){
#     nu <- file.copy(i,cull.blast.dir)
#   }
#
#   if(verbose)
#     cat("Done!\n")
#
#   nu <- file.remove(file.path(cull.blast.dir,
#                               "Orthogroups.txt"))
#   blast.files <- list.files(cull.blast.dir,
#                            full.names = T,
#                            pattern = "Blast")
#
#   for(i in blast.files){
#     subset_blast(blast.file = i,
#                  genenum.list = genenum.list,
#                  n.cores = n.cores,
#                  verbose = verbose)
#   }
# }
#
#
# #######################################################
# #######################################################
# # merge gff with gene numbers
# pull_gff <- function(gff,
#                     blk.line,
#                     gene.index){
#   x <- blk.line
#   setkey(gene.index, id)
#
#   g1 <- gff[with(gff, genome == x$genome1 &
#                   chr == x$chr1 &
#                   start <= x$end1 &
#                   end >= x$start1),]
#   g1$block.id <- x$block.id
#   setkey(g1, id)
#
#   g2 <- gff[with(gff, genome == x$genome2 &
#                   chr == x$chr2 &
#                   start <= x$end2 &
#                   end >= x$start2),]
#   g2$block.id <- x$block.id
#   setkey(g2, id)
#
#   g1 <- merge(gene.index, g1)
#   g2 <- merge(gene.index, g2)
#
#   return(list(g1, g2))
# }
#
# #######################################################
# #######################################################
# # read blast files
# read_allBlast <- function(blast.dir){
#   blast.files <- list.files(blast.dir,
#                            full.names = T,
#                            pattern = "Blast")
#   out <- rbindlist(lapply(blast.files, function(x)
#     fread(x, header = F,
#           stringsAsFactors = F,
#           check.names = F,
#           col.names = c("gene.num1", "gene.num2",
#                         "perc.iden", "align.length",
#                         "n.mismatch", "n.gapOpen",
#                         "q.start", "q.end",
#                         "s.start", "s.end",
#                         "eval", "score"))))
#   return(out)
# }
#
# #######################################################
# #######################################################
# # buid gff networks from orthogroups
#
# make_mapFromOGs <- function(gff.wNum,
#                             cull.blast.dir){
#
#   og <- readLines(file.path(cull.blast.dir,
#                             "Orthogroups.txt"))
#
#   og <- lapply(og, function(x)
#     strsplit(x, " ")[[1]])
#   ons <- sapply(og, function(x)
#     x[1])
#   names(og) <- ons
#
#   og <- lapply(og, function(x)
#     x[-1])
#   og.name <- names(og)
#   og.length <- sapply(og, length)
#
#   od <- data.table(og.id = rep(og.name, og.length),
#                   id = unlist(og),
#                   stringsAsFactors = F)
#   setkey(od, id)
#
#   gffn <- rbindlist(unlist(gff.wNum,
#                            recursive = F))
#   setkey(gffn,id)
#   gffn <- merge(od, gffn)
#
#   gffn[, complete := all(duplicated(block.id) |
#                            duplicated(block.id, fromLast = T)),
#        by = og.id]
#   setkey(gffn, id)
#   return(gffn)
# }
#
#
#
# #######################################################
# #######################################################
# # combine block information with gff of orthogroups
# make_blockMetadata <- function(cds.md,
#                                gffog,
#                                blk){
#   gffog <- merge(cds.md, gffog)
#
#   gffog.incomplete <- gffog[!gffog$complete,]
#   gffog.incomplete$unique.block <- with(gffog.incomplete,
#                                         paste0(genome,"_", block.id))
#
#   gene.block <- gffog[,list(block = unique(block.id)),
#                       by = id]
#   gene.block <- split(gene.block$block, gene.block$id)
#
#   blk.list <- split.data.table(blk, "block.id")
#   spl.gff <- split.data.table(gffog.incomplete, "og.id")
#
#   blk.md <- blk[,list(genome = c(genome1, genome2),
#                       chr = c(chr1,chr2),
#                       start = c(start1,start2),
#                       end = c(end1,end2)),
#                 by = list(block.id)]
#
#   blk.md$unique.block <- with(blk.md, paste0(genome,"_", block.id))
#   setkey(blk.md, unique.block, block.id, genome)
#
#   blk.md2 <- blk.md[,c("block.id","genome","unique.block")]
#   setkey(blk.md2, unique.block)
#
#   return(list(spl.incompleteGff = spl.gff,
#               blk.metadata = blk.md,
#               simple.blk.metadata = blk.md2))
#
# }
#
# #######################################################
# #######################################################
# # find orphan genes in blocks
#
# find_orphans <- function(spl.gff,
#                          blk.md2,
#                          n.cores){
#   y <- rbindlist(mclapply(spl.gff, mc.cores = n.cores, function(x){
#     xmd <- blk.md2[blk.md2$block.id %in% unique(x$block.id), ]
#     xmd <- xmd[!xmd$unique.block %in% x$unique.block, ]
#     xmd$gene2map <- x$id[which.max(x$length)]
#     xmd$og.id <- x$og.id[1]
#     return(xmd)
#   }))
#
#   setkey(y, unique.block, block.id, genome)
#
#   return(y)
# }
#
# #######################################################
# #######################################################
# # load annotation to memory
# load.annotations = function(genomeIDs,
#                             cds.dir,
#                             peptide.dir,
#                             assembly.dir){
#
#   cds.fastas <- do.call(c, lapply(genomeIDs, function(i)
#     readDNAStringSet(file.path(cds.dir,
#                                paste0(i, ".fa")))))
#
#   pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
#     readAAStringSet(file.path(peptide.dir,
#                               paste0(i, ".fa")))))
#
#   fais <- rbindlist(lapply(genomeIDs, function(i){
#     tmp <- read.delim(file.path(assembly.dir,
#                                 paste0(i, ".fa.fai")),
#                       header = F,
#                       stringsAsFactors = F,
#                       col.names = c("chr",  "chr.length",
#                                     "v1", "v2", "v3"))[, 1:2]
#     tmp$genome <- i
#     return(tmp)
#   }))
#
#   setkey(fais, genome, chr)
#   return(list(peptide = pep.fastas,
#               cds = cds.fastas,
#               chrlen = fais))
# }
#
