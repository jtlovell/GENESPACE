  #' @title pipe_syntenicBlocks
  #'
  #' @description
  #' \code{pipe_syntenicBlocks} pipe_syntenicBlocks
  #'
  #' @param genomeIDs Character vector of genome IDs to consider
  #' @param clean.by Character, either 'merge' or 'dbscan'. Alternatively,
  #' if NULL (default), then merge is used if the min.block.size >=5 and
  #' dbscan is used if min.block.size< 5. While dbscan method may be a
  #' bit less sensitive, the massive number of blocks when min.block.size is
  #' small necessitates a non-iterative cleaning procedure.
  #' @param dir.list Directory list, produced by `check_environment`.
  #' @param peptide.dir Directory containing peptide fasta files
  #' @param gff.dir Directory containing peptide fasta files
  #' @param tmp.dir Directory containing peptide fasta files
  #' @param blast.dir Directory containing peptide fasta files
  #' @param block.dir Directory containing peptide fasta files
  #' @param cull.blast.dir Directory containing peptide fasta files
  #' @param mcscan.dir Directory containing peptide fasta files
  #' @param results.dir Directory containing peptide fasta files
  #' @param min.block.size Numeric, minimum required block size
  #' @param plotit Logical, should plots be generated?
  #' @param min.hit.density Numeric, the proportion of hits possible
  #' in a block for that block to be retained.
  #' @param min.unique.hits The number of unique hits required in each genome
  #' for a block to be retained.
  #' @param max.hit.density The maximum hit density in a block
  #' @param initial.mergeBuffer Numeric, the buffer for an initial scan
  #' to merge overlapping blocks. Should be set to either 0 or 1, unless
  #' one is ok with sacraficing speed for accuracy, in which case,
  #' it can be set anywhere up to min.block.size - 1.
  #' @param final.mergeBuffer Numeric, the buffer to search for blocks to merge
  #' in the last merge step. Should not exceed min.block.size -1.
  #' @param max.propSmallChr2Merge Numeric, the maximum size of two blocks to be
  #' merged is set based on this parameter. If set to 1, entire chromosomes
  #' may be merged, if set to 0, no merging ever happens. Values between 0 and
  #' 1, result in a maximum block size, equal to the proportion of the number of
  #' blast hits on the smallest chromosome.
  #' @param n.mergeIter Numeric, the number of iterations that the merge should
  #' go through.
  #' @param return.initial Logical, should the initial orthofinder results
  #' be returned?
  #' @param return.mcscanRaw Logical, should the unmerged results from MCScanX
  #' be returned?
  #' @param return.overlapMerged Logical, should the overlap-merged
  #' orthofinder results be returned?
  #' @param verbose Logical, should updates be printed?
  #' @param ... Not in use yet.
  #' @details More here
  #' @return Nothing.
  #'
  #' @examples
  #' \dontrun{
  #' none yet
  #' }
  #' @import data.table
  #' @export
  pipe_syntenicBlocks <- function(genomeIDs,
                                  dir.list = NULL,
                                  peptide.dir = NULL,
                                  gff.dir = NULL,
                                  tmp.dir = NULL,
                                  blast.dir = NULL,
                                  block.dir = NULL,
                                  cull.blast.dir = NULL,
                                  mcscan.dir = NULL,
                                  results.dir = NULL,
                                  verbose = TRUE,
                                  min.block.size = 3,
                                  min.unique.hits = min.block.size,
                                  plotit = FALSE,
                                  min.hit.ratio = .2,
                                  initial.mergeBuffer = 1,
                                  n.cores = 1,
                                  final.mergeBuffer = NULL,
                                  return.initial = TRUE,
                                  return.mcscanRaw = FALSE,
                                  return.overlapMerged = FALSE,
                                  max.size2merge = 100,
                                  clean.by = NULL,
                                  gap.multiplier = 8,
                                  mcscan.m.param = NULL,
                                  merge.overlaps = F,
                                  clean.byDBscan = T,
                                  cull.byDBscan = T,
                                  cull.byMCscan = T,
                                  return.ogblast = T,
                                  ...){


    #######################################################
    all.dirs <- !any(is.null(peptide.dir),
                     is.null(gff.dir),
                     is.null(tmp.dir),
                     is.null(blast.dir),
                     is.null(cull.blast.dir),
                     is.null(mcscan.dir),
                     is.null(block.dir))
    if (is.null(dir.list) & !all.dirs) {
      stop("Either a directory list or each directory parameter must be specified\n")
    }

    if (!is.null(dir.list)) {
      peptide.dir <- dir.list$peptide
      gff.dir <- dir.list$gff
      tmp.dir <- dir.list$tmp
      blast.dir <- dir.list$blast
      cull.blast.dir <- dir.list$cull.blast
      block.dir <- dir.list$block
      results.dir <- dir.list$results
      mcscan.dir <- dir.list$mcscan
    }
    #######################################################

    #######################################################
    if (verbose)
      cat("##########\n# - Part 1: Parsing orthofinder ...\n#\tProgress:\n")

    if(!is.null(mcscan.m.param)){
      mcsp <- paste("-a -s", min.block.size,
                    "-m", mcscan.m.param,
                    "-w 2 -e 1")
    }else{
      mcsp <- paste("-a -s", min.block.size,
                    "-m", min.block.size*gap.multiplier,
                    "-w 2 -e 1")
    }

    init.results <- process_orthofinder(
      gff.dir = gff.dir,
      genomeIDs = genomeIDs,
      blast.dir = blast.dir,
      mcscan.dir = mcscan.dir,
      eps.radius = c(min.block.size*5*gap.multiplier,
                     min.block.size*2*gap.multiplier,
                     min.block.size*1*gap.multiplier),
      n.mappingWithinRadius = c(min.block.size,
                                min.block.size,
                                min.block.size),
      cull.byDBscan = cull.byDBscan,
      mcscan.param = mcsp,
      return.ogblast = return.ogblast)

    gff <- init.results$gff

    #######################################################

    #######################################################
    if (verbose)
      cat("##########\n# - Part 2: Building collinear blocks with MCScanX...\n")


    if(!is.null(mcscan.m.param)){
      mcsp <- paste("-a -s", min.block.size,
                    "-m", mcscan.m.param,
                    "-w 2")
    }else{
      mcsp <- paste("-a -s", min.block.size,
                    "-m", min.block.size*gap.multiplier,
                    "-w 2")
    }

    synteny.results <- pipe_mcs(blast = init.results$blast$map,
                                gff = gff,
                                mcscan.dir = mcscan.dir,
                                mcscan.param = mcsp)

    blk <- synteny.results$block
    map <- synteny.results$map

    synteny.results <- make_blocks(map)

    #######################################################
    #######################################################
    #######################################################

    if(!clean.byDBscan){
      merged_blk <- synteny.results
      merged_blk.out <- NULL
    }else{
      if (verbose)
        cat("##########\n# - Part 3: Cleaning / merging via dbscan\n",
            "initial n blocks / mappings =",
            nrow(synteny.results$block),
            "/",
            nrow(synteny.results$map),"\n\t")
      tmp = synteny.results$map
      tmp$unique = with(tmp, paste(genome1, genome2, chr1, chr2))
      spl.map = split.data.table(tmp, "unique")

      min.rad = min.block.size
      merged_map = rbindlist(lapply(spl.map, function(tmp){
        x <- run_dbs(y = tmp[,c("rank1","rank2"),with = F],
                     eps.radius = min.rad+.1,
                     mappings = min.block.size)
        tmp$block.id = x$cluster
        cols = sample(rainbow(length(x$cluster)),size = length(x$cluster), replace = F)
        return(tmp)
      }))
      merged_map = merged_map[merged_map$block.id != 0,]

      merged_map$block.id = with(merged_map, as.numeric(as.factor(paste(unique, block.id))))
      merged_blk = make_blocks(merged_map,rename.blocks = F, rerank = T)
      merged_blk.out = merged_blk
      if (verbose)
        cat("cleaned n blocks / mappings =",
            nrow(merged_blk$block),
            "/",
            nrow(merged_blk$map),"\n")
    }

    #######################################################
    #######################################################
    #######################################################

    if(!merge.overlaps){
      merged.overlaps <- NULL
    }else{
      if (verbose)
        cat("##########\n# - Part 4: Merging overlapping blocks\n",
            "initial n blocks / mappings =",
            nrow(merged_blk$block),
            "/",
            nrow(merged_blk$map),"\n\t")

      merged.overlaps <- merge_blocks(
        blk = merged_blk$block,
        map = merged_blk$map,
        buffer = 0,
        n.cores =  n.cores,
        max.size2merge = 1e6)

      if (verbose)
        cat("merged to n blocks / mappings =",
            nrow(merged.overlaps$block),
            "/",
            nrow(merged.overlaps$map),"\n")

    }

    out <- list(synteny.results = synteny.results,
                init.results = init.results,
                merged.dbs = merged_blk.out,
                merged.results = merged.overlaps)
    return(out)

  }
