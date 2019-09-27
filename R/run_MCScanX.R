#' @title run_MCScanX
#'
#' @description
#' \code{run_MCScanX} run_MCScanX
#'
#' @param blast data.table, containing the merged gff and blast results.
#' Unlike the 'map' object, which really just needs merged gff coordinates,
#' this must have all the blast8 columns. See details.
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#' @param MCScanX.s.param numeric of length 1, that specifies the 's'
#' (block size) parameter for MCScanX.
#' @param MCScanX.m.param numeric of length 1, that specifies the 'm'
#' (n. gaps) parameter for MCScanX.
#' @param MCScanX.path file.path, specifying the location of the
#' MCScanX program. This directory must contain the executable
#' '/MCScanX'.
#' @param silent.mcs logical, should MCScanX progress be reported?
#'
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
run_MCScanX <- function(blast,
                        genomeIDs,
                        mcscan.dir,
                        gff,
                        MCScanX.s.param,
                        MCScanX.m.param,
                        MCScanX.path,
                        overwrite.output.dir = F,
                        verbose = T){

  parse_mcs <- function(mcs.file){
    if (length(readLines(paste0(mcs.file, ".collinearity"))) < 12) {
      return(NULL)
    }else{
      mcscan.raw <- read.delim(paste0(mcs.file, ".collinearity"),
                               sep = "\t", header = F,
                               comment.char = "#", strip.white = T,
                               stringsAsFactors = F)

      fac <- as.numeric(as.factor(sapply(as.character(mcscan.raw$V1), function(x)
        strsplit(x, "-")[[1]][1])))
      genes <- with(mcscan.raw, c(V2, V3))
      out <- data.table(id1 = mcscan.raw$V2,
                        id2 = mcscan.raw$V3,
                        block.id = fac)
      out <- data.table(out)
      setkey(out, id1, id2)
      return(out)
    }
  }
  # -- Specify the parameters
  if(verbose)
    cat("Prepping data for MCScanX ... ")
  MCScanX.tool <- file.path(MCScanX.path, "MCScanX")
  mcscan.param <- paste("-a -b 2 -s", MCScanX.s.param,
                        "-m", MCScanX.m.param)

  # -- Check that MCScanX is there
  if (!file.exists(MCScanX.tool))
    stop(MCScanX.tool, "doesn't exist\n\t",
         "Make sure the MCScanX program is installed in",
         MCScanX.path,"\n")

  check_gsDir(dir2check = mcscan.dir,
              overwrite.output.dir = overwrite.output.dir)

  gf <- subset(gff, genome %in% genomeIDs)
  setkey(gf, genome, chr)
  gf[,gene.num := paste0("gene",1:nrow(gf))]
  lets <- paste0(letters,letters)[1:length(unique(gf$genome))]
  gf[,genome.num := as.numeric(factor(genome, levels = genomeIDs))]
  gf[,chr.num := as.numeric(as.factor(chr)),
     by = genome]
  gf[,genome.num := paste0(lets[genome.num],chr.num)]
  bl <- subset(blast,
               genome1 %in% genomeIDs &
                 genome2 %in% genomeIDs)
  bl <- merge_gffWithBlast(gff = gf,
                     blast = bl,
                     mirror = T,
                     include.gene.num = T,
                     include.genome.num = T)



  gff.in <- gf[,c("genome.num", "gene.num", "start", "end")]
  blast.in <- bl[,c("gn1", "gn2", "perc.iden", "align.length",
                       "n.mismatch", "n.gapOpen", "q.start", "q.end",
                       "s.start", "s.end", "eval", "score")]
  write.table(gff.in,
              file = file.path(mcscan.dir, "xyz.gff"),
              row.names = F,
              col.names = F,
              quote = F, sep = "\t")
  write.table(blast.in,
              file = file.path(mcscan.dir, "xyz.blast"),
              row.names = F,
              col.names = F,
              quote = F, sep = "\t")
  if(verbose)
    cat("Done!\nRunning MCScanX ... ")
  if (!verbose) {
    com <- paste(MCScanX.tool, mcscan.param,
                 file.path(mcscan.dir, "xyz"), "&> /dev/null")
  } else {
    com <- paste(MCScanX.tool, mcscan.param,
                 file.path(mcscan.dir, "xyz"))
  }
  system(com)

  if(verbose)
    cat("\nParsing MCScanX output ... ")
  gn.ids <- parse_mcs(mcs.file = file.path(mcscan.dir, "xyz"))
  setnames(gn.ids, 1:2, c("gn1","gn2"))

  if("block.id" %in% colnames(bl))
    bl[,block.id := NULL]

  mcs.out <- merge(gn.ids[,c("gn1","gn2")], bl, by = c("gn1", "gn2"))

  mcs.out <- subset(mcs.out, genome1 != genome2)
  mcs.out <- mcs.out[,colnames(blast)[colnames(blast) %in% colnames(mcs.out)], with = F]

  if("gn1" %in% colnames(mcs.out))
    mcs.out[,gn1 := NULL]
  if("gn2" %in% colnames(mcs.out))
    mcs.out[,gn2 := NULL]
  if("genome.num1" %in% colnames(mcs.out))
    mcs.out[,genome.num1 := NULL]
  if("genome.num2" %in% colnames(mcs.out))
    mcs.out[,genome.num2 := NULL]

  if(verbose)
    cat("Done!\n")
  nl <- with(subset(blast, genome1 != genome2), table(genome1, genome2))
  no <- with(mcs.out, table(genome1, genome2))

  if(verbose) {
    cat("Original blast hit counts:\n")
    print(nl)
    cat("\nMCScanX-culled collinear blast hit counts:\n")
    print(no)
  }

  return(mcs.out)
}
