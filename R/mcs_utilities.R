#' @title Orthofinder utility functions
#' @name of_utilities
#' @aliases prep_mcs
#' @aliases parse_mcs
#' @aliases run_mcs
#'
#' @description
#' \code{mcs_utilities} Several utilities functions meant for internal calling of MCScanX
#' @name of_utilities
#' @param blast character, directory containing the gff3 formatted annotation files
#' @param gff character, directory containing the orthofinder output
#' @param mcscan.dir character, gff file name and path
#' @param mcscan.param character, genome identifiers
#' @param silent.mcs character, genome abbreviations
#' @note \code{mcs_utilities} is a generic name for the functions documented.
#' \cr
#' If called, \code{mcs_utilities} returns its own arguments.
#'

#' @title Prep and run MCScanX
#' @description
#' \code{prep_mcs} Prep and run MCScanX
#' @rdname blk_utilities
#' @import data.table
#' @export
prep_mcs <- function(blast,
                     gff,
                     mcscan.dir,
                     mcscan.param,
                     silent.mcs){

  gff[,chr.num := frank(chr, ties.method = "dense"),
      by = genome]
  gff$rank.start = frank(gff, genome, chr, start, ties.method = "random")
  gff$rank.end = frank(gff, genome, chr, end, ties.method = "random")
  gff[,genome.num := frank(genome, ties.method = "dense")]

  lets = paste0(letters,letters)[1:length(unique(gff$genome))]
  gff$genome.abbrev = paste0(lets[gff$genome.num],1)

  gff.in = gff[,c("genome.abbrev","id","rank.start","rank.end")]
  gff.in = gff.in[with(gff.in, order(genome.abbrev, rank.start)),]
  blast.in = blast[,c("id1","id2","perc.iden","align.length",
                      "n.mismatch","n.gapOpen", "q.start", "q.end",
                      "s.start", "s.end","eval","score")]

  write.table(gff.in,
              file = file.path(mcscan.dir,"xyz.gff"),
              row.names = F,
              col.names = F,
              quote = F, sep= "\t")
  write.table(blast.in,
              file = file.path(mcscan.dir,"xyz.blast"),
              row.names = F,
              col.names = F,
              quote = F, sep= "\t")
  if(silent.mcs){
    com <- paste("MCScanX", mcscan.param, file.path(mcscan.dir,"xyz"),"&> /dev/null")
  }else{
    com <- paste("MCScanX", mcscan.param, file.path(mcscan.dir,"xyz"))
  }

  system(com)
  return(file.path(mcscan.dir,"xyz"))
}

#' @title Parse MCScanX output
#' @description
#' \code{parse_mcs} Parse MCScanX output
#' @rdname blk_utilities
#' @import data.table
#' @export
parse_mcs <- function(mcs.file){
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
  out<-data.table(out)
  setkey(out, id1, id2)
  return(out)
}

#' @title Run a single instance of MCScanX
#' @description
#' \code{run_mcs} Run a single instance of MCScanX
#' @rdname blk_utilities
#' @import data.table
#' @export
run_mcs <- function(blast,
                    gff,
                    mcscan.dir,
                    mcscan.param,
                    silent.mcs){
  mcs.file = prep_mcs(blast, gff,
                      mcscan.dir,
                      mcscan.param,
                      silent.mcs = silent.mcs)
  mcs.parsed = parse_mcs(mcs.file)
  blast = data.table(blast)
  setkey(blast, id1, id2)
  map.out <- data.table(merge(mcs.parsed, blast))
  setkey(map.out, block.id, chr1, start1)
  return(map.out)
}
