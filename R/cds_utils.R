#' @title selection utility functions
#' @description
#' \code{cds_utils} Functions that allow selection stat calculation in GENESPACE
#' @name syn_utils
#'
#' @param pep.file file path to peptide fasta file, with the sequences to align
#' @param cds.file file path to cds fasta file, with the sequences to align
#' @param tmp.dir file path to temporary directory
#' @param pal2nal.tool file path to the pal2nal tool, including the executable
#' @param msa.clu file path to the clustall-formatted MSA
#' @param cds.msa.fa file path to the MSA in fasta format.
#' @param codeml.output.file file path to the codeml output file
#' @param codeml.msa.file paml-formatted CDS MSA.
#' @param codeml.cntr.file file path to codeml control file
#' @param mafft.params parameters to pass to mafft
#' @param base1 first base (ATCG)
#' @param base2 second  base (ATCG)
#' @param cds.msa dna StringSet object containing the cds multiple alignment
#' @param round2digits integer, how many places to round stats to.
#' @param ... not currently in use
#'
#' @note \code{cds_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{cds_utils} returns its own arguments.
#'
#' @examples
#' \dontrun{
#' calc_selectionStats(pep.file = "~/Downloads/test.pep.fa",
#'     cds.file = "~/Downloads/test.cds.fa",
#'     tmp.dir = "/Users/jlovell/Downloads",
#'     pal2nal.tool = "/Users/jlovell/Documents/comparative_genomics/programs/pal2nal.v14/pal2nal.pl")
#' }
#'
#'
#' @title Calculate selection stats
#' @description
#' \code{calc_selectionStats} Do alignments and calculate selection stats for
#' a set of sequences
#' @rdname cds_utils
#' @import data.table
#' @export
calc_selectionStats <- function(pep.file,
                                cds.file,
                                geneIDs = NULL,
                                cds.dir = NULL,
                                peptide.dir = NULL,
                                tmp.dir,
                                codeml.msa.file = file.path(tmp.dir,"codon.aln"),
                                msa.clu = file.path(tmp.dir,"msa.clu"),
                                cds.msa.fa = file.path(tmp.dir,"msa.fa"),
                                mafft.params = "--retree 1 --quiet",
                                pal2nal.tool){

  # -- Do the alignments, returning both fasta and paml formats
  align.files <- align_cds(pep.file = pep.file,
                           cds.file = cds.file,
                           tmp.dir = tmp.dir,
                           codeml.msa.file = codeml.msa.file,
                           msa.clu = msa.clu,
                           cds.msa.fa = cds.msa.fa,
                           mafft.params = mafft.params,
                           pal2nal.tool = pal2nal.tool)

  # -- Calculate codeml stats
  codeml.stats <- calc_kaks(codeml.msa.file = align.files$codeml.msa.file,
                            tmp.dir = tmp.dir)

  # -- Calculate 4-fold site stats
  fdtv.stats <- pairwise_4dtv(cds.msa.fa = align.files$cds.msa.fa)

  # -- return data set. If there are codeml results, merge, otherwise populate with NAs.
  if(is.na(codeml.stats)[1] & length(codeml.stats) == 1){
    out.stats <- fdtv.stats
    for(i in c("t","S","N","KaKs","Ka","Ks")) out.stats[,i] <- NA
  }else{
    out.stats <- merge(fdtv.stats, codeml.stats, by = c("id1","id2"))
  }
  return(list(stats = out.stats, files = align.files))
}


#' @title Calculate stats from codeml
#' @description
#' \code{calc_kaks} Calculate stats from codeml
#' @rdname cds_utils
#' @import data.table
#' @export
calc_kaks <- function(codeml.msa.file, tmp.dir){

  codeml.cntr.file <- file.path(tmp.dir,"tmp.cdml.ctl")
  codeml.output.file <- file.path(tmp.dir,"tmp.cdml")
  make_codemlCtl(codeml.msa.file = codeml.msa.file,
                 codeml.output.file = codeml.output.file,
                 codeml.cntr.file = codeml.cntr.file)

  run_cmd <- paste("codeml", basename(codeml.cntr.file))
  rCdml <- system(run_cmd, intern=TRUE)

  if(strsplit(rCdml," ")[[1]][1] == "CODONML"){
    kaks <- parse_Codeml(codeml.output.file = codeml.output.file)
  }else{
    kaks <- NA
  }
  return(kaks)
}

#' @title Parse codeml output
#' @description
#' \code{parse_Codeml} Parse codeml output to obtain dN/dS and other values
#' @rdname cds_utils
#' @import data.table
#' @export
parse_Codeml <- function(codeml.output.file){

  # -- Take all rows with stats and toss out non-numeric observations, splitting by =
  parse_codemlStats <- function(x)
    rbindlist(lapply(x, function(y){
      tmp <- data.table(t(as.numeric(sapply(strsplit(y,"=")[[1]][-1],  function(z)
        gsub("[^0-9\\.]", "", z)))))
      setnames(tmp, c("t","S","N","KaKs","Ka","Ks"))
      return(tmp)
    }))

  # -- Take all rows with gene IDs, toss parentheses and split by whitespace
  parse_codemlIDs <- function(x)
    rbindlist(lapply(x, function(y){
      tmp <- data.table(t(sapply(strsplit(y," ")[[1]][c(2,5)], function(z)
        gsub(")","",gsub("(", "", z, fixed = T), fixed = T))))
      setnames(tmp, c("id1","id2"))
      return(tmp)
    }))

  # -- read in the results
  cdmlInfo <- readLines(codeml.output.file)

  # -- find the start of the pairwise statistics
  wh.start <- grep("pairwise comparison, codon frequencies", cdmlInfo, fixed = T)
  s <- cdmlInfo[wh.start:length(cdmlInfo)]

  # -- find and parse the gene IDs
  s.ids <- parse_codemlIDs(s[grep("...", s, fixed = T)])
  s.dat <- parse_codemlStats(s[grep("^t=", s)])
  return(cbind(s.ids, s.dat))
}

#' @title make codeml control file
#' @description
#' \code{make_codemlCtl} make codeml control file
#' @rdname cds_utils
#' @export
make_codemlCtl <- function(codeml.msa.file,
                           codeml.output.file,
                           codeml.cntr.file){
  codeml_template.ctl <- "seqfile = input.msa\n    outfile = output.cdml\n       noisy = 0\n     verbose = 0\n     runmode = -2\n   cleandata = 1\n     seqtype = 1\n   CodonFreq = 2\n       model = 2\n     NSsites = 0\n       icode = 0\n       Mgene = 0\n   fix_kappa = 0\n       kappa = 2\n   fix_omega = 0\n       omega = 1\n   fix_alpha = 1\n       alpha = .0\n      Malpha = 0\n       ncatG = 4\n       clock = 0\n       getSE = 0\nRateAncestor = 0\n      method = 0"
  codeml_template.ctl <- sub("input.msa", codeml.msa.file, codeml_template.ctl)
  codeml_template.ctl <- sub("output.cdml", codeml.output.file, codeml_template.ctl)
  writeLines(codeml_template.ctl, con = codeml.cntr.file)
}

#' @title conduct multiple sequence alignments
#' @description
#' \code{align_cds} call mafft and pal2nal to conduct multiple sequence alignments
#' @rdname cds_utils
#' @export
align_cds <- function(pep.file,
                     cds.file,
                     tmp.dir,
                     pal2nal.tool,
                     codeml.msa.file = file.path(tmp.dir,"codon.aln"),
                     msa.clu = file.path(tmp.dir,"msa.clu"),
                     cds.msa.fa = file.path(tmp.dir,"msa.fa"),
                     mafft.params = "--retree 1 --quiet"){

  # -- run msa
  mafft.com <- sprintf("mafft %s --clustalout %s > %s",
                 mafft.params,
                 pep.file,
                 msa.clu)
  msatmp <- system(mafft.com)

  # -- run pal2nal
  pal.com <- sprintf("perl %s %s %s -output fasta -nogap > %s",
                     pal2nal.tool, msa.clu, cds.file, cds.msa.fa)
  paltmp <- system(pal.com)

  pal.com <- sprintf("perl %s %s %s -output paml -nogap > %s",
                 pal2nal.tool, msa.clu, cds.file, codeml.msa.file)
  paltmp <- system(pal.com)

  return(list(cds.file = cds.file,
              pep.file = pep.file,
              msa.clu = msa.clu,
              cds.msa.fa = cds.msa.fa,
              codeml.msa.file = codeml.msa.file))
}

#' @title calculate 4-fold stats on msa
#' @description
#' \code{pairwise_4dtv} calculate 4-fold stats on msa
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readDNAStringSet
#' @export
pairwise_4dtv <- function(cds.msa.fa){

  # -- Read in the fasta formatted MSA
  cds.msa <- readDNAStringSet(cds.msa.fa)

  # -- Find all combinations of gene IDs
  eg <- data.table(expand.grid(names(cds.msa), names(cds.msa), stringsAsFactors = F))
  setnames(eg, c("id1","id2"))
  eg <- eg[with(eg, id1 != id2),]

  # -- For each pairwise combination, subset the alignment and calculate stats
  out <- rbindlist(apply(eg, 1, function(x){
    cds.tmp <- cds.msa[x]
    stats.4dtv <- data.table(id1 = x[1],
                             id2 = x[2],
                             calc_4dtv(cds.tmp))
    return(stats.4dtv)
  }))
  return(out)
}

#' @title calculate 4-fold stats
#' @description
#' \code{calc_4dtv} calculate 4-fold stats
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readDNAMultipleAlignment DNAStringSet
#' @importFrom IRanges successiveViews
#' @export
calc_4dtv <- function(cds.msa,
                      round2digits = 5){

  if (class(cds.msa) != "DNAStringSet")
    stop("cds.msa must be a DNAStringSet object\n")

  if (length(cds.msa) != 2)
    stop("MSA must have exactly two sequences\n")

  #######################################################
  fd.AA <- rep(c("L","V","S","P","T","A","R","G"), each = 4)
  fd.DNA <- c("CTT", "CTC", "CTA", "CTG",
              "GTT", "GTC", "GTA", "GTG",
              "TCT", "TCC", "TCA", "TCG",
              "CCT", "CCC", "CCA", "CCG",
              "ACT", "ACC", "ACA", "ACG",
              "GCT", "GCC", "GCA", "GCG",
              "CGT", "CGC", "CGA", "CGG",
              "GGT", "GGC", "GGA", "GGG")
  #######################################################

  # -- function to convert ss to character vector of codons
  ss2codChar <- function(x)
    as.character(unlist(data.frame(
      successiveViews(x, rep.int(3L, length(x) %/% 3L))
    )))

  # -- convert two CDS sequences into codon data.table
  cod <- data.table(cds1 = ss2codChar(cds.msa[[1]]),
                    cds2 = ss2codChar(cds.msa[[2]]),
                    stringsAsFactors = F)
  cod$base1 <- toupper(substring(cod$cds1, 3, 3))
  cod$base2 <- toupper(substring(cod$cds2, 3, 3))

  # -- subset to 4-fold sites and add in amino acid data
  FD <- subset(cod,
               cds1 %in% fd.DNA &
                 cds2 %in% fd.DNA)
  FD$AA1 <- fd.AA[match(FD$cds1, fd.DNA)]
  FD$AA2 <- fd.AA[match(FD$cds2, fd.DNA)]
  FD <- subset(FD, AA1 == AA2)

  # -- Calculate the whether mutations are transitions and transversions
  FD$Tv <- with(FD, is.transversion(base1, base2))
  FD$Ts <- with(FD, is.transition(base1, base2))

  # -- Number of 4-fold sites
  FDsites <- nrow(FD)

  # -- Number of transitions and transversions
  Tv <- sum(FD$Tv)
  Ts <- sum(FD$Ts)

  # -- ratio of transitions and transversions
  FDTV <- round(Tv/FDsites, round2digits)
  FDTS <- round(Ts/FDsites, round2digits)

  # -- Corrections of FDTV
  FDTVc <- ifelse(FDTV < 0.5,
                  (-0.5) * log(1 - 2 * FDTV),
                  99.0)
  FDTVc <- round(FDTVc, round2digits)

  out <- data.table(FDsites = FDsites,
                    Ts = Ts,
                    Tv = Tv,
                    FDTS = FDTS,
                    FDTV = FDTV,
                    FDTVc = FDTVc)

  return(out)
}

#' @title Check for base transition
#' @description
#' \code{is.transition} Check for base transition
#' @rdname cds_utils
#' @import data.table
#' @export
is.transition <- function(base1,
                          base2){

  transition <- list(A = c("G"),
                     C = c("T"),
                     G = c("A"),
                     T = c("C"))

  if (length(base1) != length(base2) || length(base1) == 0L)
    stop("base1 and base2 must be character vectors of equal length")
  out <- sapply(seq(base1), function(x){
    return(base1[x] %in% transition[[base2[x]]])
  })
  return(out)
}

#' @title Check for base transversion
#' @description
#' \code{is.transversion} Check for base transversion
#' @rdname cds_utils
#' @import data.table
#' @export
is.transversion <- function(base1,
                            base2){

  transversion <- list(A = c("T","C"),
                       C = c("A","G"),
                       G = c("T","C"),
                       T = c("A","G"))

  if (length(base1) != length(base2) || length(base1) == 0L)
    stop("base1 and base2 must be character vectors of equal length")
  out <- sapply(seq(base1), function(x){
    return(base1[x] %in% transversion[[base2[x]]])
  })
  return(out)
}

#' @title Converts MSA file to AXT format
#' @description
#' \code{write_axt} Converts MSA file to AXT format
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readDNAMultipleAlignment DNAStringSet
#' @export
write_axt <- function(cds.msa.fa,
                      tmp.dir = NULL){

  cds.msa <- DNAStringSet(
    readDNAMultipleAlignment(cds.msa.fa, format = "fasta"))

  axt <- paste(unlist(lapply(seq(cds.msa), function(x) {
    toString(cds.msa[[x]])
  })), collapse = "\n")

  axt <- paste0(paste(names(cds.msa), collapse = "-"),"\n", axt)

  if (is.null(tmp.dir))
    tmp.dir <- dirname(cds.msa.fa)

  out_axt <- sprintf("%s/%s.axt", tmp.dir, basename(cds.msa.fa))
  writeLines(axt, con = out_axt)
  return(out_axt)
}
