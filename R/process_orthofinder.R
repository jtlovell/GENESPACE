#' @title Prepare orthofinder for block construction.
#' @description
#' \code{process_orthofinder} Pipe to process the results of Orthofinder.
#'
#' @param gff.dir character, directory containing the gff3 formatted annotation files
#' @param blast.dir character, directory containing the orthofinder output
#' @param mcscanx.input.dir character, directory where MCScanX temporary files should be stored
#' @param genomeIDs character, genome identifiers
#' @param MCScanX.param character, parameters to be passed to MCScanX
#' @param n.mappingWithinRadius numeric, number of hits required to be in the radius
#' @param eps.radius numeric, size of the radius
#' @param pairs.only logical, should only pairs of hits in orthofinder output be retained
#' @param str2drop character, string in attribute column of gff file to be dropped
#' @param str2parse character, string in attribute column of gff file to use as the separator
#' @param whichAttr numeric, which attribute should be returned in the
#' gff attribute column
#' @param verbose logical, should updates be reported?
#'
#' @details This is the primary blast-processing function of GENESPACE. For additional
#' information, see documentation for the utilities - `?of_utilites`. Here, we
#' run the following:
#'
#' \enumerate{
#'    \item{**Import gff annotation**: calls the internal function `import_gff`,
#'    which subsets the annotations to just entries where the third column
#'    is 'gene'. It then parses the 9th 'attribute' column in the gff, to the
#'    `whichAttr` element of the `str2parse` separated attribute entry. Then,
#'    drops the `str2drop` string from the parsed attribute. For example,
#'    take the following attribute entry:
#'    'ID=Bradi1g00200.v3.1;Name=Bradi1g00200;ancestorIdentifier=Bradi1g00200.v2.1'.
#'    In this case, the default will take the 2nd ';" separated entry and drop the
#'    'Name=' string, returning 'Bradi1g00200'.}
#'    \item{**Import orthofinder results**: calls the internal function
#'    `import_ofResults`, which uncompresses blast results, reads species
#'    ID-name dictionary, organizes by the order of entries in the `genomeIDs`
#'    vector and reads the orthofinder geneID-name dictionary. Last, it reads
#'    the orthofinder orthogroup networks and parses them into a simple named
#'    list}
#'    \item{**Reads and parses blast results**: For each reciprocal pairwise mapping,
#'    calls the internal function `import_blast`, which does the following: (1)
#'    reads in raw blast results, (2) concatenates reciprocal blasts, (3) orders
#'    by decreasing bit score, (4) drops duplicated gene pairs, (5) drops any
#'    genes that are not in orthogroups, (6), if `pairs.only`, reduces to
#'    only gene-pairs that are in the same orthogroups}
#'    \item{**Subsets to clustered hits**: Searches for physically proximate
#'    clusters of blast hits by both 2d density-based clutering (DBScan) and
#'    blast-hit runs (MCScanX). DBScan is run first, retaining hits where
#'    there are `n.mappingWithinRadius` hits within the specified `eps.radius`.
#'    By default, these parameters become increasingly stringent, thereby
#'    first dropping pure noise points, then dropping sparse-hitting regions
#'    (e.g. ancient duplications). Finally, the data is processed by MCScanX.
#'    The `MCScanX.param` specified should be very lax. By default, we require
#'    blocks to have 10 genes, but allow for 25 gaps in the alignments.}
#' }
#'
#' @return Returns a list of length 3, containing:
#'  \enumerate{
#'    \item{`gff` - The concatenated and parsed annotations}
#'    \item{`ortho.info` - A list of length 4, containing the
#'    orthogroups, species mappings, species indices and gene
#'    indices from orthofinder}
#'    \item{`blast` - A data.table with concatenated blast results.}
#' }
#' @import data.table
#' @export
process_orthofinder <- function(gff.dir,
                               genomeIDs,
                               blast.dir,
                               mcscan.dir,
                               pairs.only = T,
                               str2drop = "Name=",
                               str2parse = ";",
                               whichAttr = 2,
                               n.mappingWithinRadius = c(2,2,2),
                               eps.radius = c(50,20,10),
                               mcscan.param = "-a -s 2 -m 10 -w 2 -e 1",
                               verbose = T){

  gff <- import_gff(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    verbose = verbose,
    str2drop = str2drop,
    str2parse = str2parse,
    whichAttr = whichAttr)

  of.blast <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = blast.dir,
    verbose = verbose)

  blast <- import_blast(
    species.mappings = of.blast$species.mappings,
    genomeIDs = genomeIDs,
    orthogroups = of.blast$orthogroups,
    gff = gff,
    gene.index = of.blast$gene.index,
    verbose = verbose)

  cull.dbs <- cull_blastByDBS(blast = blast,
                              n.mappingWithinRadius = n.mappingWithinRadius,
                              eps.radius = eps.radius,verbose = T)

  cull.mcs <- pipe_mcs(blast = cull.dbs,
                       gff = gff,
                       mcscan.dir = mcscan.dir,
                       mcscan.param = mcscan.param)

  return(list(gff = gff,
              ortho.info = of.blast,
              blast = cull.mcs))
}
