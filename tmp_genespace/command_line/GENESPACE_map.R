#!/usr/bin/env Rscript


# -- Set up the command line arguments
###########
suppressPackageStartupMessages(library("argparser"))

tmp.dir = "/Users/jlovell/Desktop/Pvirgatum_v5genespace/tmp"
results.dir = "/Users/jlovell/Desktop/Pvirgatum_v5genespace/results"
blast.dir = "/Users/jlovell/Desktop/Pvirgatum_v5genespace/blast"
input.dir = "/Users/jlovell/Desktop/Pvirgatum_v5genespace/genomes"
mcscan.dir = "/Users/jlovell/Desktop/Pvirgatum_v5genespace/mcscanx"
raw_annot.dir = "/Users/jlovell/Desktop/Pvirgatum_v5genespace/raw_annotations"

peptide.dir = file.path(input.dir,"peptide")
cds.dir = file.path(input.dir,"cds")
assembly.dir = file.path(input.dir,"assembly")
gff.dir = file.path(input.dir,"gff")

genomeIDs = c("Pvirgatum",
              "PhalliiHAL",
              "Sviridis",
              "Phallii",
              "Sbicolor",
              "Bdistachyon")
ploidy = c(4,2,2,2,2,2)
abbrevs = c("Pv","Ph","Pa","Sv","Sb","Bd")

p <- arg_parser(description = "GENESPACE step #1 - Generate collinear blocks and gene-mappings")
p <- add_argument(p, arg = "--directory", default="NULL",
                  help="the directory containing the following subfolders: raw_annotations, genomes")
p <- add_argument(p, arg = "--counts", default="matched_counts.csv",
                  help="the name of the transcript.counts file")
p <- add_argument(p, arg = "--lib_info", default="matched_info.csv",
                  help="the library information dataset.")
p <- add_argument(p, arg = "--baseFormula", default="~1",
                  help="The full formula with all factors to test")
p <- add_argument(p, arg = "--testName", default="NULL",
                  help="The name to append to columns of stats")
p <- add_argument(p, arg = "--redFormula",  default="~1",
                  help="the formula that lacks the term to test")
p <- add_argument(p, arg = "--subsetColumn",  default="NULL",
                  help="the name of the column that can be subset. If specified, so must --factor2keep")
p <- add_argument(p, arg = "--libs2keep",  default="NULL",
                  help="a vector of the library ids to be retained")
p <- add_argument(p, arg = "--factors2keep",  default="NULL",
                  help="a vector of factor levels that will be retained. If specified, so must --subsetColumn")
p <- add_argument(p, arg = "--outputfilename",  default="de.stats",
                  help="name of stats file name")

args <- parse_args(p)
