# Pipeline

## Overview:

Here we extend basic annotation-annotation alignments by searching for un-annotated sequences
that appear to be orthologous to CDS regions annotated in related species. To accomplish this goal,
we constrain our search within collinear blocks that share an evolutionary origin between species. We search
for, then add the un-annotated sequences into the fasta files, then fina

We start by employing the standard blast-orthofinder workflow to find pairwise orthology networks between species.

# 1. Setup
### Folders for temporary, results, blast and input files
tmp.dir = "/Users/jlovell/Documents/test_compareGenespace/tmp"
results.dir = "/Users/jlovell/Documents/test_compareGenespace/results"
blast.dir = "/Users/jlovell/Documents/test_compareGenespace/blast"
input.dir = "/Users/jlovell/Documents/test_compareGenespace/genomes"
mcscan.dir = "/Users/jlovell/Documents/test_compareGenespace/mcscanx"

### The input directory must contain the following sub-directories:
peptide.dir = file.path(input.dir,"peptide")
cds.dir = file.path(input.dir,"cds")
assembly.dir = file.path(input.dir,"assembly")
gff.dir = file.path(input.dir,"gff")

### Metadata must include:
genomeIDs = c("PhalliiHAL","Phallii","Sviridis","Sbicolor","Pvirgatum")
ploidy = c(2,2,2,2,4)
abbrevs = c("Ph","Pa","Sv","Sb","Pv")

### The following programs must be installed and placed into the path:
#### "bedtools","MCScanX","Diamond","orthofinder","exonerate"

### Also, ensure that the following R packages are installed:
#### "data.table", "biostrings", "dbscan"


# 2. Run the program, step-by-step
## 2.1 -- Check that everything is right ...
check_environment()

## 2.2 -- Run initial pairwise BLASTs, using orthofinder and the Diamond blast engine
blast.results = runParse_orthofinder(
  peptide.dir = peptide.dir,
  tmp.dir = tmp.dir,
  blast.dir = blast.dir,
  gff.dir = gff.dir,
  ploidy = ploidy,
  min.propMax = .3,
  min.score = 100,
  nmapsPerHaplotype = 1,
  eps.radius = c(100,50,20),
  n.mappingWithinRadius = c(10,10,10),
  runOF = T,
  fasta.pattern = "*.fa",
  verbose = T)

single.genes = find_genesWithoutOrthos(blast.dir, gff.dir)

## 2.3 -- Make preliminary blocks via MCScanX
mcscan.raw = run_MCScanX(
  blast.results = blast.results,
  abbrevs = abbrevs,
  mcscanx.input.dir = mcscan.dir,
  MCScanX.params = "-a -s 5 -m 25")
mcscan.blocks = make_blocks(mcscan.raw)

## 2.4 -- Merge blocks that overlap (2x)
merge1 = merge_overlappingBlocks(
  map = mcscan.blocks$map,
  blk = mcscan.blocks$block,
  buffer = 1.5,
  verbose = T)
merged = merge_overlappingBlocks(
  map = merge1$map,
  blk = merge1$block,
  buffer = 1.5,
  verbose = T)

merged$block$block.id = paste0("block_", merged$block$block.id)

## 2.5 -- Rerun ortho-finder within blocks to define genes without orthologs
ogs.bygroup = run_orthoFinderInBlock(
  blk = merged$block,
  tmp.dir = tmp.dir,
  buffer = 50000,
  peptide.dir = peptide.dir,
  gff.dir = gff.dir,
  ncores = 8,
  onlyParseOF = F)

## 2.6 -- Find un-annotated sequences by mapping with exonerate

blk = merged$block
genes2ignore = single.genes

blat_protein2Block = function(peptide.dir,
                              genome.dir,
                              tmp.dir,
                              results.dir){

}
...

## 2.7 -- Pull unannotated sequences, re-run orthofinder, with these included, extracting alignments.
...

## 2.8 -- Annotate alignments
...


############################################################################
############################################################################
# 1.2 Get metadata organized
## Analyses are orgnaized by the names of the genomes. These need to be a unique, single word.


## We need to know the ploidy of each genome, this lets us decide the number of blast hits to retain


## Also, for some analyses, we need a two-letter abbreviation for each.


############################################################################
############################################################################
# 1.3 Build the output and temporary data structures












