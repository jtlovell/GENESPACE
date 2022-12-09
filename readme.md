# This is **beta** version 1.0.4 of the GENESPACE R package. **USE WITH CAUTION** 

This is (hopefully) the last beta release of v1 GENESPACE. We expect a full release to the main branch in early December 2022. If you use this beta version, thanks! Your testing will help improve this and future releases. Please report any problems you run into and we will address them ASAP. 

## GENESPACE: an R package for synteny- and orthology-constrained comparative genomics. 

This is GENESPACE v1. GENESPACE is an analytical pipeline that lets you compare syntenic regions across multiple genomes. The manuscript describing GENESPACE is now pubished in eLife. [Find the article here](https://elifesciences.org/articles/78526). Please cite this if you use GENESPACE. **NOTE** There are some significant changes to the structure of GENESPACE in V1 (but few changes to the underlying algorithms). These changes are detailed below in '3.1: Changes to GENESPACE in v1'. 

To run GENESPACE, you need:

1. A valid installation of OrthoFinder, MCScanX and R (see part 1: software installation)
2. A peptide (protein) fasta and a bed-formatted gene model annotation file for 2 or more genomes (see part 2: formatting your annotations)

With these data in hand, you initialize a GENESPACE run, then run the pipeline from R:

```
library(GENESPACE)
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = "~/Desktop/programs/MCScanX/")
gpar <- run_genespace(gsParam = gpar) 
```

`run_genespace` completes the full pipeline, which includes:

1. Tandem array discovery
2. Syntenic block coordinate calculation
3. Synteny-constrained orthogroups
4. Pairwise dotplots
5. Syntenic position interpolation of all genes
6. Pan-genome annotation construction
7. Multi-genome riparian plotting

You can then explore the results using `query_genespace` (see section 5: exploring the results). For example, this lets you pull out genes in regions of interest across all genomes

## Running with example data

We recommend doing a test run from real data before applying GENESPACE to your own datasets. To help facilitate this, run these commands in R (assuming valid installation)

```
# -- specify these parameters so that they are valid paths in your system 
genomeRepo <- "~/Desktop/gs_v1_runs/rawGenomeRepo/"
wd <- "~/Desktop/gs_v1_runs/test4readme2"
path2mcscanx <- "~/Desktop/programs/MCScanX/"

# -- download files to your system from NCBI 
# (also can be accomplished with) download_exampleData()
dir.create(file.path(genomeRepo, "human"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz",
  destfile = file.path(genomeRepo, "human", "GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
  destfile = file.path(genomeRepo, "human", "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"))

dir.create(file.path(genomeRepo, "chicken"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz",
  destfile = file.path(genomeRepo, "chicken", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz",
  destfile = file.path(genomeRepo, "chicken", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz"))

# -- reformat the annotations
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo, 
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi", 
  genespaceWd = wd)

# -- intialize the genespace run
gpar <- init_genespace(
  wd = wd, nCores = 1,
  path2mcscanx = path2mcscanx)

# -- run the genespace pipeline
out <- run_genespace(gpar)
```

## 1. Software installation

Installation for v1 is identical to v0.9.4. Use the same conda (or other) installation environment as before. To install from scratch, you need R, a few 3rd party programs and a few R packages. Here is how to install these:

#### 1.1 Install R

GENESPACE is meant to be run interactively in the R environment for statistical computing. So, you need to have R installed. See [CRAN](https://www.r-project.org/) for the most recent release. 

#### 1.2 Install orthofinder

GENESPACE also requires third-party software that can be installed as follows.

`OrthoFinder` (which includes `DIAMOND2`) is most simply installed via conda (in the shell, not R). 

```
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```

If conda is not available on your machine, you can install orthofinder from a number of other sources. See [orthofinder documentation](https://github.com/davidemms/OrthoFinder) for details.  
Regardless of how `OrthoFinder` is installed, ensure that you have `OrthoFinder` version >= 2.5.4 and `DIAMOND` version >= 2.0.14.152.

#### 1.3 Install MCScanX

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 

#### 1.4 Install GENESPACE

Once the above 3rd party dependencies are installed, get into R. If you made a conda environment, its useful to open R directly from that environment so that OrthoFinder stays in the path. 

```
conda activate orthofinder
open -na rstudio # if using rstudio, otherwise, simply `R`
```

Once in R, the easiest way to install GENESPACE uses the package `devtools` (which may need to be installed separately):

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE@dev")
```

If they are not yet installed, this will install a few dependencies directly (ggplot2, igraph, dbscan, R.utils, parallel). However, you will need to install the bioconductor packages separately:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)
```

## 2. Formatting your annotations

This can be the trickiest part of running GENESPACE. If you are having trouble, open an issue. 

#### 2.1 Convert a previous run and re-run GENESPACE

Earlier versions of GENESPACE used 'simplified gff3' annotation files. To make GENESPACE compatible with existing bioinformatic software, we have opted to use the bed format for v1 onward. If you have a previous GENESPACE run, you can easily convert the input to v1 as follows (this also copies over orthofinder results so those can be re-used too): 

```
wd <- "/place/to/store/new/run"
convert_input2v1(
  v1Dir = wd, 
  existingDir = "/path/to/existing/run")
```

#### 2.2 GENESPACE-readable annotation format

For each genome, GENESPACE needs:

1. bed formatted coordinates of each gene (chr, start, end, name), other fields are allowed, but will be ignored by GENESPACE
2. peptide sequences in fasta format, where the header exactly matches the "name" (4th) bed column

Example bed file:

```
cat /bed/genome1.bed | head -2
chr1  1500  2000  gene1
chr1  2500  2200  gene2
```

Example peptide fasta file: 

```
cat /peptide/genome1.fa | head -4
>gene1
MGASGRGAGEQQSPSSTG
>gene2
MLVMSECKGRDRSPSSSM
```

#### 2.3 Manual parsing of raw annotation files

It's possible to format all files using custom scripts. If you opt for this route, re-name and place the genespace-formatted annotations in your working directory. For example for two genomes ("genome1", "genome2"):

```
/workingDirectory
└─ peptide
    └─ genome1.fa
    └─ genome2.fa
└─ bed
    └─ genome1.bed
    └─ genome2.bed
```

#### 2.4 Naming and storing raw annotation files

Its usually preferable to have a static repository containing the raw versions of all of your genome annotations, lets call this directory `/genomeRepo`. Each unique genome annotation should be stored in its own subdirectory with an informative name. Here, each genome is given a species id, genotype id and genome version. The annotation files of interest are the gff3-formated protein coding gene features and the primary transcript translated cds (peptides). These are the two most common types of annotation files available. Download them from your favorite repositories and place them in a directory structure like below. The file names within each genome directory can be anything as long as they have file names that can be found through regular expression matching (see `parse_annotations(gffString, faString)`).

```
/genomeRepo
└─ species1_genoX_v1.0_NCBI
    └─ peptidesGenoX.fa
    └─ genesGenoX.gff3
└─ species2_genoY_v1.0_NCBI
    └─ peptidesSpecies2.fa
    └─ genes.gff3
└─ species3_genoW_v1.0_phytozome
    └─ peptides.fa
    └─ genesSpecies2v1.gff3
└─ species4_genoZ_v1.0_otherRepo
    └─ peptides.fa
    └─ genes.gff3
  ...
```

#### 2.5 Partially-automated annotation parsing

Assuming the directory structure from 2.4, we want to parse the raw annotations into the directory structure outlined in 2.3. The convenience function `parse_annotations` does this. Lets take the simple case where we have translatedCDS and gene.gff3 files from NCBI. In this case, we can specify the "ncbi" present which will automatically parse these files correctly:

```
parsedPaths <- parse_annotations(
  rawGenomeRepo = "/genomeRepo", 
  genomeDirs = c("species1_genoX_v1.0_NCBI", "species2_genoY_v1.0_NCBI"),
  genomeIDs = c("species1", "species2"),
  gffString = "gff3",
  faString = "fa",
  presets = "ncbi", 
  genespaceWd = "/path/to/GENESPACE/workingDir")
```

This will parse the two NCBI-formatted annotations into matched bed and fasta files and place those in the /bed and /peptide subdirectories of the GENESPACE working directory. Importantly, it will re-name the annotations with the genomeIDs. This allows for informative naming conventions in the raw file repo, but simpler names for plotting and analysis. 

#### 2.6 Non-standard annotation formats

Most repositories have their own annotation formats. You may want to combine annotations from NCBI with those from phytozome (plant genome repo) or non-standard repositories. To do this, you simply re-run parse_annotations with different genomes. For example add a phytozome-formatted annotation:

```
parsedPaths <- parse_annotations(
  rawGenomeRepo = "/genomeRepo", 
  genomeDirs = "species3_genoW_v1.0_phytozome",
  genomeIDs = "species3",
  gffString = "gff3",
  faString = "fa",
  presets = "phytozome",
  genespaceWd = "/path/to/GENESPACE/workingDir")
```

Then add a non-standard annotation where there is only a simple geneID header on the peptide and the gene names are in the "GeneID" ';'-delimited 9th gff column. 

```
parsedPaths <- parse_annotations(
  rawGenomeRepo = "/genomeRepo", 
  genomeDirs = "species4_genoZ_v1.0_otherRepo",
  genomeIDs = "species4",
  gffString = "gff3",
  faString = "fa",
  headerEntryIndex = 1, 
  gffIdColumn = "GeneID",
  genespaceWd = "/path/to/GENESPACE/workingDir")
```

#### 2.7 Additional parameterization

There are lots of parameters that can be combined to parse every annotation format. See the help file for details. The main parameters to get your annotation parsed correctly are: 

  - `gffIdColumn`: the name (or index) of the ';' separated field in the attributes column contains the geneID. 
  - `headerEntryIndex`: the name (or index) of the `headerSep`-separated field that contains the geneID in the fasta header
  - `gffStripText`: character string to strip off the geneID in the gff3 attribute column
  - `headerStripText`: character string to strip off the geneID in the fasta header

Parsing can take some troubleshooting, which can be aided by setting `troubleshoot = TRUE`. This prints the first 10 lines of the raw and parsed gff and fasta headers. For some annotation types, you may need to give a fairly complex regular expression to gffStripText or headerStripText. For example, if the fasta name has some additional unwanted info after the last "." in the name: `headerStripText = '[.][^.]+$'`. 


## 3. Running GENESPACE

#### 3.1. Considerations for parameterization


#### 3.2 One line to run it all


#### 3.3 Running orthofinder separately


#### 


## 3. The GENESPACE pipeline

#### 3.1 Changes to GENESPACE v0.9.4 --> v1.0.2+

**MAJOR CHANGES** GENESPACE is now structured differently. This breaks backwards compatibility to < v1, but was required to reduce the issues with raw data input. 

  1. The user must now specify pre-formatted bed files instead of raw gff3 formatted annotations. This allows GENESPACE to catch issues with data formatting before running orthofinder etc. We are truly sorry to alter the specification and input so much, but it is unavoidable -- we want GENESPACE to integrate with existing bioinformatic pipelines, and too much internal customization was required, which produced hard-to-resolve errors.
  2. Plots are now ALWAYS printed to file and never to the R graphics panel: xy dotplots and riparian plots are now saved as vector pdfs in the new /dotplot and /riparian directories respectively. 
  3. Better organization of results output. Previously all output went into /results. /results now only contains the parsed orthofinder results and annotated bed files. All other results are placed in /syntenicHits, /dotplots, /riparian, and /pangenome. 
  4. All GENESPACE functions are now integrated into `run_genespace()`, which completes all steps internally. While the user can still run each step independently, we no longer provide step-by-step documentation and suggest that all runs use this pipeline.
  5. There is a new parameter: `blkRadius`. This controls the size and precision of syntenic blocks. Smaller values will pick up all minor inversions, but will produce much larger output files and potentially less informative global syntenic plots. Larger values will aggregate inversions smaller than blkRadius, but will also make more visually appealing plots. The previous splitting of 'regions' and 'blocks' is no longer applicable. Accordingly, syntenic hit files no longer have the columns 'regID' and 'lgBlkID'. 
  6. There are two changes to the actual algorithm. First, 'inBuffer' hits can no longer be outside of the syntenic block coordinates. This minor change proved important and allowed us to simplify the pangenome construction step since we were better able to estimate the syntenic positions of all genes. Second, we now handle proximate orthogroups (tandem arrays) a bit differently: large dispersed arrays are flagged and ignored (dramatically improving syntenic accuracy in poorly annotated regions like the pericentromere). 
  
**MINOR CHANGES**

  7. code cleanup and better documentation, both in line and in help files
  8. much better data QC and handing of errors in parameter specification
  9. ggplot2-integrated graphics and R object output for potential downstream interactive viz.
  10. internal position calculation for riparian plots and inferred syntenic positions transferred to a new function `integrate_synteny()`. This speeds up riparian plotting by 10x or more. 
  11. `plot_riparian()` function generalization that allows the user to specify a block coordinate matrix and nothing else (allowing plotting of data from other packages). 
  12. replace `orthofinderMethod = "fast"` with `onewayOF = TRUE` to match the new `orthofinder -1` specification. 
  13. a new function `query_genespace()`, allows positional or orthology based queries of pan-genome, pairwise hits, or raw orthogroups/orthologues. 
  14. Use of phylogenetically hierarchical orthogroups (HOGs) instead of default orthogroups (in most cases) to be inline with orthofinder best practices

#### 3.2 Running orthofinder

By default, GENESPACE runs the full OrthoFinder program internally from R. If an orthofinder run is already available (and uses exactly the same genomeIDs as the GENESPACE run), GENESPACE will detect this and not run OrthoFinder. 

For very large runs, we highly recommend running OrthoFinder outside of R. You can do this by using `init_genespace(..., path2orthofinder = NA)`. When you use `run_genespace` with the resulting gs parameters, GENESPACE will prepare the input files and print the OrthoFinder command to the console before killing the run. You can then run that command in your shell or on a server separately. 

Once orthofinder is run, you can point to the output directory via `init_genespace(..., rawOrthofinderDir = "/path/to/run")`, or just put it in the /orthofinder directory in your GENESPACE wd. GENESPACE will find the necessary results, QC them, and if all looks good, move them to the /results directory and begin the pipeline. 

#### 3.3 Aggregating orthofinder information

The bed files are concatenated and a genomeID column is added so that all gene position is in one file: /results/combBed.txt. OrthoFinder ofIDs (unique gene identifiers), OGs (orthogroup.tsv) and HOGs (phylogenetically hierarchical orthogroups, N0.tsv) are parsed and added as new columns to the combBed file. This provides all the raw information needed for all further steps. 

Problematic genes and chromosomes are flagged by two rules: chromosomes must have > `blkSize` (default = 5) unique orthogroups. Orthogroups must hit < `ploidy * maxOgPlaces` (default maxOgPlaces = 8) unique positions in the genome (separated by > `synBuffer`, default = 150). Genes in either the problematic orthogroups or chromosomes are flagged as `noAnchor`, meaning they will still be considered for the pangenome steps, but cannot be used as syntenic anchors or to infer sytenic positions across genomes. 

With HOGs (or OGs if `useHOGs = FALSE`) in hand, tandem arrays are calculated as follows: 1) any OG with multiple members on a single genome-by-chromosome combination are 'potential tandem arrays'; 2) potential tandem arrays with a gap between members > synBuffer are split into their own clusters; 3) clusters with > 1 member are tandem arrays; 4) genes in the most physically central position in the array are 'arrayReps'. Gene order is re-calculated for just the array reps and steps 1-4 are run iteratively until convergence. Genes that are not array representatives are also falgged as `noAnchor`. 

The information from the combBed file are added to the blast files to produce annotated blast files, aka synHits. Of particular importance is whether both the query and target genes are in the same orthogroup "sameOg" and are not flagged as `noAnchor`. Reciprocal query/target blast hits are aggregated into a single file where the larger genome is the query and the smaller the target. 

#### 3.4 Calculating pairwise synteny

The information stored in the combined/annotated bed file is merged with the blast results to make synHits matrices. The function `synteny` annotates each unique pairwise synHits file with three columns: "isAnchor", "inBuffer" and "blkID". If blkID is NA and inBuffer == FALSE, that blast hit is not syntenic. Otherwise, the hit is syntenic and may be used for direct syntenic position interpolation if isAnchor == TRUE. 

These flags are determined through a three step process. First, potential syntenic anchors are extracted from the synHits as !noAnchor, (and if onlyOgAnchors = TRUE, sameOg) and are the topN scoring hits for each gene. topN = ploidy of the alternate genome. These potential anchor hits are fed into MCScanX. Those in collinear blocks are potential syntenic anchors. Second, potential syntenic anchors are clustered into large regions with dbscan where the search radius is the `synBuffer`. Potential syntenic anchors are then re-called by re-running MCScanX on the anchors within a region. Finally, syntenic anchors are clustered into blocks using dbscan.

#### 3.5 Mapping syntenic positions across genomes

## 4. Output file formats

GENESPACE produces four results stored as text files. 

#### 4.1 The synHits output format

The main output of GENESPACE are the annotated blast hits, called synHits. These files are pairwise between two genomes (or within a genome) and collapse reciprocal hits into a single file (to save space + time). 

The query (genome 1) and target (genome 2) gene position information are taken directly from the combBed.txt file. These data, which have names corresponding exactly to the combBed.txt header, fill the first 16 columns of the synBlast.txt.gz files. The next 10 columns are the standard blast8 fields (minus the query and target IDs). Finally, there are 6 columns added by GENESPACE synteny pipeline: 

1. 'sameOg': TRUE/FALSE, are the query and target genes in the same orthogroup (defined by the og column in the combBed file)
2. 'noAnchor': TRUE/FALSE, is either the query or the target flagged as noAnchor in the combBed file?
3. 'isAnchor': TRUE/FALSE, is this blast hit a syntenic anchor, and thereby a valid coordinate to connect positions between genomes?
4. 'inBuffer': TRUE/FALSE, is this blast hit within a buffered radius of a syntenic anchor?
5. 'blkID': unique syntenic block identifier
6. 'sameInblkOg': like sameOg, but just for inBlkOrthofinder

#### 4.2 The pangenome output format

The function `pangenome` converts orthogroup and position data into a single matrix. This is stored in the /pangenome directory. This is the raw data that goes into the pangenome-annotation output (position x genome matrix). The refPangenomeAnnot.txt file contains 11 columns:

1. 'pgID': the unique position-by-orthogroup combination identifier. Each row in the pangenome-annotation is a unique value of this field. 
2. 'pgGenome': the reference genome
3. 'pgChr': the chromosome (syntenic or actual) on the reference genome
4. 'pgOrd': the gene-rank order position (interpolated or actual) on the reference genome
5. 'genome': genome of the focal gene (id column)
6. 'og': orthogroup ID taken from the combBed.txt file
7. 'isRep': TRUE/FALSE, is this gene the representative gene for this pgID?
8. 'ofID': orthofinder gene ID
9. 'id': bed-specified gene ID
10. 'isNSOrtho': TRUE/FALSE, is this gene a non-syntenic ortholog to the representative gene?  If so, flag with "*" in the pangenome-annotation. 
11. 'isArrayRep': TRUE/FALSE, is this gene a non-representative tandem array member? If so, flag with "+" in the pangenome-annotation. 

#### 4.3 The combBed.txt file

Here are the first two lines of the combBed.txt file for the human-chicken (shown just whitespace delimited here).

```
chr start end id ofID pepLen ord genome arrayID isArrayRep globOG globHOG synOG inblkOG noAnchor og
1 7227 96035 TBCE 1_17482 578 922 human NoArr17298 TRUE OG0011962 NA 7295 NA FALSE 7295
1 7273 18867 MTX1 1_11094 466 923 human NoArr15947 TRUE OG0002010 NA 1124 NA FALSE 1124
```

Column definitions:

1. 'chr': chromosome identifier copied from the bed file
2. 'start': gene start position copied from the bed file
3. 'end': gene end position copied from the bed file
4. 'id': gene identifier copied from the bed file
5. 'ofID': unique orthofinder ID taken from SequenceIDs.txt
6. 'pepLen': number of amino acids in the predicted peptide
7. 'ord': gene rank order position along the entire genome
8. 'genome': unique genome ID, taken from the bed file name
9. 'arrayID': tandem array identifier, if the field begins with "NoArr", then this gene is not part of a tandem array
10. 'isArrayRep': TRUE/FALSE specifying whether a gene is the representative (closest physically to the median position of the array)
11. 'globOG': global orthogroup ID, taken from orthogroups.tsv
12. 'globHOG': global phylogenetically hierarchical orthogroup, taken from N0.tsv
13. 'synOG': syntenic OG (or HOG if useHOGs = TRUE) where only genes that are members of syntenic blast hits are allowed to be in the same subgraph (OG)
14. 'inblkOG': if inBlkOrthofinder = TRUE, this is the re-calculated HOG within each syntenic block
15. 'noAnchor': TRUE/FALSE whether the gene can be considered for synteny
16. 'og': changes depending on the stage of the run. Initially it is the globHOG or globOG. Then after synteny, it is populated with synOG. Lastly, may be populated with inblkOG. 

#### 4.4 The syntenic block coordinates

This file, stored in /results/blkCoords.txt, contains pairwise block coordinates for each pair of genomes (genome1, genome2; 1 or 2 appended to the column name indicates that this data is associated with either genome1 or genome2 respectively). For each genome, there are 15 columns:

1. 'chr': chromosome ID
2. 'minBp': lowest basepair (bp) position
3. 'maxBp': greatest basepair (bp) position
4. 'minOrd': lowest gene rank order (ord) position
5. 'maxOrd': highest gene rank order (ord) position
6. 'minGene': lowest position gene ID in the block
7. 'maxGene': higest position gene ID in the block
8. 'nHits': number of hits
9. 'startBp': begining bp position (taking into account orientation)
10. 'endBp': end bp position (taking into account orientation)
11. 'startOrd': begining ord position (taking into account orientation)
12. 'endOrd': end ord position  (taking into account orientation)
13. 'firstGene': first gene in the block (taking into account orientation)
14. 'lastGene': last gene in the block (taking into account orientation)
