# This is **beta** version 1.0.3 of the GENESPACE R package. **USE WITH CAUTION** 

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

The function `synteny` annotates each unique pairwise synHits file with three columns: "isAnchor", "inBuffer" and "blkID". If blkID is NA and inBuffer == FALSE, that blast hit is not syntenic. Otherwise, the hit is syntenic and may be used for direct syntenic position interpolation if isAnchor == TRUE. These flags are determined through a three step process: First, potential syntenic anchors are extracted from the synHits as !noAnchor, (and if onlyOgAnchors = TRUE, sameOg) and are the topN scoring hits for each gene. topN = ploidy of the alternate genome. These potential anchor hits are fed into MCScanX. Those in collinear blocks are potential syntenic anchors. Second, potential syntenic anchors are clustered into large regions with dbscan where the search radius is the `synBuffer`. Potential syntenic anchors are then re-called by re-running MCScanX on the anchors within a region. Finally, syntenic anchors are clustered into blocks using dbscan.

#### 3.5 Mapping syntenic positions across genomes


#### 3.6 The synHits output format


#### 3.7 The pangenome output format


#### 3.8 The combBed.txt file



The only inputs required for genespace are a directory containing matched annotations in /peptide and /bed subdirectories and the path to the MCScanX installation.

```
gpar <- init_genespace(
    wd = "/path/to/GENESPACE/workingDir", 
    path2mcscanx = "/path/to/MCScanX")
```

This initializes a run where each genome is named following the file names in the peptide and /bed subdirectories and assumes that all genomes are haploid. 

### 5.2 Specifying custom paths to diamond and orthofinder

Installation of orthofinder (and diamond) represents one of the most common issues getting GENESPACE up and running. If you cannot install into a conda environment and open R from there (or put orthofinder in the path), you will to specify a custom path to the orthofinder and diamond program installations through `init_genespace(..., path2orthofinder = "path/to/orthofinder", path2diamond = "path/to/diamond")`. 

### 5.3 Specifying ploidy

Genome ploidy represents the most common required customization. If all of your genomes have the same ploidy, they can be set by specifying a single value. For example for all diploid genomes (likely representing inbred tetraploid genomes) `init_genespace(..., ploidy = 2)`. 

However, you may only have a couple higher ploidy genomes. In this case, its preferable to specify both the genomeIDs and ploidy to ensure you have the right ploidy matching the right genome. For example: `init_genespace(..., ploidy = c(1, 2, 1), genomeIDs = c("species1", "species2", "species3"))`, where two genomes are haploid and one is an outbred diploid. 

#### 3.4 Main parameters to `init_genespace`

There are many ways to customize the GENESPACE run. Check the helpfile for `init_genespace` for details about them all. Here are the other most common:

- `useHOGs` - should hierarchical orthogroups be used instead of global orthogroups. HOGs tend to split arrays and better reflect orthologs. Global orthogroups better capture paralogous regions. Default (NA) chooses whichever matches the ploidy of the genomes best. 
- `outgroup` or `ignoreTheseGenomes` - genomeIDs that should only be used in the orthofinder run but not considered for synteny etc. If not specified, all genomes are used. `outgroup` is a bit of a missnomer and the argument is kept here for now for better backwards compatibility to <v1. 
- `blkSize` - default is 5. This is the minimum number of anchored blast hits needed to call a syntenic block. 
- `nGaps` - default is 5. Passed to `MCScanX -m` as the maximum number of gaps in a syntenic block. 
- `synBuff` - default is 100. This is the maximum euclidean distance (in gene-rank order units) from a syntenic anchor for a gene to be considered in the syntenic buffer. Also the minimum distance between two genes to split a tandem array into two clusters.   
- `nCores` - the number of parallel processes to run
- `orthofinderInBlk` - should OrthoFinder be re-run within each large syntenic region? Default is TRUE if any genome ploidy > 1. 





  

### 1.3 If you are starting from scratch ... 

Previously GENESPACE operated on raw gff and peptide annotations. The internal conversion to a simpler format caused most of the reported errors. To simplify this and reduce errors, GENESPACE now only takes processed .bed and .fasta annotations. See below, but in short, these are standard 4-column bed files where the 4th column in the gene ID, which exactly matches the header in the fasta file. 

Before running GENESPACE you need to either build these files yourself (and store them in /bed and /peptide) subdirectories. Alternatively, with the right parameterization, `parse_annotations()` can do the formatting for you.

Below is an example of how to run the full GENESPACE pipeline between human and chicken genomes from NCBI:

First off, make a directory to store the raw genome annotation files

```
genomeRepo <- file.path("~/Downloads/testGenespace/rawRepo")
if(!dir.exists(genomeRepo))
  dir.create(genomeRepo)
dir.create(file.path(genomeRepo, "human"))
dir.create(file.path(genomeRepo, "chicken"))
```

Then download human annotation from NCBI

```
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz",
  destfile = file.path(genomeRepo, "human", "GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
  destfile = file.path(genomeRepo, "human", "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"))
```

And chicken from NCBI

```
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz",
  destfile = file.path(genomeRepo, "chicken", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz",
  destfile = file.path(genomeRepo, "chicken", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz"))
```

Then parse the annotations into the simplified GENESPACE bed and fasta. See 4.3 (automatic parsing of annotations) for more details. 

```
wd <- "~/Downloads/testGenespace/verts"
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo, 
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi", 
  genespaceWd = wd)
```

Now, we can initialize the GENESPACE run

```
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = "/path/to/MCScanX/")
```

Finally, with all the parameters set and the annotations formatted, we can run the entire GENESPACE pipe with one function call:

```
gpar <- run_genespace(gsParam = gpar)
```





### 3.3 Open an interactive R session

You can open R using the the GUI, Rstudio, or as an interactive session in the terminal. If you are using a conda environment or specifying the path to orthofinder in the `$PATH`, you need to enter an R interactive session from that shell environment, which can be done by either calling `R` (for command line interface) `open -na rstudio` if using Rstudio. 

### 3.4 Install GENESPACE R package

Once in R, GENESPACE can be installed directly from github via:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")
```

### 3.5 Install R dependencies

The above command will install the CRAN-sourced dependencies (`data.table`, `dbscan` and `R.utils`). The bioconductor dependencies (`rtracklayer` and `Biostrings`) need to be installed separately via:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
```



## 5. Setting and checking GENESPACE parameters

### 5.1 Initalizing the GENESPACE run with defaults

