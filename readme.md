---
title: "A guide to GENESPACE"
author: "JT Lovell"
date: "5/13/2021"
output:
  html_document:
    df_print: paged
    toc: TRUE
---

# 0. Overview

GENESPACE is a comparative genomics framework implemented in the R environment for statistical computing. The premise is that, when analyzing high-quality genome assemblies and annotations, we can improve the confidence of evolutionary inference by combining two sources of evidence for common ancestry: synteny/collinearity of gene order on chromosomes and coding sequence similarity (homology). In addition to providing a second line of evidence beyond sequence similarity, synteny offers a mechanism by which gene order can be predicted across genomes, even when presence-absence variation means that not all genomes are represented in every orthology network ('orthogroup'). Further, since all inference is constrained within syntenic regions, genespace is agnostic to ploidy, segmental duplicates, inversions and other chromosomal complexities. 

GENESPACE outputs a synteny-constrained and -anchored orthogroup pan-genome annotation among multiple genomes. This simple text file allows for extraction and exploration of regional (e.g. a QTL interval) gene-level variation, a necessary step to integrate comparative and quantitatve genomic goals. 

# 1. Setup

### 1.1 Installing Dependencies

GENESPACE is designed only for use MacOSX and Linux. Testing on OSX 10.16 (Big Sur) and __ (). 
Requires the following third party software to run:

- orthofinder (v2.5.2 or later)
- diamond (v2.0.8.146 or later)
- MCScanX (MCScanX_h)
- R (v4.0.3 or later)

And R package dependencies:

- R/data.table (v1.13.6 or later)
- R/igraph (v1.2.6 or later)
- R/Biostrings (v2.58.0 or later)
- R/dbscan (v1.1-5 or later)

Orthofinder and diamond2 can be installed from conda/miniconda: `conda install -c bioconda orthofinder`. MCScanX must be [installed](https://github.com/wyp1125/MCScanX) independently. R can also be installed via conda or from [CRAN](https://cran.r-project.org/mirrors.html). 

### 1.2 Install GENESPACE

To install GENESPACE, enter R, either by typing `R` into the command line or `open -na rstudio` for RStudio. 

If not installed, install the `devtools` package via `install.packages("devtools")`. This will permit one-line installation from the gitlab repo. 

Install genespace: `devtools::install.gitlab("https://code.jgi.doe.gov/plant/genespace-r")`

This should also install the dependencies, `data.table`, `igraph`, and `dbscan`. For full functionality, also install the following:
`Biostrings` and `Rutils`. 

If devtools is not available or you'd prefer, you can install from the provided GENESPACE.tar.gz file, saved in a location where you have execute privileges. If running on the command line, prior to entering R, run: `R cmd install path/to/GENESPACE.tar.gz`. Alternatively, if running GENESPACE from RStudio, open Rstudio once the orthofinder conda environment is activated, then install GENESPACE by clicking on the 'install' icon under the 'packages' tab. Select 'install from package archive file' and navigate to the location of the saved GENESPACE.tar.gz file. 

### 1.3 Input data format

GENESPACE requires two inputs for each genome: primary transcript annotation peptide fasta files and gene feature coordinates. The sequence headers in the fasta must be able to parsed to match the attributes column of the gff3 gene annotation. 

If data comes from phytozome, you want the 'gene.gff3' and 'protein_primaryTranscriptOnly.fa' files. 

If data comes from NCBI, you to get the data from the 'FTP directory for RefSeq assembly' file tree and download the 'genomic.gff.gz' and 'translated_cds.faa.gz' files. 

If getting data from other sources, you may need to ensure that your gene IDs can be properly matched. It is not necessary to have only the primary transcript, as GENESPACE will take the longest uniquely named sequence. However, in many cases the longest transcript is not the primary, so care should be taken. 


# 2. Pipeline overview

### 2.1 Set the parameters for the run

The first step in GENESPACE is to specify all downstream parameters. This approach is preferable to parameterization of individual functions for consistency across calls. Further, most major issues with data structure, software dependencies or directory permissions will be discovered straight away. 

A fully parameterized `setup_genespace` might look like this:

```
gsPars <- setup_genespace(
  genomeIDs = c("Switchgrass", "Rice"),
  orthofinderMethod = "fast",
  outgroup = NULL,
  speciesIDs = c("Pvirgatum", "Osativa"),
  versionIDs = c("v5.1", "v7.0"),
  rawGenomeDir = "/rawGenomeRep",
  ploidy = c(2,1),
  nCores = 4,
  path2orthofinder = "orthofinder", # assuming installed via conda and in env
  path2diamond = "diamond", # assuming installed via conda and in env
  path2mcscanx = "programs/MCScanX",
  gffString = "gene.gff",
  pepString = "pep",
  verbose = TRUE)
```

- `genomeIDs`: what you want to call each genome. These need to be unique even if the underlying data is duplicated. If one of these genomes should be an outgroup and used only for orthogroup inference and not synteny, this genomeID should be specied under `outgroup`. By default, there are no outgroups. 
- `ploidy`: the expected ploidy of each genome. This is the genome ploidy, which does not necesarily match the organisms cytotype (e.g. most diploid species have haploid genomes)
- `speciesIDs`: only used to detect the raw annotation files (see below)
- `versionIDs`: only used to detect the raw annotation files (see below)
- `rawGenomeDir`: a path to a directory containing all raw gzipped annotation files
    - In this directory, make a subdirectory for each species you'd like to run
    - The species directories contain a subdirectory for each unique genome version to use. 
    - Each genome version subdirectory contains a subdirectory 'annotation', which contains the raw genome annotations. For example:
    - gff: genespace_runs/rawGenomeRepo/Pvirgatum/v5.1/annotation/gene.gff3.gz
    - peptide: genespace_runs/rawGenomeRepo/Pvirgatum/v5.1/annotation/pep.fa.gz
    - This is the raw output structure of data downloaded from phytozome. 
    
This is also where the paths for the 3rd party software are checked. If not in the path when R was opened (or from a conda env), the paths need to be specified, for example:

- `path2mcscanx = "/programs/MCScanX"`
- `path2diamond = "/programs/diamond"`
- `path2orthofinder = "/programs/orthofinder"`

### 2.2 Convert to GENESPACE file format

The basic gff and peptide header conversions are accomplished by `match_gffFasta` or the more user-friendly `parse_ncbi` and `parse_phytozome` functions. 

This is how `match_gffFasta` can be parameterized. It should allow for pretty much any format of fasta and gff3.
```
match_gffFasta(
  genomeIDs, # genomeIDs to parse
  path2rawGff3, # where the raw gff of these genome are
  path2rawPeptide, # where the raw peptide fasta for these genomes are
  gffEntryType = "gene", # column 3 string to match
  gffIdColumn = "Name", # attribute column ; separated field to use
  headerEntryIndex = 4, # field to use in the fasta header
  headerSep = " ", # field separator i the fasta header
  headerStripText = "locus=", # text to drop in the fasta header
  troubleshoot = TRUE) # print the pre and post-parsed head of the files
```

The different parsers can be combined by specifying `overwrite = TRUE`. For example, 
to parse genomes from NCBI, phytozome and a third party with a different format:

```
parse_phyotozome(gsParam = gsp, genomeIDs = "Rice") # gsp from setup_genespace
parse_ncbi(gsParam = gsp, genomeIDs = "Zmays", overwrite = T)
with(gsp, match_gffFasta(
  genomeIDs = genomeIDs[1],
  path2rawGff3 = gffRaw[1],
  headerEntryIndex = 3,
  headerStripText = "locus=",
  path2rawPeptide = peptideRaw[1],
  overwrite = T, troubleshoot = T))
  
# running without overwrite = TRUE after parsing each just returns the paths
# which you will need later
gsa <- parse_phyotozome(gsParam = gsp) 
```

### 2.3 Build the orthofinder database

GENESPACE relies on a database of blast-like diamond2 hits that are subsequently clustered into orthogroups by `orthofinder`. We offer two distinct methods for generation of these data, which are specified in `setup_genespace`

1. `orthofinderMethod = "fast"`: This approach runs `diamond2 --fast` on all unique combinations of genomes (not all pairwise combinations) and then stops after orthogroup inference via `orthofinder -og`. While this approach is roughly 10x faster than the default method, it is an approximation and should be used ONLY IN THE FOLLOWING CASES:
    - the genomes that are being analyzed (excluding the outgroup) are very closely related (e.g. within species pangenomes).
    - the user does not care much about the presence of highly-diverged gene arrays or other forms of gene duplications (e.g. WGD).
    - the user only wants to visually explore the patterns of synteny 
2. `orthofinderMethod = "default"`: This simply runs `orthofinder -f` which, by default runs `diamond2 --sensitive` and uses gene trees for dissection of orthogroups. While slower, this method will be far better at detecting ancient whole-genome duplications or orthology among diverged genomes. 

To build the orthofinder database from R:

```
blastDB <- build_OFDb(
  gsParam = gsp, # from setup_genespace
  gsAnnot = gsa) # from parse_phytozome (or whatever)
```

### 2.4 Set parameters for synteny constraint

GENESPACE relies on high-confidence, but broadly defined syntenic regions to vizualize data and build a pangenome. This step uses a priori defined ploidy and user-specified synteny parameters to split diamond2 hits into syntenic and non-syntenic regions. 

First, the user must specify the synteny parameters for the run. Here is an example, but there are other parameters that can be specified for polyploids etc. 

```
synp <- set_syntenyParams(
  gsParam = gsp, # from setup_genespace
  gsAnnot = gsa, # from parse_phytozome (or whatever)
  nGaps = 5, # passed to MCScanX -m parameter
  blkSize = 5, # passed to MCScanX -s parameter
  nSecondHits = 0, # number of hits beyond what is expected by ploidy
  synBuff = 100, # gene-rank order distance around MCScanX-defined 'anchors'
  onlyOg = TRUE) # only consider hits between members of the same orthogroup?
```

The default specification is very robust when using moderately evenly phylogenetically distributed genomes. However, this step is separate from overall genespace parameterization because, in some cases deciding how GENESPACE looks for synteny must be more carefully specified. For example, if running 10 grass genomes and pineapple, the user probably expects there to be far smaller and less contiguous blocks of synteny between the grasses and pineapple than within the grasses. As such, the user could edit the synp object to increase the number of gaps in searches involving pineapple as follows:

```
which2change <- with(synp, 
  genome1 != genome2 & (genome1 == "Pineapple" | genome2 == "Pineapple"))
synp$nGaps[which2change] <- 25
```

### 2.5 split hits into syntenic and non-syntenic regions

Once the synteny parameters are specified, they are passed to `extract_synteny`, which uses the rules to decide synteny with a simple command:

```
synh <- extract_synteny(
  synParam = synp, # from set_syntenyParams
  gsParam = gsp, # from setup_genespace
  gsAnnot = gsa) # from parse_phytozome (or whatever)
```

This automatically produces a low-resolution pdf with the positions of syntenic regions for each genome comparison, stored in /results/syntenicRegionHeatmap.pdf. We recommend that the user looks at this plot to ensure that the number and positions of syntenic regions looks right. If not, re-parameterize and re-run. 

### 2.6 make high-confidence syntenic blocks

The key to making high-confidence syntenic block coordinates is five-fold:

1. Restrict the search to broad syntenic regions 
2. Use only high-confidence hits; GENESPACE only 'sees' reciprocal-best-hits and hits where both target and query are in the same orthogroup. 
3. Use within-region gene-rank order as the x-y positions for block construction (reranked for only high confidence hits) 
4. Use MCScanX to cull hits to collinear regions (but not for block assignment)
5. Use a density-based block assigner, like `dbscan`. 

Genespace accomplishes this via `pull_synOgs`. 

```
synOgFiles <- pull_synOgs(
  synParam = synp, # from set_syntenyParams
  gsParam = gsp, # from setup_genespace
  gsAnnot = gsa) # from parse_phytozome (or whatever)
```

This produces all the files needed for GENESPACE plotting routines (see below), including the blockCoordinates. 

### 2.7 Construct the syntenic pangenome database

The final step of GENESPACE is to decode the syntenic hits, block coordinates and inferred syntenic positions of missing genes into a pangenome. This is accomplished by `build_pgDb`. 

```
pgFiles <- build_pgDb(
  synParam = synp, # from set_syntenyParams
  gsParam = gsp, # from setup_genespace
  gsAnnot = gsa, # from parse_phytozome (or whatever)
  synOgFiles = synOgFiles) # from pull_synOgs
```

This produces a pangenome annotation and saves it to a text file. If working in a fairly small space, the entire pangenome can be read into memory with:

```
pg <- read_pangenome(
  gsParam = gsp, 
  wholeThing = TRUE,
  referenceGenome = gsp$genomeIDs[1])
```

This converts a 'long' formatted pangenome text file into a more human-readable format with a column for each genome and a row for each unique position in the pangenome, anchored against physical positions of genes in the specified reference genome. 

For most cases, we expect the user to seek specific regions. This can be accomplished through several methods of querying. For example, to pull all elements of the pangenome in the first 10Mb of Chr01 for the first genome:

```
pg <- read_pangenome(
  gsParam = gsp, 
  referenceGenome = gsp$genomeIDs[1],
  chr = "Chr01",
  regionStart = 0, 
  regionEnd = 10e6)
```

# 3. Plotting

### 3.1 Dotplots

### 3.2 Riparian plots

### 3.3 Pangenome PAV plots



# 4. Legal

GENESPACE R Package (GENESPACE) Copyright (c) 2021,HudsonAlpha Institute for Biotechnology. All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.

NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative  works, and perform publicly and display publicly, and to permit others to do so.


*** License Agreement ***

GENESPACE R Package (GENESPACE) Copyright (c) 2021, HudsonAlpha Institute for Biotechnology. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the HudsonAlpha Institute for Biotechnology nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.
