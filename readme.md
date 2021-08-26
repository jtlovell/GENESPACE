---
title: "A guide to GENESPACE v0.7"
author: "JT Lovell"
date: "25-Aug 2021"
---

# 1. Overview

GENESPACE is a comparative genomics framework implemented in the R environment for statistical computing. The premise is that, when analyzing high-quality genome assemblies and annotations, we can improve the confidence of evolutionary inference by combining two sources of evidence for common ancestry: synteny/collinearity of gene order on chromosomes and coding sequence similarity (homology). In addition to providing a second line of evidence beyond sequence similarity, synteny offers a mechanism by which gene order can be predicted across genomes, even when presence-absence variation means that not all genomes are represented in every orthology network ('orthogroup'). Further, since all inference is constrained within syntenic regions, genespace is agnostic to ploidy, segmental duplicates, inversions and other chromosomal complexities. 

GENESPACE outputs a synteny-constrained and -anchored orthogroup pan-genome annotation among multiple genomes. This simple text file allows for extraction and exploration of regional (e.g. a QTL interval) gene-level variation, a necessary step to integrate comparative and quantitatve genomic goals. 

# 2. Quick start

Assuming `orthofinder`, `MCScanX`, `GENESPACE` and dependencies are installed, and genomes being used are all from NCBI, here is the full pipeline to build syntenic pan-genome annotations and make plots. 

## 2.1 Parameterize the run

Here, we are going to use human, chimpanzee and mouse genomes from NCBI. Each genome is a haploid representation of the diploid genotypes, so ploidy = 1. 

```
gpar <- init_genespace(
  genomeIDs = c("human","chimp","mouse"),
  speciesIDs = c("Homo_sapiens","Pan_troglodytes","Mus_musculus"),
  versionIDs = c("GRCh38.p13", "Clint_PTRv2", "GRCm39"),
  ploidy = c(1,1,1),
  orthofinderMethod = "default",
  wd = "/mammalGenespace/run1",
  nCores = 4, 
  overwrite = T,
  gffString = "gff", 
  pepString = "translated_cds",
  path2orthofinder = "NA",
  path2mcscanx = "/programs/MCScanX",
  rawGenomeDir = "/mammalGenespace/rawRepo"))
```

## 2.2 Parse the annotations

Since these are from NCBI, there is a simple parser that formats the annotations for use in GENESPACE. 

```
parse_ncbi(gsParam = gpar)
```

## 2.3 Set the desired synteny parameters.

 There are a lot of parameters, see `?set_syntenyParams`, but the defaults should be fine for most situations.

```
gpar <- set_syntenyParams(
  gsParam = gpar)
```

To see the full parameterization, print the synteny parameters as `print(gpar$params$synteny)`. This can be edited manually if needed to fine-tune the synteny search to diverse patterns of genome evolution

## 2.4 Get the command to run orthofinder

If orthofinder is installed and in the path, this command will run orthofinder from R, but this is not necessary, or even preferable. Instead, it will print the command and you can run that in an environment with orthofinder in the path. 

```
gpar <- run_orthofinder(gsParam = gpar, overwrite = F)
```

All results from the orthofinder run are stored externally to R, but downstream GENESPACE functions know where to look for them.

## 2.5 constrain blast results to synteny

The primary engine of GENESPACE is `synteny`, which parses an orthofinder run into syntenic orthogroups. This operates on the syntenyParameters generated from `set_syntenyParam`. 

```
blks <- synteny(gsParam = gpar)
```

The synteny-constrained hits and corresponding dotplots are stored by default in the /results directory. However, this can be changed by alterning the path to results in the gsParams. For example `gpar$paths$results <- "/results2"; dir.create(gpar$paths$results)`. For convenience, the syntenic block coordinates are returned, but not needed for downstream analyses. 

## 2.6 build a pan-genome annotation

Here, the synteny constrained hits are used to predict their positions against a chosen reference genome. Then the orthology networks are decoded into a tabular format. 

```
pg <- pangenome(gsParam = gpar, refGenome = "human")
```

## 2.7 Plotting

In addition to the dotplots made within `synteny`, GENESPACE can generate two other types of plots. The first is a riparian, which tracks synteny across multiple genomes, color coding and ordering chromosomes by synteny to a reference genomes. For example:

```
ripd  <- plot_riparian(
  gsParam = gpar, 
  plotRegions = T,
  useOrder = T, 
  genomeIDs = c("human","chimp","mouse"), 
  refGenome = "chimp")
```

The degree of presence absence variation in the global any synteny-constrained orthogroups can be plotted as:

```
plot_pav(
  gsParam = gPar, 
  genomeIDs = c("human","chimp","mouse"), 
  maxComb2plot = 50, 
  excudePrivate = F) 
```


# 3. Setup

### 3.1 Installing Dependencies

The GENESPACE R package is designed only for use MacOSX and Linux. Testing on OSX 10.16 (Big Sur) and __ (). The R package requires that `MCScanX` is [installed](https://github.com/wyp1125/MCScanX) as well as the following R dependencies:

- R/data.table (v1.13.6 or later)
- R/igraph (v1.2.6 or later)
- R/Biostrings (v2.58.0 or later)
- R/dbscan (v1.1-5 or later)

GENESPACE cannot be run without results from `orthofinder`. This can be installed on a separate server or locally. The simplest method is with conda: `conda install -c bioconda orthofinder`, but other methods are available. See [orthofinder documentation](https://github.com/davidemms/OrthoFinder). 


### 3.2 Install GENESPACE

To install GENESPACE, enter R, either by typing `R` into the command line or `open -na rstudio` for RStudio. If not installed, install the `devtools` package via `install.packages("devtools")`. This will permit one-line installation from the gitlab repo. 

Install genespace: `devtools::install_git("https://code.jgi.doe.gov/plant/genespace-r")`

This should also install the dependencies, `data.table`, `igraph`, and `dbscan`. For full functionality, also install the following:
`Biostrings` and `Rutils`. 

If devtools is not available or you'd prefer, you can install from the provided [tarball file](https://code.jgi.doe.gov/plant/genespace-r/-/archive/master/genespace-r-master.tar.gz), saved in a location where you have execute privileges. If running on the command line, prior to entering R, run: `R cmd install path/to/GENESPACE.tar.gz`. Alternatively, if running GENESPACE from RStudio, open Rstudio once the orthofinder conda environment is activated, then install GENESPACE by clicking on the 'install' icon under the 'packages' tab. Select 'install from package archive file' and navigate to the location of the saved GENESPACE.tar.gz file. 

### 3.3 Input data format

GENESPACE requires two inputs for each genome: primary transcript annotation peptide fasta files and gene feature coordinates. The sequence headers in the fasta must be able to parsed to match the attributes column of the gff3 gene annotation. 

If data comes from phytozome, you want the 'gene.gff3' and 'protein_primaryTranscriptOnly.fa' files. 

If data comes from NCBI, you to get the data from the 'FTP directory for RefSeq assembly' file tree and download the 'genomic.gff.gz' and 'translated_cds.faa.gz' files. 

If getting data from other sources, you may need to ensure that your gene IDs can be properly matched. It is not necessary to have only the primary transcript, as GENESPACE will take the longest uniquely named sequence. However, in many cases the longest transcript is not the primary, so care should be taken. 


# 4. Pipeline details

### 4.1 init_genespace parameters

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

### 4.2 Convert to GENESPACE file format

The basic gff and peptide header conversions are accomplished by `parse_annotations` or the more user-friendly `parse_ncbi` and `parse_phytozome` functions. 

This is how `parse_annotations` can be parameterized. It should allow for pretty much any format of fasta and gff3.
```
parse_annotations(
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
with(gsp, parse_annotations(
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


# 5. Legal

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
