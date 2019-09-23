---
title: GENESPACE readme
author: John T. Lovell
date: 19-September 2019
---

## Introduction

As chromosome-scale whole genome assemblies and annotations become more common across many taxonomic groups, it is crucial to be able to make evolutionary inference regarding sequences that are derived from common ancestors. There are two major lines of evidence to test whether genes are related: 1) sequence similarity (i.e. homology) and 2) similarity of physically proximate genes (i.e. synteny). If genes from two species are homologous, it is possible they arose from duplications (i.e. paralogs), are members of conserved similar gene familes but are unrelated, and arose from a speciation event and share a common ancestor (i.e. orthologs). The same is true for synteny, which can be related to ancient duplications or retained gene order since speciation. However, our confidence in the type of evolutionary relationship between sequences can be increased by testing for the presence of *syntenic* orthologs and paralogs. 

To date, a number of programs can infer synteny and homology, and a number of studies have combined the two principals. The R package, 'GENESPACE' does this explicitly by first inferring homologous gene networks ('orthogroups') through the `orthofinder` program, then parses orthogroups into syntenic blocks via `MCScanX`. Finally, blast hits are constrained to syntenic regions and classified as orthologous or paralogous through `orthofinder`. The three output files are 1) a bedmap with syntenic block coordinates, 2) an annotated blast file with all syntenic homologs, and 3) a bedmap with the expected syntenic blocks of all genes.

## What can GENESPACE do for me?

GENESPACE is an R package designed to improve the speed, sensitivity, utility and visualization of comparative genomics in the post-genomic world. It is best run interactively in R, but also can be called as a pipeline from the command line. 

The primary functionality of GENESPACE is **to develop high-confidence mappings between orthologous genes from multiple related genomes**. Depending on the system and parameters, GENESPACE can also effectively map paralogs among genomes that have undergone both recent and ancient whole genome duplication (WGD) events (e.g. how much synteny is retained, age of WGD, pattern of gene retention, etc.). GENESPACE also can produce highly-customizable publication-quality plots. It outputs datasets that can be used to infer a number of downstream attributes, including tandem arrays, pseudogenization, ancestral state reconstructions and selection / neutrality staatistics. 

##  Data Requirements

GENESPACE is written for default parameters to work directly off the JGI phytozome geneome annotation format (https://phytozome.jgi.doe.gov/pz/portal.html); however, nearly all formatting can be accomodated. 

**NOTE** Depending on the RAM available, GENESPACE can run up to 12 diploid genomes (6 tetraploid genomes, etc.). However, we recommend doing smaller runs (<= 8 genomes) unless synteny is at least marginally conserved across all species. 

To run GENESPACE, you need:

1. complete gene annotations from >= 2 genomes. 
2. gff3-formatted gene annotations 
3. fasta files containing peptide sequence for the primary transcript of each gene.
4. the gff3 and fasta files need to have gene identifiers that can be used to link the two.
5. there must be some synteny between genomes. GENESPACE cannot be used to find related genes between very diverged species (i.e. plants and animals), nor should it be used to look at very ancient whole genome duplications. 

**NOTE** Genomes must have decent contiguity. The required level of contiguity depends on the minimum desired syntenic block size. But, the larger the syntenic blocks, the better the results. So, if using genomes with small or un-ordered contigs, take the results lightly. 

## System Requirements

1. MacOSX or LINUX 
2. R v3.5.2 or later
3. OrthoFinder version 2.3.3 or later
4. MCScanX downloaded from http://chibba.pgml.uga.edu/mcscan2/


## Running GENESPACE interactively

With the initial v0.5 release, the GENESPACE pipeline must be run interactively. Later releases will have a commaand line one-liner to call the whole thing. 

### Data formatting 

GENESPACE requires an orthofinder-formatted blast search result, where the gff3 'attribute' columns match the gene IDs in the fasta headers. We can achieve these results via the sata preparatation pipeline:

1. `convert_genomes`: Format / simplify fasta head and match the gene ID in gff3. 
2. `run_orthofinder`: Wrapper to run orthofinder and aggregate results
3. `remake_ofInput`: Subset blast results (i.e. to top hits by gene)
4. `rerun_orthofinder`: Re-run orthofinder on subsetted blast results.

### Forming and finalizing syntenic blocks

Determining sequence similarity is fairly trival. Programs like orthofinder can parse such blast results. 

5. `pipe_mcscanx`: an R wrapper for the MCScanX program. Clusters hits into collinear blocks and drops non-collinear blast hits. 
6. `clean_blocks`: uses fixed-radius 2d density clustering to retain hits in high-quality blocks. 
7. `complete_graph`: fills out incomplete complete networks (e.g. for genes A,B,C: [A-B,B-C] --> [A-B,B-C,C-A])
8. `extend_blocks`: searches for any blast hits within a fixed radius of previously defined syntenic blocks. Crucial for under-retained duplicate regions, as these may get broken up into small pieces by MCScanX and clean_blocks. Should be followed by a final pipe_mcscanx run. 

### Classifying hits as orthologous, paralogous or other

9. `assign_synHomologs`: searches for all syntenic blast hits, feeds these to orthofinder, builds gene trees and returns an annotated blast file with a column specifying the hit as ortholog, paralog, or syntenic homolog.

### An example with mammals. 







