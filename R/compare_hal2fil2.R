tmp.dir = "/Users/jlovell/Documents/comparative_genomics/hal2_vs_fil2/tmp"
results.dir = "/Users/jlovell/Documents/comparative_genomics/hal2_vs_fil2/results"
blast.dir = "/Users/jlovell/Documents/comparative_genomics/hal2_vs_fil2/blast"
input.dir = "/Users/jlovell/Documents/comparative_genomics/hal2_vs_fil2/genomes"
mcscan.dir = "/Users/jlovell/Documents/comparative_genomics/hal2_vs_fil2/mcscanx"

peptide.dir = file.path(input.dir,"peptide")
cds.dir = file.path(input.dir,"cds")
assembly.dir = file.path(input.dir,"assembly")
gff.dir = file.path(input.dir,"gff")

### Metadata must include:
genomeIDs = c("PhalliiHAL","Phallii","Sviridis","Sbicolor")
ploidy = c(2,2,2,2)
abbrevs = c("Ph","Pa","Sv","Sb")

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

## 2.3 -- Make preliminary blocks via MCScanX and some post-processing
mcscan.raw = run_MCScanX(
  blast.results = blast.results,
  abbrevs = abbrevs,
  mcscanx.input.dir = mcscan.dir,
  MCScanX.params = "-a -s 5 -m 10")
tmp = make_blocks(mcscan.raw)
tmp = add2_blocks(blast.results = blast.results,
                    map = tmp$map,
                    blk = tmp$block,
                    buffer = 1, n.iter = 10)
tmp = find_newBlocks(blk = tmp$block, map = tmp$map,
                      blast.results = blast.results, min.unique.genes = 5)
tmp = add2_blocks(blast.results = blast.results,
                     map = tmp$map,
                     blk = tmp$block,
                     buffer = 1, n.iter = 10)

# 2.4 -- Re-run MCScanX with just the blast hits that passed 2.3
tmp.map = tmp$map[,c("id1","id2","genome1","genome2")]
blast.cull = merge(blast.results, tmp.map, by = c("id1","id2","genome1","genome2"))
mcscan.cull = run_MCScanX(
  blast.results = blast.cull,
  abbrevs = abbrevs,
  mcscanx.input.dir = mcscan.dir,
  MCScanX.params = "-a -s 10 -m 15")
final = make_blocks(mcscan.cull)

for(i in 1:3){
  final = merge_blocks(
    map = final$map,
    blk = final$block,
    buffer = .5,
    verbose = T)
}
final = add2_blocks(blast.results = blast.cull,
                    map = final$map,
                    blk = final$block,
                    buffer = 0.5, n.iter = 2)

combined = merge_blockBreakPoints(map = final$map,
                                  genomeIDs = genomeIDs,
                                  checkOvl.bp = 2e5,
                                  checkOvl.rank = 4)

tosplit = finalize_mappings(blk = combined,
                  map = final$map,
                  gff.dir = gff.dir,
                  genomeIDs = genomeIDs)

