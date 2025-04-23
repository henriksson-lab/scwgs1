
### Perform alignment
system("echo START BascetAlignToReference >> time.txt; echo `date +%s` >> time.txt")
BascetAlignToReference(
  bascetRoot,
  useReference="/home/m/mahogny/mystore/atrandi/bwa_ref/ref10/all.fa",
  numLocalThreads=10
)
system("echo END BascetAlignToReference >> time.txt; echo `date +%s` >> time.txt")



### Generate fragments BED file suited for quantifying reads/chromosome using Signac later -- this is a wrapper for mapshard
system("echo START BascetBam2Fragments >> time.txt; echo `date +%s` >> time.txt")
BascetBam2Fragments(
  bascetRoot,
  runner=SlurmRunner(bascet_runner, ncpu="2")
)
system("echo END BascetBam2Fragments >> time.txt; echo `date +%s` >> time.txt")


### Count reads per chromosome
system("echo START BascetCountChrom >> time.txt; echo `date +%s` >> time.txt")
BascetCountChrom(
  bascetRoot,
  runner=SlurmRunner(bascet_runner, ncpu="2")
)
system("echo END BascetCountChrom >> time.txt; echo `date +%s` >> time.txt")



