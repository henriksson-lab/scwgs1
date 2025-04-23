bascetRoot <- getwd()

if(TRUE){
  
  setwd("~/github/zorn")
  
  source("R/job_general.R")
  source("R/job_local.R")
  source("R/job_slurm.R")
  source("R/bascet_file.R")
  source("R/zorn.R")
  source("R/shell.R")
  source("R/zorn_aggr.R")
  source("R/count_kmer.R")
  source("R/countsketch.R")
  source("R/refgenome.R")
  source("R/kraken.R")
  source("R/container.R")
  source("R/ext_tools.R")
  
} else {
  
  library(Zorn)
  
}


bascet_instance.default <- getBascetSingularityImage(store_at="~/mystore/")
#bascet_runner.default <- SlurmRunner(account="hpc2n2025-074", ncpu="10")
bascet_runner <- LocalRunner(direct = TRUE, show_script=FALSE)

setwd(bascetRoot)





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



