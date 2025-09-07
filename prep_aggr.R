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


bascetInstance.default <- getBascetSingularityImage(storeAt="~/mystore/")
#bascet_runner.default <- SlurmRunner(account="hpc2n2025-074", ncpu="10")
bascet_runner <- LocalRunner(direct = TRUE, showScript=FALSE)

setwd(bascetRoot)


if(file.exists("quast.1.zip")) {
  system("echo START BascetAggregateQUAST >> time.txt; echo `date +%s` >> time.txt")
  print("quast")
  BascetCacheComputation(
    bascetRoot,
    "cache_quast",
    BascetAggregateQUAST(
      bascetRoot
    )
  )
  system("echo END BascetAggregateQUAST >> time.txt; echo `date +%s` >> time.txt")
}


if(file.exists("abricate.1.zip")) {
  system("echo START BascetAggregateAbricate >> time.txt; echo `date +%s` >> time.txt")
  print("abricate")
  BascetCacheComputation(
    bascetRoot,
    "cache_abricate",
    BascetAggregateAbricate(
      bascetRoot
    )
  )
  system("echo END BascetAggregateAbricate >> time.txt; echo `date +%s` >> time.txt")  
}





if(file.exists("fastqc.1.zip")) {
  system("echo START BascetAggregateFASTQC >> time.txt; echo `date +%s` >> time.txt")
  print("fastqc")
  BascetCacheComputation(
    bascetRoot,
    "cache_fastqc",
    BascetAggregateFASTQC(
      bascetRoot, verbose=TRUE
    )
  )
  system("echo END BascetAggregateFASTQC >> time.txt; echo `date +%s` >> time.txt")  
}  





