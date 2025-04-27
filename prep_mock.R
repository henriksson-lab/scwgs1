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
bascet_runner.default <- SlurmRunner(account="hpc2n2025-074", ncpu="10")
#####bascet_runner <- LocalRunner(direct = TRUE, show_script=TRUE)


setwd(bascetRoot)

###################################################
################################################### debarcode
###################################################

rawmeta_dir <- readLines("rawdata.txt")
rawmeta <- DetectRawFileMeta(rawmeta_dir)

barcode_error <- NULL

if(file.exists("barcode_error.txt")) {
  barcode_error <- readLines("barcode_error.txt")
}



system("echo START BascetGetRaw >> time.txt; echo `date +%s` >> time.txt")
BascetGetRaw(
    bascetRoot,
    rawmeta
)
system("echo END BascetGetRaw >> time.txt; echo `date +%s` >> time.txt")

###################################################
################################################### shardify
###################################################

### Decide cells to include
h <- ReadHistogram(bascetRoot,"debarcoded")
#includeCells <- h$cellid[h$count>10]       ########### 10 for miseq
includeCells <- h$cellid[h$count>100]       ########### 10 for miseq
length(includeCells)

### Shardify i.e. divide into multiple sets of files for parallel processing.
# this command will spawn many processes, as it does random I/O on the input files
system("echo START BascetShardify >> time.txt; echo `date +%s` >> time.txt")
BascetShardify(
  bascetRoot,
  includeCells = includeCells,
  #num_output_shards = 5, ############### 10,
  #num_output_shards = 10,
  num_output_shards = 20,
  runner=SlurmRunner(bascet_runner.default, ncpu="4")  #not much CPU needed. increased for memory demands
)
system("echo END BascetShardify >> time.txt; echo `date +%s` >> time.txt")



### Get reads in fastq format for BWA
system("echo START asfq >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  inputName="filtered",
  outputName="asfq",
  out_format="R1.fq.gz",
  runner=SlurmRunner(bascet_runner.default, ncpu="4")  #in and out processes
)
system("echo END asfq >> time.txt; echo `date +%s` >> time.txt")





### Get reads in fastq format for BWA
system("echo START fastp >> time.txt; echo `date +%s` >> time.txt")
BascetRunFASTP(
  bascetRoot,
  numLocalThreads=10,
  inputName="asfq",
  outputName="fastp"
)
system("echo END fastp >> time.txt; echo `date +%s` >> time.txt")



system("echo START totirp >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  "fastp",
  "new_filtered",
  out_format="tirp.gz"
)
system("echo END totirp >> time.txt; echo `date +%s` >> time.txt")




################################################################################
################## Alignment ###################################################
################################################################################



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
  runner=SlurmRunner(bascet_runner.default, ncpu="2") ################################# todo save in temporary pos, then move
)
system("echo END BascetCountChrom >> time.txt; echo `date +%s` >> time.txt")







################################################################################
################## Preprocessing with KRAKEN ###################################
################################################################################


### Run Kraken on each cell
system("echo START BascetRunKraken >> time.txt; echo `date +%s` >> time.txt")
BascetRunKraken(
  bascetRoot,
  useKrakenDB="/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8",
  numLocalThreads=10,
  inputName="fastp",
  runner=SlurmRunner(bascet_runner.default, ncpu="14")  #needs a lot of memory!
)
system("echo END BascetRunKraken >> time.txt; echo `date +%s` >> time.txt")

system("echo START BascetMakeKrakenCountMatrix >> time.txt; echo `date +%s` >> time.txt")
BascetMakeKrakenCountMatrix(
  bascetRoot,
  numLocalThreads=10
)
system("echo END BascetMakeKrakenCountMatrix >> time.txt; echo `date +%s` >> time.txt")



################################################################################
################## count sketching #############################################
################################################################################


system("echo START BascetComputeCountSketch >> time.txt; echo `date +%s` >> time.txt") 
BascetComputeCountSketch(
  bascetRoot,
  inputName="new_filtered"
)
system("echo END BascetComputeCountSketch >> time.txt; echo `date +%s` >> time.txt")



system("echo START BascetGatherCountSketch >> time.txt; echo `date +%s` >> time.txt")
BascetGatherCountSketch(
  bascetRoot
)
system("echo END BascetGatherCountSketch >> time.txt; echo `date +%s` >> time.txt")




################################################################################
################## Informative KMER ############################################
################################################################################



### Compute minhashes for each cell
system("echo START BascetComputeMinhash >> time.txt; echo `date +%s` >> time.txt")
BascetComputeMinhash(
  bascetRoot,
  inputName = "new_filtered",
)
system("echo END BascetComputeMinhash >> time.txt; echo `date +%s` >> time.txt")


### Gather minhashes into a single histogram
system("echo START BascetMakeMinhashHistogram >> time.txt; echo `date +%s` >> time.txt")
BascetMakeMinhashHistogram(
  bascetRoot
)
system("echo END BascetMakeMinhashHistogram >> time.txt; echo `date +%s` >> time.txt")


### Pick KMERs
#kmer_hist <- BascetReadMinhashHistogram(bascetRoot)
#useKMERs <- kmer_hist$kmer[kmer_hist$cnt>5]

## Build count table by looking up selected KMERs in per-cell KMER databases
if(file.exists("use_kmers.txt")){
 useKMERs <- readLines("use_kmers.txt")

 system("echo START BascetQueryFq >> time.txt; echo `date +%s` >> time.txt")
 BascetQueryFq(
   bascetRoot,
   useKMERs=useKMERs
 )
 system("echo END BascetQueryFq >> time.txt; echo `date +%s` >> time.txt")
}






################################################################################
################## De novo Assembly ############################################
################################################################################



### Assemble all genomes
system("echo START skesa >> time.txt; echo `date +%s` >> time.txt")
BascetMapCell(
  bascetRoot,
  withfunction = "_skesa",
  inputName = "new_filtered",
  outputName = "skesa"
)
system("echo END skesa >> time.txt; echo `date +%s` >> time.txt")




### quast QC
system("echo START quast >> time.txt; echo `date +%s` >> time.txt")
BascetMapCell(
  bascetRoot,
  withfunction = "_quast",
  inputName = "skesa",
  outputName = "quast"
)
system("echo END quast >> time.txt; echo `date +%s` >> time.txt")




### abricate QC
system("echo START abricate >> time.txt; echo `date +%s` >> time.txt")
BascetMapCell(
  bascetRoot,
  withfunction = "_abricate",
  inputName = "skesa",
  outputName = "abricate",
  args = list(DATABASE_DIR="ncbi")
)
system("echo END abricate >> time.txt; echo `date +%s` >> time.txt")


### FastQC
system("echo START BascetMapCellFASTQC >> time.txt; echo `date +%s` >> time.txt")
BascetMapCellFASTQC(
  bascetRoot,
  inputName = "new_filtered"  # todo make mapcell work on fastq as input too
)
system("echo END BascetMapCellFASTQC >> time.txt; echo `date +%s` >> time.txt")



