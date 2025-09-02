
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

bascetInstance.default <- getBascetSingularityImage(store_at="~/mystore/")
bascetRunner.default <- SlurmRunner(account="hpc2n2025-074", ncpu="10")
#####bascet_runner <- LocalRunner(direct = TRUE, show_script=TRUE)

setwd(bascetRoot)

###################################################
################################################### debarcode
###################################################


rawmeta_dir <- readLines("rawdata.txt")
rawmeta <- DetectRawFileMeta(rawmeta_dir)

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
  num_output_shards = 20,
  runner=SlurmRunner(bascet_runner.default, ncpu="4")  #not much CPU needed. increased for memory demands
)
system("echo END BascetShardify >> time.txt; echo `date +%s` >> time.txt")


### Get reads in fastq format for BWA
system("echo START asfq >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  "filtered",
  "asfq",
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



################################################################################
################## Alignment ###################################################
################################################################################

### Perform alignment
system("echo START BascetAlignToReference >> time.txt; echo `date +%s` >> time.txt")
BascetAlignToReference(
  bascetRoot,
  inputName="fastp",
  useReference="/home/m/mahogny/mystore/atrandi/bwa_ref/combo_hx/all.fa.gz",   #for saliva -- should filter out these reads
  numLocalThreads=10
)
system("echo END BascetAlignToReference >> time.txt; echo `date +%s` >> time.txt")




if(FALSE){

### Generate fragments BED file suited for quantifying reads/chromosome using Signac later -- this is a wrapper for mapshard
system("echo START BascetBam2Fragments >> time.txt; echo `date +%s` >> time.txt")
BascetBam2Fragments(
  bascetRoot,
  runner=SlurmRunner(bascet_runner.default, ncpu="2")
)
 
system("echo END BascetBam2Fragments >> time.txt; echo `date +%s` >> time.txt")

}

### Count reads per chromosome
system("echo START BascetCountChrom >> time.txt; echo `date +%s` >> time.txt")
BascetCountChrom(
  bascetRoot,
  min_matching=50,
  runner=SlurmRunner(bascet_runner.default, ncpu="10")
)
system("echo END BascetCountChrom >> time.txt; echo `date +%s` >> time.txt")


################################################################################
################## Host filtering ##############################################
################################################################################


system("echo START BascetFilterAlignment >> time.txt; echo `date +%s` >> time.txt")
BascetFilterAlignment(  
    bascetRoot, 
    numLocalThreads=1,
    inputName="unsorted_aligned", 
    outputName="nohost_aligned",
    keep_mapped=FALSE,
    runner=SlurmRunner(bascet_runner.default, ncpu="1") #foo
)
system("echo END BascetFilterAlignment >> time.txt; echo `date +%s` >> time.txt")




################################################################################
################## back to tirp ################################################   back to TIRP			if later commands could take BAM, this would not be needed  TODO
################################################################################


system("echo START totirp >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  inputName="nohost_aligned",
  outputName = "new_filtered",
  out_format="tirp.gz"
)
system("echo END totirp >> time.txt; echo `date +%s` >> time.txt")



######## and a new FASTQ file

system("echo START tofastq_filtered >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  inputName="nohost_aligned",
  outputName = "new_filtered_fa",
  out_format="R1.fq.gz"
)
system("echo END tofastq_filtered >> time.txt; echo `date +%s` >> time.txt")




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
################## Preprocessing with KRAKEN ###################################
################################################################################


### Run Kraken on each cell
system("echo START BascetRunKraken >> time.txt; echo `date +%s` >> time.txt")
BascetRunKraken(
  bascetRoot,
  #inputName = "fastp",
  inputName = "new_filtered_fa",
  useKrakenDB="/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8",
  numLocalThreads=10,
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
################## de novo kmer analysis #######################################
################################################################################




# /home/m/mahogny/mystore/atrandi/kraken_ref/standard-8


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


### Pick KMERs -- do in R
#kmer_hist <- BascetReadMinhashHistogram(bascetRoot)
#useKMERs <- kmer_hist$kmer[kmer_hist$cnt>5]


## Build count table by looking up selected KMERs in per-cell KMER databases
if(file.exists("use_kmers.txt")){

 useKMERs <- readLines("use_kmers.txt")

 system("echo START BascetQueryFq >> time.txt; echo `date +%s` >> time.txt")
 BascetQueryFq(
   bascetRoot,
   inputName = "new_filtered",
   useKMERs=useKMERs
 )
 system("echo END BascetQueryFq >> time.txt; echo `date +%s` >> time.txt")

}




################################################################################
################## Assembly ####################################################
################################################################################



### Assemble all genomes
system("echo START skesa >> time.txt; echo `date +%s` >> time.txt")
BascetMapCell(
  bascetRoot,
  withfunction = "_skesa",
  inputName = "new_filtered",
  outputName = "skesa",
  runner=SlurmRunner(bascet_runner.default, ncpu="16") #needs a lot of cpu. oddly at the end??
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













# 
# 
# how many % humna, % xantham, other?
# barplots
#
# b-an01 [~/mystore/atrandi/v4_wgs_saliva1]$ ls *nohost_aligned* -lh
# -rw-rw----+ 1 mahogny folk  1,0G apr 18 20:02 nohost_aligned.1.bam
# -rw-rw----+ 1 mahogny folk  1,5G apr 18 20:02 nohost_aligned.10.bam
# -rw-rw----+ 1 mahogny folk  1,4G apr 18 20:02 nohost_aligned.11.bam
# -rw-rw----+ 1 mahogny folk  1,5G apr 18 20:02 nohost_aligned.12.bam
# -rw-rw----+ 1 mahogny folk  880M apr 18 20:01 nohost_aligned.13.bam
# -rw-rw----+ 1 mahogny folk  1,2G apr 18 20:02 nohost_aligned.14.bam
# -rw-rw----+ 1 mahogny folk  1,1G apr 18 20:01 nohost_aligned.15.bam
# -rw-rw----+ 1 mahogny folk  1,5G apr 18 20:03 nohost_aligned.16.bam
# -rw-rw----+ 1 mahogny folk  1,3G apr 18 20:03 nohost_aligned.17.bam
# -rw-rw----+ 1 mahogny folk  1,2G apr 18 20:02 nohost_aligned.18.bam
# -rw-rw----+ 1 mahogny folk  1,5G apr 18 20:03 nohost_aligned.19.bam
# -rw-rw----+ 1 mahogny folk  1,1G apr 18 20:02 nohost_aligned.2.bam
# -rw-rw----+ 1 mahogny folk  1,2G apr 18 20:02 nohost_aligned.20.bam
# -rw-rw----+ 1 mahogny folk  1,3G apr 18 20:02 nohost_aligned.3.bam
# -rw-rw----+ 1 mahogny folk  1,3G apr 18 20:02 nohost_aligned.4.bam
# -rw-rw----+ 1 mahogny folk  1,6G apr 18 20:03 nohost_aligned.5.bam
# -rw-rw----+ 1 mahogny folk  1,7G apr 18 20:03 nohost_aligned.6.bam
# -rw-rw----+ 1 mahogny folk  1,2G apr 18 20:02 nohost_aligned.7.bam
# -rw-rw----+ 1 mahogny folk 1012M apr 18 20:01 nohost_aligned.8.bam
# -rw-rw----+ 1 mahogny folk  1,5G apr 18 20:02 nohost_aligned.9.bam
# b-an01 [~/mystore/atrandi/v4_wgs_saliva1]$ ls unsorted_aligned.* -lh
# -rw-rw----+ 1 mahogny folk 2,1G apr  8 15:51 unsorted_aligned.1.bam
# -rw-rw----+ 1 mahogny folk 1,9G apr  8 15:49 unsorted_aligned.10.bam
# -rw-rw----+ 1 mahogny folk 2,2G apr  8 15:51 unsorted_aligned.11.bam
# -rw-rw----+ 1 mahogny folk 2,4G apr  8 15:54 unsorted_aligned.12.bam
# -rw-rw----+ 1 mahogny folk 3,0G apr  8 16:00 unsorted_aligned.13.bam
# -rw-rw----+ 1 mahogny folk 2,9G apr  8 16:00 unsorted_aligned.14.bam
# -rw-rw----+ 1 mahogny folk 2,2G apr  8 15:51 unsorted_aligned.15.bam
# -rw-rw----+ 1 mahogny folk 1,9G apr  8 15:50 unsorted_aligned.16.bam
# -rw-rw----+ 1 mahogny folk 2,7G apr  8 15:58 unsorted_aligned.17.bam
# -rw-rw----+ 1 mahogny folk 2,6G apr  8 15:56 unsorted_aligned.18.bam
# -rw-rw----+ 1 mahogny folk 2,5G apr  8 15:54 unsorted_aligned.19.bam
# -rw-rw----+ 1 mahogny folk 2,7G apr  8 15:53 unsorted_aligned.2.bam
# -rw-rw----+ 1 mahogny folk 1,6G apr  8 15:48 unsorted_aligned.20.bam
# -rw-rw----+ 1 mahogny folk 2,3G apr  8 15:51 unsorted_aligned.3.bam
# -rw-rw----+ 1 mahogny folk 2,0G apr  8 15:50 unsorted_aligned.4.bam
# -rw-rw----+ 1 mahogny folk 2,8G apr  8 15:58 unsorted_aligned.5.bam
# -rw-rw----+ 1 mahogny folk 2,5G apr  8 15:56 unsorted_aligned.6.bam
# -rw-rw----+ 1 mahogny folk 2,2G apr  8 15:52 unsorted_aligned.7.bam
# -rw-rw----+ 1 mahogny folk 2,8G apr  8 15:58 unsorted_aligned.8.bam
# -rw-rw----+ 1 mahogny folk 2,3G apr  8 15:53 unsorted_aligned.9.bam



