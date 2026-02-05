### preprocessing of data done in separate file, prep_*.R
if(FALSE) {
  #remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
  BiocManager::install("scDblFinder")
}

##################### TODO read (base) mahogny@beagle:/husky/henriksson/atrandi/mock2$ samtools idxstats aligned.1.bam > alstats.tsv

if(FALSE){
  source("/home/mahogny/github/zorn/R/job_general.R")
  source("/home/mahogny/github/zorn/R/job_local.R")
  source("/home/mahogny/github/zorn/R/job_slurm.R")
  source("/home/mahogny/github/zorn/R/bascet_file.R")
  source("/home/mahogny/github/zorn/R/zorn.R")
  source("/home/mahogny/github/zorn/R/shell.R")
  source("/home/mahogny/github/zorn/R/zorn_aggr.R")
  source("/home/mahogny/github/zorn/R/aggr_functions.R")
  source("/home/mahogny/github/zorn/R/count_kmer.R")
  source("/home/mahogny/github/zorn/R/countsketch.R")
  source("/home/mahogny/github/zorn/R/refgenome.R")
  source("/home/mahogny/github/zorn/R/kraken.R")
  source("/home/mahogny/github/zorn/R/container.R")
  source("/home/mahogny/github/zorn/R/ext_tools.R")
  
} else {
  library(Zorn)
}


list_color_mock <- c(khroma::color("muted")(n = 9),"#DDDDDD")

library(Seurat)
library(Signac)
library(sqldf)
library(stringr)
library(ggplot2)
library(Matrix)
library(GenomicRanges)

library(future)
plan("multicore", workers = 10)
#plan("multiprocess", workers = 20)

## Mock community: which chromosome is which strain
#mapSeq2strain <- read.csv("~/github/scwgs/map_seq2strain.csv")
mapSeq2strain <- read.csv("/husky/henriksson/atrandi/bwa_ref/atcc/map_seq2strain.csv")
mapSeq2strain$id <- stringr::str_replace_all(mapSeq2strain$id,"_","-")
strain_genomesize <- sqldf::sqldf("select 0 as len, strain from mapSeq2strain group by strain")

bascet_runner <- LocalRunner(direct = TRUE)

list_datasets <- c(
  "mock2"
)



#dataset_name <- "mock2"
dataset_name <- "v6_251128_jyoti_mock_bulk"
dataset_name <- "v6_251205_saliva2_mda"
dataset_name <- "v6_260116_mock_pta"
dataset_name <- "v7_260116_mock_pta"
#bascetRoot = "/husky/henriksson/atrandi/v2_wgs_miseq2/"  #for development
#bascetRoot <- "/big/henriksson/atrandi_mock_prev" #17gb
bascetRoot <- "/husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk"
bascetRoot <- "/husky/henriksson/atrandi/v6_251205_saliva2_mda"
bascetRoot <- "/husky/henriksson/atrandi/v6_260116_mock_pta"
bascetRoot <- "/husky/henriksson/atrandi/v6_simulated4"
bascetRoot <- "/husky/henriksson/atrandi/v7_260116_mock_pta"
# or? /big/henriksson/atrandi_mock #271 gb

#bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)
plotDir <- file.path("~/github/scwgs/plots/",dataset_name)
if(!file.exists(plotDir)){
  dir.create(plotDir, recursive = TRUE)
  print("made plot dir")
}

plotDirAll <- file.path("~/github/scwgs/plots/all")
if(!file.exists(plotDirAll)){
  dir.create(plotDirAll, recursive = TRUE)
  print("made plot dir all")
}


#100 for miseq??
min_nCount_RNA <- 1000

setwd("/home/mahogny/github/zorn") #to put SQL in the right place

################################################################################
################################################################################ mock2 alignment stats
################################################################################

# reference sequence name, sequence length, # mapped reads and # unmapped reads.
alstats <- read.table("/husky/henriksson/atrandi/mock2/alstats.tsv")
df <- data.frame(
  name=alstats$V1,
  cnt=rowSums(alstats[,3:4])
)
sum(df$cnt[str_starts(df$name,"NZ")])/sum(df$cnt)

23121616/sum(alstats[,4])  #98% do not map



##################### TODO read (base) mahogny@beagle:/husky/henriksson/atrandi/mock2$ samtools idxstats aligned.1.bam > alstats.tsv


################################################################################
##################### new duplication analysis #################################
################################################################################

cnt_with_dup <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)
cnt_no_dup <- ReadBascetCountMatrix(bascetRoot,"chromcount_dedup", verbose=FALSE)

df <- merge(
  data.frame(
    cell = rownames(cnt_with_dup$X),
    cnt_with_dup = rowSums(cnt_with_dup$X)
  ),
  data.frame(
    cell = rownames(cnt_no_dup$X),
    cnt_no_dup = rowSums(cnt_no_dup$X)
  )
)


ggplot(df, aes(cnt_with_dup, cnt_no_dup)) + geom_point()
ggplot(df, aes(cnt_no_dup/cnt_with_dup, cnt_no_dup)) + geom_point()

hist(df$cnt_no_dup/df$cnt_with_dup)

#all(rownames(cnt_with_dup$X)==rownames(cnt_no_dup$X))




################################################################################
##################### Kraken first analysis ####################################
################################################################################

# for(i in 1:20){
#   onemat <- ReadBascetCountMatrix_one(paste0("/husky/henriksson/atrandi//2050905_saliva/kraken.",i,".h5"))
#   print(dim(onemat$X))
#   #print(table(rownames(onemat$X))[table(rownames(onemat$X))>1])
#   print(i)
#   #some cells in more than one matrix
#   print("BASCET__A1_B5_B8_A10" %in% rownames(onemat$X))  #15 and 4
#   print("BASCET__H3_H6_H9_H11" %in% rownames(onemat$X))  #15 and 4
#   #all empty??
#   print("")
#   #add QC, check names unique?
# }


mat <- ReadBascetCountMatrix(bascetRoot,"kraken_mat", verbose=FALSE)
dim(mat$X)
table(rownames(mat$X))[table(rownames(mat$X))>1] ##many are in twice!!!


### Read matrix, rename 
rownames(mat$X) <- stringr::str_remove(rownames(mat$X),stringr::fixed("BASCET_"))
mat$obs$`_index` <- stringr::str_remove(mat$obs$`_index`,stringr::fixed("BASCET_"))

table(rownames(mat$X))[table(rownames(mat$X))>1] ##many are in twice!!!

taxid_ob <- CreateAssayObject(t(mat$X))  ## replace any BASCET_? or keep? seurat: ('_'), replacing with dashes ('-')
adata <- CreateSeuratObject(counts = taxid_ob, project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!

## Add KRAKEN consensus taxonomy to metadata
kraken_taxid <- KrakenFindConsensusTaxonomy(mat$X)
rownames(kraken_taxid) <- kraken_taxid$cell_id
kraken_taxid <- kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
adata@meta.data <- cbind(
  adata@meta.data,
  kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
) #AddMetaData behaves weirdly!!


#### Save metadata for overlay in other sections (support merging different objects too)
saveRDS(kraken_taxid, file.path(bascetRoot,"kraken_metadata.RDS"))
saveRDS(adata@meta.data, file.path(bascetRoot,"adata_kraken_metadata.RDS"))

## Subset data
sum(adata$nCount_RNA>min_nCount_RNA)
adata <- adata[,adata$nCount_RNA>min_nCount_RNA]
adata

#adata <- adata[,adata$species!="Homo sapiens"]  #clearly background!
adata$log_cnt <- log10(1+adata$nCount_RNA)
#adata$perc_human <- adata@assays$RNA$counts["Homo sapiens",]/colSums(adata@assays$RNA$counts)
adata$perc_human <- adata@assays$RNA$counts["taxid-9607",]/colSums(adata@assays$RNA$counts) ###### TODO fix!

## Dimensional reduction using kraken
if(TRUE){
  #ATAC-seq style. think not the best way here
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap")  ## depth seems to be less of a problem here
} else {
  
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:20, reduction.name = "kraken_umap")
}

#saveRDS(adata, "/home/mahogny/github/rbiscvi/adata.RDS")



################################################################################
################## kraken - doublet analysis ###################################
################################################################################

## Using scDblFinder; need to use bioconductors classes
library(SingleCellExperiment)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts=adata@assays$RNA$counts)) #cols should be cell names, rows gene names

library(scDblFinder)
sce <- scDblFinder::scDblFinder(sce)
#9.4% doublets
adata$scDblFinder.score <- sce$scDblFinder.score
adata$scDblFinder.class <- sce$scDblFinder.class

saveRDS(adata, file.path(bascetRoot,"cache_adata_kraken.RDS"))
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
  
FeaturePlot(adata, features = c("scDblFinder.score")) + xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_doublet_score.pdf"), width=7, height=7, limitsize=FALSE) 

DimPlot(adata, group.by = c("scDblFinder.class")) + xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_doublet_class.pdf"), width=7, height=7, limitsize=FALSE) 

VlnPlot(adata, "log_cnt", group.by = "scDblFinder.class")
ggsave(file.path(plotDir, "kraken_doublet_violin.pdf"), width=3, height=7, limitsize=FALSE)


#### Save metadata for overlay in other sections (support merging different objects too)
#saveRDS(kraken_taxid, file.path(bascetRoot,"kraken_metadata.RDS"))
#saveRDS(adata@meta.data, file.path(bascetRoot,"adata_kraken_metadata.RDS"))  ### could subset all cells to these?


################################################################################
################## kraken - plots ##############################################
################################################################################

#####
##### UMAP plots
#####

## UMAP of phylum abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "phylum", reduction = "kraken_umap") + 
#  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_phylum.pdf"), width=7, height=7)

## UMAP of genus abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_genus.pdf"), width=7, height=7)

#For grants etc
DimPlot(object = adata, group.by = "genus", reduction = "kraken_umap") + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")


## UMAP of species abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_species.pdf"), width=7, height=7)



#Compare with depth. "human" cells got few counts and end up in the middle
if(stringr::str_detect(dataset_name,"saliva")){
  FeaturePlot(adata, features = "taxid-339") ## "Xanthomonas campestris" 
  ggsave(file.path(plotDir, "umap_kraken_xanthomonas.pdf"), width=7, height=7, limitsize=FALSE) 
}

FeaturePlot(adata, features = c("log_cnt")) + xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_logcnt.pdf"), width=7, height=7, limitsize=FALSE)

FeaturePlot(adata, features = c("perc_human")) + xlab("KRAKEN1") + ylab("KRAKEN2")        ## it really seems gone in saliva data!
ggsave(file.path(plotDir, "umap_kraken_human.pdf"), width=7, height=7, limitsize=FALSE)


#####
##### Other plots
#####

############# Barplot of species abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select species, count(*) as cnt from df group by species")
df <- df[order(df$cnt),]
df$species <- factor(df$species, levels = df$species)
ggplot(df[df$cnt>5,], aes(species, cnt)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  ylab("Cell count")
ggsave(file.path(plotDir, "kraken_hist_species.pdf"), width=7, height=7, limitsize=FALSE)

############# Barplot of phylum abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select phylum, count(*) as cnt from df group by phylum")
df <- df[order(df$cnt),]
df$phylum <- factor(df$phylum, levels = df$phylum)
ggplot(df[df$cnt>5,], aes(phylum, cnt)) + 
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  ylab("Cell count")
ggsave(file.path(plotDir, "kraken_hist_phylum.pdf"), width=7, height=7, limitsize=FALSE)

############# Barplot of genus abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select genus, count(*) as cnt from df group by genus")
df <- df[order(df$cnt),]
df$genus <- factor(df$genus, levels = df$genus)
ggplot(df[df$cnt>5,], aes(genus, cnt)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() +
  ylab("Cell count")
ggsave(file.path(plotDir, "kraken_hist_genus.pdf"), width=7, height=7, limitsize=FALSE) 

############### Kneeplots for kraken
showNumSpecies <- 10
if(stringr::str_detect(bascetRoot,"saliva")){
  showNumSpecies <- 30
}

KrakenKneePlot(adata, groupby = "phylum", showNumSpecies=showNumSpecies) ############################### why does it not work? NA-value...   STILL NOT WORKING
ggsave(file.path(plotDir, "kraken_kneeplot_per_phylum.pdf"), width=7, height=7, limitsize=FALSE)

KrakenKneePlot(adata, groupby = "genus", showNumSpecies=showNumSpecies)
ggsave(file.path(plotDir, "kraken_kneeplot_per_genus.pdf"), width=7, height=7, limitsize=FALSE)

KrakenKneePlot(adata, groupby = "species", showNumSpecies=showNumSpecies)
ggsave(file.path(plotDir, "kraken_kneeplot_per_species.pdf"), width=7, height=7, limitsize=FALSE)



######## TODO need to redo all KrakenKneePlot!!!

## quick and dirty barplot
df <- data.frame(
  cnt = sort(colSums(adata@assays$RNA@counts), decreasing = TRUE)
)
df$ind <- 1:nrow(df)
ggplot(df, aes(ind, cnt)) + geom_line() + scale_x_log10() + scale_y_log10()


################################################################################
################## Alignment-based analysis (mock only) ########################
################################################################################

if(!str_detect(bascetRoot,"saliva")){
  
  ##########
  ########## Read the aligned counts
  ##########
  
  cnt <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)
  adata <- CreateSeuratObject(
    counts = CreateAssayObject(t(cnt$X)), ## according to anndata standard, there should be a transpose here, unless we abstract this step 
    project = "proj", min.cells = 0, min.features = 0, assay = "chrom_cnt"
  ) 
  
  #min_nCount_RNA <- 10
  ## Dimensional reduction using kraken
  sum(adata$nCount_chrom_cnt>min_nCount_RNA)
  keep_cells <- adata$nCount_chrom_cnt > min_nCount_RNA
  #keep_cells <- adata$nCount_chrom_cnt > sort(adata$nCount_chrom_cnt, decreasing = TRUE)[10000]
  sum(keep_cells)
  adata <- adata[,keep_cells] ## Reduce to sensible number
  adata

  ########## Produce a count matrix on strain level
  adata[["species_cnt"]] <- ChromToSpeciesCount(adata, mapSeq2strain)  #gives warning. coerce ourselves to dgCMatrix
  
  #Figure out which species has most reads in which cell
  cnt <- adata@assays$species_cnt$counts
  adata$species_aln <- rownames(cnt)[apply(cnt, 2, which.max)]
  adata$log_cnt <- log10(1+adata$nCount_chrom_cnt)
  adata$species_aln_short <- stringr::str_split_i(adata$species_aln,stringr::fixed(" ("),1)
  
  
  ##########
  ########## Doublet detection using scDblFinder
  ##########
  
  library(SingleCellExperiment)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=adata@assays$species_cnt$counts)) #cols should be cell names, rows gene names
  
  library(scDblFinder)
  sce <- scDblFinder::scDblFinder(sce)
  #14% doublets, in mock1
  adata$scDblFinder.score <- sce$scDblFinder.score  ### performs not so great
  adata$scDblFinder.class <- sce$scDblFinder.class
  
  
  ##########
  ########## Generate UMAP
  ##########
  
  if(FALSE){
    #ATAC-seq style. think not the best way here
    adata <- RunTFIDF(adata)
    adata <- FindTopFeatures(adata, min.cutoff = 'q0')
    adata <- RunSVD(adata)
    DepthCor(adata)
    adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:(nrow(adata)-1), reduction.name = "chrom_umap")  ## depth seems to be less of a problem here
  } else {
    #RNA-seq style
    adata <- NormalizeData(adata)
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = nrow(adata))
    adata <- ScaleData(adata, features = rownames(adata))
    adata <- RunPCA(adata, features = VariableFeatures(object = adata))
    adata <- RunUMAP(adata, dims = 1:50, reduction.name = "chrom_umap")  #  dims = 1:(nrow(adata)-1),
  }

  ### Cache results  
  saveRDS(adata, file.path(bascetRoot,"cache_adata_align.RDS"))
  adata <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))
  
  if(FALSE){

    #Bacillus
    FeaturePlot(adata, features = mapSeq2strain[mapSeq2strain$strain=="Bacillus pacificus",]$id)
    mean(adata@meta.data$species_aln_short=="Bacillus pacificus")
    #Bacillus pacificus (ATCC 10987) CP086329.1.fasta
    #b-cn1613 [~/mystore/atrandi/bwa_ref/ref10/separate]$ 
    #main genome: 75044
    #plasmid: 4884 --- look for genes that should be on plasmid. check annot
    
    #/husky/henriksson/atrandi/bwa_ref/atcc/ATCC-AGP-24087822e79a496f/genome_ATCC_10987_687931d9b06b4cb4/assembly/
    
    FeaturePlot(adata, features = c(
      #'Bacillus pacificus (ATCC 10987)
      "CP086328.1",
      "CP086329.1" ######## this is the plasmid!
    ))
    
    #nova 3; these do not co-localize as expected
    FeaturePlot(adata, features = c(
      "NC-017316.1",  #Enterococcus faecalis OG1RF, complete sequence
      "NZ-OQ434560.1" #Enterococcus faecalis OG1RF plasmid cfr/cyl, complete sequence
    ))
    
    FeaturePlot(adata, features = c(
        #      ==> Deinococcus radiodurans (ATCC 13939) chr1.fasta <==
        "CP150840.1", # Deinococcus radiodurans R1 = ATCC 13939 = DSM 20539 strain ATCC 13939 chromosome 1, complete sequence #extends
        "CP150841.1", # Deinococcus radiodurans R1 = ATCC 13939 = DSM 20539 strain ATCC 13939 chromosome 2, complete sequence #2 clusters, subset of 1
        "CP150843.1", # Deinococcus radiodurans R1 = ATCC 13939 = DSM 20539 strain ATCC 13939 plasmid pCP1, complete sequence
        "CP150842.1"  # Deinococcus radiodurans R1 = ATCC 13939 = DSM 20539 strain ATCC 13939 plasmid pMP1, complete sequence      
    ))
    
    
    FeaturePlot(adata, features = c(
        #==> Cereibacter sphaeroides (ATCC 17029) chr1.fasta <==
        "NC-009049.1", #Cereibacter sphaeroides ATCC 17029 chromosome 1, complete sequence
        "NC-009050.1", #Cereibacter sphaeroides ATCC 17029 chromosome 2, complete sequence
        "NC-009040.1"  #Cereibacter sphaeroides ATCC 17029 plasmid pRSPH01, complete sequence
    ))
    
    
    
    FeaturePlot(adata, features = c(
      rownames(adata)
    ))

    #adata <- FindNeighbors(adata)#, dims = 1:10)
    #adata <- FindClusters(adata, resolution = 0.05)
    
    #
    #>NC_017316.1 Enterococcus faecalis OG1RF, complete sequence
    #Enterococcus faecalis (ATCC 47077) plasmid.fasta:>NZ_OQ434560.1 Enterococcus faecalis OG1RF plasmid cfr/cyl, complete sequence
    
    #Escherichia coli (ATCC 700926).fasta:>NZ_LR881938.1 Escherichia coli str. K-12 substr. MG1655 strain K-12 chromosome MG1655, complete sequence
  }
  
  
  
  
  ### Plotting - max species per cell
  DimPlot(object = adata, label = TRUE, group.by = "species_aln_short") + 
    xlab("BWA1") + 
    ylab("BWA2") + ggtitle("") #+ NoLegend()
  ggsave(file.path(plotDir,"umap_alignment_species.pdf"), width = 7, height = 5)

  ### Plotting - max species per cell, doublets removed
  DimPlot(object = adata[,adata$scDblFinder.class=="singlet"], label = TRUE, group.by = "species_aln_short") + 
    xlab("BWA1") + 
    ylab("BWA2") + ggtitle("") #+ NoLegend()
  ggsave(file.path(plotDir,"umap_alignment_nodoublet.pdf"), width = 7, height = 5)
  
  # Save distribution for simulation etc
  cnt <- adata@assays$species_cnt$counts
  df_align_bc_dist <- data.frame(
    species=adata$species_aln, #rownames(cnt)[apply(cnt, 2, which.max)],
    cnt=apply(cnt, 2, max)
  )
  saveRDS(df_align_bc_dist, file.path(bascetRoot,"align_metadata.RDS"))
  

  ##########
  ########## Kneeplot for all
  ##########

  plot_one_kneeplot_all <- function(adata){
    df <- data.frame(
      cnt=sort(colSums(adata@assays$species_cnt$counts), decreasing = TRUE)
    )
    df$index <- 1:nrow(df) 
    ggplot(df, aes(index, cnt)) + geom_line() + scale_x_log10() + scale_y_log10() + theme_bw() + ylab("Read count") + xlab("Cell index")
  }
  plot_one_kneeplot_all(adata)
    
  ggsave(file.path(plotDir,"alignment_kneeplot_all.pdf"), width = 6, height = 4)
  
  ##########
  ########## Kneeplot per species
  ##########
  
  DefaultAssay(adata) <- "species_cnt"
  KneeplotPerSpecies(adata)
  ggsave(file.path(plotDir,"alignment_kneeplot.pdf"), width = 6, height = 4)

  ##########
  ########## Barnyard plot
  ##########
  
  #BarnyardPlotMatrix()  this is an option
  
  ### Xanthomonas campestris -- need to remove

  rownames(adata@assays$species_cnt[1:10,])
  adata@assays$species_cnt$counts[1:10,]
  
  bp <- data.frame(
    maxc = MatrixGenerics::colMaxs(adata@assays$species_cnt$counts[1:10,]), ### remove xanthomonas etc
    totc = colSums(adata@assays$species_cnt$counts[1:10,])  ### remove xanthomonas etc
  )
  bp$restc <- bp$totc - bp$maxc
  mean(bp$maxc/bp$totc>0.8)
  
  ggplot(bp, aes(maxc,restc)) + 
    geom_point() +
    theme_bw() +
    xlab("Dominant species count") +
    ylab("Other species count") #+ scale_x_log10() + scale_y_log10()
  ggsave(file.path(plotDir,"alignment_barnyard.pdf"), width = 5, height = 5)
  
  hist(log10(bp$maxc/bp$totc), breaks=100)
  
  ### better? 
  ggplot(bp, aes(maxc,restc)) + 
    geom_point() +
    theme_bw() +
    scale_x_log10()+
    scale_y_log10()+
    xlab("Dominant species count") +
    ylab("Other species count")
  
  #bacillus in all top ones
  adata@meta.data$species_aln[bp$totc>50000]

  adata@assays$species_cnt$counts[,bp$totc>50000]  
  adata@meta.data$species_aln[bp$totc>50000]
  
  
  vscor <- round(digits = 2,cor(t(as.matrix(adata@assays$species_cnt$counts[1:10,]))))  #### any species more than others?
  
  
  #Could store "otherness" as a score to plot
  
  ##########
  ########## Bias plot
  ##########
  
  df <- adata@meta.data
  df <- sqldf("select count(*) as cnt, species_aln_short from df group by species_aln_short")
  df$dataset <- dataset_name
  df
  saveRDS(df, file.path(bascetRoot,"count_per_species_aln.RDS"))
  

}


################################################################################
################## informative KMER-analysis - distribution ####################
################################################################################

p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")
if(file.exists(p_minh)) {
  kmerHist <- BascetReadMinhashHistogram(bascetRoot)
  #kmerHist <- kmerHist[order(kmerHist$cnt, decreasing=TRUE),]
  if(!file.exists(p_use_kmers)) {
    
    kmerHist$rand_index <- sample(1:nrow(kmerHist)) #add a tie breaker
    kmerHist <- kmerHist[order(kmerHist$cnt, kmerHist$rand_index, decreasing=TRUE),]

    #Pick the most common KMERs    
    picked_kmers <- kmerHist$kmer[1:100000]

#     picked_kmers <- ChooseInformativeKMERs( ### 0.002 => 11212 this killed conda!
#       kmerHist,
# #      minFreq = 0.2 #simulated2; much higher cutoff needed for 34k 
# #      minFreq = 0.002 #novaseq1, really few!
# #      minFreq = 0.0015 #saliva, really few!
#       minFreq = 0.001 #novaseq2, really few!
# #      minFreq = 0.0005 #novaseq3, really few!
#       ) ################################################### option of max count?
    writeLines(picked_kmers, p_use_kmers)
    bascetRoot
  }
  
  
  
  ### KMER count histogram
  kmerHist$rank <- 1:nrow(kmerHist)
  ggplot(kmerHist[sample(1:nrow(kmerHist),min(30000, nrow(kmerHist))),], aes(rank, cnt)) +   # not>2!!
    #ggplot(kmerHist[sample(which(kmerHist$cnt>10), 5000),], aes(rank, cnt)) +   # not>2!!
    geom_point() + 
    scale_x_log10() + 
    scale_y_log10() +
    theme_bw() +
    ylab("Count")+
    xlab("Rank")
  ggsave(file.path(plotDir,"info_kmerHist.pdf"), width = 5, height = 5)
  ggsave(file.path(plotDir,"info_kmerHist.png"), width = 5, height = 5)
  
}



################################################################################
################## informative KMER-analysis - umap etc ########################
################################################################################


### Read count matrix
cnt <- ReadBascetCountMatrix(  #  x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix  ---- do this conversion manually
  bascetRoot,
  "kmer_counts"
) 
cnt <- cnt$X
dim(cnt)
#rownames(cnt) <- paste0("_",rownames(cnt))  #### todo fix naming

if(FALSE){
  ### Order by frequency of feature
  cnt <- cnt[,order(colSums(cnt), decreasing = TRUE)]
  
  ### Subset features (speed up, and to test how many we need)
  cnt <- cnt[,1:10000]
}



###### Comparison of abundance, histogram vs query
if(FALSE){
  sum(colSums(cnt)>100000)
  
  ck_hist <- kmerHist[kmerHist$kmer %in% picked_kmers,]
  ck_hist <- kmerHist#[kmerHist$kmer %in% picked_kmers,]
  colnames(ck_hist) <- c("kmer","cnt_hist")
  
  ck_query <- data.frame(
    kmer=colnames(cnt),
    cnt_q = colSums(cnt),
    cnt_q0 = colSums(cnt>0)
  )
  
  ck_merge <- merge(ck_query, ck_hist)
  ggplot(ck_merge, aes(cnt_q, cnt_hist)) + geom_point() + xlim(0,1e8) #+ ylim(0,1e8)  ######### how can correlation be so bad??
  ggplot(ck_merge, aes(cnt_q0, cnt_hist)) + geom_point() ######### how can correlation be so bad??  ... likely because we don't full depth
}


if(TRUE) {

  cnt <- cnt[order(rowSums(cnt), decreasing = TRUE),]
  #rowSums(cnt)>20000
  
  ## Load subset of real cells
  #keep_cells <- rowSums(cnt)>20000
  #sum(keep_cells)
  adata <- CreateSeuratObject(
    counts = CreateAssayObject(t(cnt[1:10000,])),  #cutoff on cell count
    #counts = CreateAssayObject(t(cnt[keep_cells,])),  #cutoff on kmer count
    assay = "infokmer"
  )
  
} else {
  
  ### Alternative subsetting, and binarization
  sum(rowSums(cnt)>20000)
  adata <- CreateSeuratObject(
    counts = CreateAssayObject(t(cnt[rowSums(cnt)>20000,]>0)),
    assay = "infokmer"
  )
  adata
  
  plot(sort(as.integer(colSums(cnt))))
  
  sum(rowSums(cnt)>20000)
  adata <- CreateSeuratObject(
    counts = CreateAssayObject(t(cnt[
      rowSums(cnt)>20000,
      colSums(cnt)<100000
    ]>0)),
    assay = "infokmer"
  )
  adata  
}


##### Add kraken metadata
kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
adata@meta.data <- cbind(
  adata@meta.data,
  kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
) #AddMetaData behaves weirdly!!

##### Perform ATAC-seq style analysis
adata$log10_nCount_infokmer <- log10(adata$nCount_infokmer)
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)

DepthCor(adata)
ggsave(file.path(plotDir,"infokmer_depthcor.pdf"), width = 4, height = 3)

adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "infokmers_umap")  ## dim 1 affected plenty if binarizing matrix

saveRDS(adata, file.path(bascetRoot,"cache_adata_infokmer.RDS"))

##### 
##### Plots
##### 

DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDir,"infokmer_umap_genus.pdf"), width = 5, height = 5)

DimPlot(object = adata, label = TRUE, group.by = "phylum", reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDir,"infokmer_umap_phylum.pdf"), width = 5, height = 5)

FeaturePlot(adata, "log10_nCount_infokmer", reduction = "infokmers_umap") + xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDir,"infokmer_umap_kmercount.pdf"), width = 5, height = 5)


#DimPlot(object = adata[,adata$genus=="Bacillus"], label = TRUE, group.by = "genus", reduction = "infokmers_umap") + NoLegend()
#DimPlot(object = adata, label = TRUE, reduction = "infokmers_umap") + NoLegend()
#table(adata$genus) #MOCK: cereibacter dominates plenty!!

DimPlot(
  object = adata[,adata$genus %in% names(which(table(adata$genus)>50))], 
  label = TRUE, group.by = "genus", reduction = "infokmers_umap")  + xlab("IK1") + ylab("IK2")
  #NoLegend()
ggsave(file.path(plotDir,"infokmer_umap_genus_dominant.pdf"), width = 7, height = 7)



################################################################################
################## Count sketch-based analysis #################################
################################################################################

#Load data as seurat object 
adata <- BascetLoadCountSketchMatrix(bascetRoot,inputName = "countsketch_mat.csv") ############# TODO: load full set of shards. negative counts!!
#adata <- BascetLoadCountSketchMatrix(bascetRoot,inputName = "countsketch.0.tsv") ############# update!
adata$log10_celldepth <- log10(adata$celldepth)

#Subset cells
min_kmers <- 400000  #40k for saliva => 8k cells; option: pick top N cells for comparability
min_kmers <- 1000  #40k for saliva => 8k cells; option: pick top N cells for comparability
keep_cells <- adata$celldepth > min_kmers
sum(keep_cells)
adata <- adata[,keep_cells]
adata <- adata[,1:5000]

#Add kraken metadata
kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
adata@meta.data <- cbind(
  adata@meta.data,
  kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
) #AddMetaData behaves weirdly!!

#Non-linear dimensional reduction
adata <- RunUMAP(adata, dims = 1:ncol(adata@reductions$kmersketch@cell.embeddings), reduction = "kmersketch")  #Searching Annoy index using 1 thread, search_k = 3000 ; can do more

#any(is.na(adata@reductions$kmersketch@cell.embeddings))



saveRDS(adata, file.path(bascetRoot,"cache_adata_cs.RDS"))

##### 
##### Plotting
##### 

DimPlot(object = adata, label = TRUE) + NoLegend() + xlab("CS1") + ylab("CS2")

DimPlot(object = adata, label = TRUE, group.by = "phylum") + NoLegend() + xlab("CS1") + ylab("CS2")
ggsave(file.path(plotDir,"umap_countsketch_phylum.pdf"), width = 7, height = 7)

DimPlot(object = adata[,adata$genus!="Streptococcus"], label = TRUE, group.by = "phylum") + NoLegend() + xlab("CS1") + ylab("CS2")
ggsave(file.path(plotDir,"umap_countsketch_phylum_nostrepto.pdf"), width = 7, height = 7)

DimPlot(object = adata, label = TRUE, group.by = "genus") + NoLegend() + xlab("CS1") + ylab("CS2")
ggsave(file.path(plotDir,"umap_countsketch_genus.pdf"), width = 7, height = 7)

DimPlot(object = adata, label = TRUE, group.by = "species") + NoLegend() + xlab("CS1") + ylab("CS2")
ggsave(file.path(plotDir,"umap_countsketch_species.pdf"), width = 7, height = 7)

FeaturePlot(adata, features = "log10_celldepth") + xlab("CS1") + ylab("CS2")
ggsave(file.path(plotDir,"umap_countsketch_logkmer.pdf"), width = 7, height = 7)




################################################################################
######### Analysis of assembled genomes: quast #################################
################################################################################

#### 
#### Load Quast and merge with KRAKEN analysis
#### 

adata_krakenmeta <- readRDS(file.path(bascetRoot,"adata_kraken_metadata.RDS"))
#rownames(adata_krakenmeta) <- paste0("_",rownames(adata_krakenmeta)) ##need to fix
#nrow(adata_krakenmeta)

#TODO: put MapListAsDataFrame into quast aggr 
quast_aggr.df <- MapListAsDataFrame(readRDS(file.path(bascetRoot,"cache_quast.RDS")))
paste(nrow(quast_aggr.df),"vs",nrow(adata_krakenmeta)) #check that overlap is similar

## perform below in there too
quast_aggr.df$tot_length <- as.integer(quast_aggr.df$`Total length`)
quast_aggr.df$largest_contig <- as.integer(quast_aggr.df$`Largest contig`)
quast_aggr.df$num_contigs <- as.integer(quast_aggr.df$`Number of contigs`)
quast_aggr.df$N50 <- as.integer(quast_aggr.df$N50)

## Include seq depth, from KRAKEN analysis
quast_aggr.df$depth <- adata_krakenmeta[rownames(quast_aggr.df),]$nCount_RNA ### poor overlap for saliva data. weird!!

## Only consider cells with full info
sum(is.na(quast_aggr.df$depth)) #should be 0 or close to
quast_aggr.df <- quast_aggr.df[!is.na(quast_aggr.df$depth),]

#### 
#### Plotting
#### 

p1 <- ggplot(quast_aggr.df, aes(tot_length, largest_contig, color=log10(depth))) +  
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  xlab("Total contig length") + 
  ylab("Largest contig")
p1
ggsave(plot = p1, file.path(plotDir,"quast_totlen_VS_largestContig.png"), width = 5, height = 5)

p2 <- ggplot(quast_aggr.df, aes(tot_length, num_contigs, color=log10(depth))) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  xlab("Total contig length") + 
  ylab("Number of contigs")
p2
ggsave(plot = p2, file.path(plotDir,"quast_totlen_VS_numContigs.png"), width = 5, height = 5)

p3 <- ggplot(quast_aggr.df, aes(N50, num_contigs, color=log10(depth))) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  xlab("N50") +  #N50 ~ weighted median statistic such that 50% of the entire assembly is contained in contigs or scaffolds equal to or larger than this value.
  ylab("Number of contigs")
p3
ggsave(plot = p3, file.path(plotDir,"quast_n50_VS_numContigs.png"), width = 5, height = 5)



ptot <- egg::ggarrange(p1,p2,p3, nrow=1)
ggsave(plot = ptot, file.path(plotDir,"quast_allqc.png"), width = 12, height = 3)
ggsave(plot = ptot, file.path(plotDir,"quast_allqc.svg"), width = 12, height = 3)

################################################################################
######### Analysis of assembled genomes: abricate ##############################
################################################################################


#### Load UMAP coordinates
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
#colnames(adata)

#### Load abricate data
abricate <- readRDS(file.path(bascetRoot,"cache_abricate.RDS"))
#rownames(abricate)

### Put into adata metadata
#cellid_abricate <- abricate[rownames(adata@meta.data),] #there is no _ for saliva!!!
#cellid_abricate <- abricate[paste0("_",rownames(adata@meta.data)),] #TODO quick hack
cellid_abricate <- abricate[paste0(rownames(adata@meta.data)),] #subset
cellid_abricate[is.na(cellid_abricate)] <- 0
#colSums(cellid_abricate) #few left -- note, returns 0 if matrix is sparse(!!)

for(i in 1:ncol(cellid_abricate)) {
  print(i)
  adata@meta.data[,paste0("abricate_", colnames(cellid_abricate)[i])] <- as.character(cellid_abricate[,i])
}

### Add collapsed category. Will emphasize one category right now; concatenate? put names instead?
collapsed_abrigate <- colnames(cellid_abricate)[apply(cellid_abricate,1,function(x) which(x>0)[1])]
collapsed_abrigate[is.na(collapsed_abrigate)] <- "_"
adata@meta.data$collapsed_abrigate <- collapsed_abrigate


###############################################
#' TODO generalize this
BascetAddSeuratMetadataDF <- function(adata, newmeta, prefix=NULL, default_val="NA") {
  
  colnames(newmeta) <- stringr::str_replace_all(colnames(newmeta), stringr::fixed(" "),"_")
  newmeta <- newmeta[rownames(adata@meta.data),]
  newmeta[is.na(newmeta)] <- default_val  ### hmmmm...
  
  if(is.null(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix,"_")
  }

  print("Adding:")
  print(paste0(prefix, colnames(newmeta)))
  for(i in 1:ncol(newmeta)) {
    adata@meta.data[,paste0(prefix, colnames(newmeta)[i])] <- newmeta[,i]
#    adata@meta.data[,paste0(prefix, colnames(newmeta)[i])] <- as.character(newmeta[,i]) ### change
  }
  
  adata
}


df <- data.frame(
  x=adata@reductions$kraken_umap@cell.embeddings[,1],
  y=adata@reductions$kraken_umap@cell.embeddings[,2],
  abricate=adata$collapsed_abrigate
)
df <- df[order(df$abricate),]
ggplot(df, aes(x,y,color=abricate)) + 
  geom_point() + 
  theme_bw() +
  scale_colour_manual(values = c("lightgray","blue","red","black","yellow","magenta","darkgreen","cyan","darkgrey","lightblue")) +
  xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir,"abricate_collapsed.pdf"), width = 5, height = 5)


#Summary of where things are
table(adata@meta.data[,c("collapsed_abrigate","genus")])
table(adata@meta.data$collapsed_abrigate)


################################################################################
######### Analysis of reads: FASTQC ############################################
################################################################################

#### Load UMAP coordinates
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))

#### Load FASTQC data and impose on UMAP object
aggr_fastqc <- readRDS(file.path(bascetRoot,"cache_fastqc.RDS"))  #saliva is broken

#rownames(fastqc_stat_passfail_r1) <- stringr::str_sub(rownames(fastqc_stat_passfail_r1),2)  #hack for saliva
#rownames(fastqc_stat_passfail_r2) <- stringr::str_sub(rownames(fastqc_stat_passfail_r2),2)  #hack for saliva

fastqc_basic_stats <- GetFASTQCbasicStats(aggr_fastqc,2)  ### something is off for saliva. why no READ1 ????????????????????????
adata <- BascetAddSeuratMetadataDF(adata, fastqc_basic_stats, prefix = "fastqc")

if(FALSE){
  #Not the most useful stats
  fastqc_stat_passfail_r1 <- GetFASTQCpassfailStats(aggr_fastqc, "1")  #saliva got _ here; need to remove  
  fastqc_stat_passfail_r2 <- GetFASTQCpassfailStats(aggr_fastqc, "2")  #saliva got _ here; need to remove
  adata <- BascetAddSeuratMetadataDF(adata, fastqc_stat_passfail_r1, prefix = "fastqc_r1")
  adata <- BascetAddSeuratMetadataDF(adata, fastqc_stat_passfail_r2, prefix = "fastqc_r2")
}


#### 
#### Plotting
#### 

### fix. broken from BascetAddSeuratMetadataDF
adata$fastqc_gc <- as.numeric(adata$fastqc_gc)
adata$fastqc_num_seq_poor_quality <- as.numeric(adata$fastqc_num_seq_poor_quality)
adata$fastqc_seqlen_from <- as.numeric(adata$fastqc_seqlen_from)
adata$fastqc_seqlen_to <- as.numeric(adata$fastqc_seqlen_to)

FeaturePlot(adata, features = "fastqc_gc", reduction = "kraken_umap", label = TRUE) + xlab("KRAKEN1") + ylab("KRAKEN2") + ggtitle("GC%")
ggsave(file.path(plotDir,"umap_kraken_fastqc_GC.pdf"), width = 5, height = 5)










################################################################################
########################## Bias plot ###########################################
################################################################################

df <- rbind(
  readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq1/count_per_species_aln.RDS"),
  readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/count_per_species_aln.RDS")
)

#Normalize cell count
for(ds in unique(df$dataset)){
  df$cnt[df$dataset==ds] <- df$cnt[df$dataset==ds]/sum(df$cnt[df$dataset==ds])*100
}

ggplot(df, aes(species_aln_short, cnt, fill=dataset)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() + 
  theme_bw() +
  ylab("Fraction (%)") +
  xlab("")+
  scale_x_discrete(limits=rev)
ggsave(file.path(plotDirAll, "biascomparison.svg"), width = 7, height = 3)










################################################################################
##################### saliva per species kneeplot ##############################
################################################################################

adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))

KrakenKneePlot(adata, groupby = "genus", showNumSpecies=20)
ggsave(file.path(plotDirAll, "kraken_saliva_kneeplot_per_genus.pdf"), width=7, height=7, limitsize=FALSE)


KrakenKneePlot(adata, groupby = "species", showNumSpecies=10)
#Streptococcus gordonii


sort(unique(adata$species))



