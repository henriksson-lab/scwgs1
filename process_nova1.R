### preprocessing of data done in separate file, prep_*.R

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


library(Seurat)
library(Signac)
library(sqldf)
library(stringr)
library(ggplot2)
library(Matrix)

################################################################################
################## Postprocessing with Bascet/Zorn #############################
################################################################################

dataset_name <- "v4_wgs_novaseq1"
dataset_name <- "v4_wgs_saliva1"



bascet_runner <- LocalRunner(direct = TRUE)
bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)
plotDir <- file.path("~/github/scwgs/plots/",dataset_name)
if(!file.exists(plotDir)){
  dir.create(plotDir, recursive = TRUE)
  print("made plot dir")
}

#bascetRoot = "/husky/henriksson/atrandi/v4_wgs_novaseq1/"
#bascetRoot = "/husky/henriksson/atrandi/v4_wgs_novaseq3/"
#bascetRoot = "/husky/henriksson/atrandi/v4_wgs_saliva1//"
#bascetRoot = "/husky/henriksson/atrandi/v2_wgs_miseq2/"  #for development

#100 for miseq??
min_nCount_RNA <- 1000

################################################################################
################## Kraken-based analysis #######################################  
################################################################################

setwd("/home/mahogny/github/zorn") #to put SQL in the right place

### Read matrix, rename 
mat <- ReadBascetKrakenMatrix(bascetRoot, "kraken")
colnames(mat) <- stringr::str_remove(colnames(mat),stringr::fixed("BASCET_"))

### Compress the representation to avoid trouble with some tools
compressed_mat <- SetTaxonomyNamesFeatures(mat) ########################## didnt load until seurat was loaded?? rowSums in Matrix package!
### works when calling code outside zorn!!


taxid_ob <- CreateAssayObject(compressed_mat)  ## replace any BASCET_? or keep? seurat: ('_'), replacing with dashes ('-')
adata <- CreateSeuratObject(counts = taxid_ob, project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!

map_cellid_depth <- data.frame(
  row.names = colnames(adata),
  cellid=colnames(adata),
  depth=adata$nCount_RNA
)


## Add KRAKEN consensus taxonomy to metadata
kraken_taxid <- KrakenFindConsensusTaxonomy(mat)
rownames(kraken_taxid) <- kraken_taxid$cell_id
kraken_taxid <- kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]) #AddMetaData behaves weirdly!!


#### How many species?
##################KrakenSpeciesDistribution(adata)  ### Could use metadata column! TODO   rewrite

saveRDS(adata@meta.data, file.path(bascetRoot,"kraken_metadata.RDS"))


## Dimensional reduction using kraken
#DefaultAssay(adata) <- "kraken"
sum(adata$nCount_RNA>min_nCount_RNA)
adata <- adata[,adata$nCount_RNA>min_nCount_RNA] ## Reduce to sensible number
#sum(adata$nCount_RNA>1000)
#adata <- adata[,adata$nCount_RNA>1000] ## Reduce to sensible number
adata

#adata <- adata[,adata$species!="Homo sapiens"]  #clearly background!
adata$log_cnt <- log10(1+adata$nCount_RNA)
adata$perc_human <- adata@assays$RNA$counts["Homo sapiens",]/colSums(adata@assays$RNA$counts)

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
#DimPlot(object = adata, label = TRUE) + NoLegend()


############# Barplot of species abundance according to KRAKEN

DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap")
ggsave(file.path(plotDir, "kraken_umap_genus.pdf"), width=15, height=5)

DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap")
ggsave(file.path(plotDir, "kraken_umap_species.pdf"), width=15, height=5)


############# Barplot of species abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select species, count(*) as cnt from df group by species")
df <- df[order(df$cnt),]
df$species <- factor(df$species, levels = df$species)
ggplot(df[df$cnt>5,], aes(species, cnt)) + geom_bar(stat="identity") + coord_flip() + theme_bw()
ggsave(file.path(plotDir, "hist_species.pdf"), width=7, height=10, limitsize=FALSE) ### did we not remove xanthomonas by alignment, saliva?? TODO

#saveRDS(adata, "/husky/henriksson/atrandi/wgs_saliva1/kraken.RDS")


#Compare with depth. "human" cells got few counts and end up in the middle

FeaturePlot(adata, features = "Xanthomonas campestris")
ggsave(file.path(plotDir, "kraken_xanthomonas.pdf"), width=7, height=10, limitsize=FALSE)
FeaturePlot(adata, features = c("log_cnt","perc_human"))
ggsave(file.path(plotDir, "kraken_perc_human.pdf"), width=14, height=10, limitsize=FALSE)


## Picks the wrong ones
KneeplotPerSpecies(adata, max_species = 10) ##does not work for saliva?? TODO 66666666666666666666666666666666666666666666
ggsave(file.path(plotDir, "kraken_kneeplot_per_species.pdf"), width=7, height=10, limitsize=FALSE)

#NOTE: 3 cereibacter??? can we collapse? TODO


################################################################################
################## kraken - doublet analysis ################################### TODO
################################################################################



























################################################################################
################## Alignment-based analysis ####################################  
################################################################################


cnt <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)

adata <- CreateSeuratObject(counts = CreateAssayObject(cnt), project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!


## Dimensional reduction using kraken
#DefaultAssay(adata) <- "kraken"
sum(adata$nCount_RNA>min_nCount_RNA)
adata <- adata[,adata$nCount_RNA>min_nCount_RNA] ## Reduce to sensible number
adata

adata$log_cnt <- log10(1+adata$nCount_RNA)

if(FALSE){
  #ATAC-seq style. think not the best way here
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:(nrow(adata)-1), reduction.name = "chrom_umap")  ## depth seems to be less of a problem here
} else {
  
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = nrow(adata))
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:(nrow(adata)-1), reduction.name = "chrom_umap")
}
DimPlot(object = adata, label = TRUE) + NoLegend()




### TODO label cell by maximum species.
### TODO kneeplot
### TODO doublet removal
### TODO barnyard plot





################################################################################
################## informative KMER-analysis ################################### plotting of distribution
################################################################################

p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")
if(file.exists(p_minh)) {
  kmer_hist <- BascetReadMinhashHistogram(bascetRoot)
  if(!file.exists(p_use_kmers)) {
    picked_kmers <- ChooseInformativeKMERs(kmer_hist)
    writeLines(picked_kmers, p_use_kmers)
  }
}


#kmer_hist <- BascetReadMinhashHistogram(bascetRoot)


### Figure 4b
kmer_hist$rank <- 1:nrow(kmer_hist)
ggplot(kmer_hist[sample(1:nrow(kmer_hist),min(30000, nrow(kmer_hist))),], aes(rank, cnt)) +   # not>2!!
  #ggplot(kmer_hist[sample(which(kmer_hist$cnt>10), 5000),], aes(rank, cnt)) +   # not>2!!
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("Count")+
  xlab("Rank")
ggsave(file.path(plotDir,"info_kmer_hist.pdf"), width = 5, height = 5)
ggsave(file.path(plotDir,"info_kmer_hist.png"), width = 5, height = 5)






################################################################################
################## Informative KMER-based analysis #############################  umap and stuff
################################################################################








################################################################################
################## Analysis of assembled genomes ###############################
################################################################################


#TODO: put MapListAsDataFrame into quast aggr 
quast_aggr.df <- MapListAsDataFrame(readRDS(file.path(bascetRoot,"cache_quast.RDS")))

## perform below in there too
quast_aggr.df$tot_length <- as.integer(quast_aggr.df$`Total length`)
quast_aggr.df$largest_contig <- as.integer(quast_aggr.df$`Largest contig`)
quast_aggr.df$num_contigs <- as.integer(quast_aggr.df$`Number of contigs`)
quast_aggr.df$N50 <- as.integer(quast_aggr.df$N50)

## Include seq depth
quast_aggr.df$depth <- map_cellid_depth[rownames(quast_aggr.df),]$depth

ggplot(quast_aggr.df, aes(tot_length, largest_contig, color=log10(depth))) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw()
ggsave(file.path(plotDir,"quast_totlen_VS_largestContig.pdf"), width = 5, height = 5)

ggplot(quast_aggr.df, aes(tot_length, num_contigs, color=log10(depth))) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw()
ggsave(file.path(plotDir,"quast_totlen_VS_numContigs.pdf"), width = 5, height = 5)

ggplot(quast_aggr.df, aes(N50, num_contigs, color=log10(depth))) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw()
ggsave(file.path(plotDir,"quast_n50_VS_numContigs.pdf"), width = 5, height = 5)



###############
############### Abricate
###############


abricate <- readRDS(file.path(bascetRoot,"cache_abricate.RDS"))

### Put into adata metadata
#cellid_abricate <- abricate[rownames(adata@meta.data),] #there is no _ for saliva!!!
cellid_abricate <- abricate[paste0("_",rownames(adata@meta.data)),] #TODO quick hack
cellid_abricate[is.na(cellid_abricate)] <- 0
colSums(cellid_abricate) #few left
#cellid_abricate

for(i in 1:ncol(cellid_abricate)) {
  print(i)
  adata@meta.data[,paste0("abricate_", colnames(cellid_abricate)[i])] <- as.character(cellid_abricate[,i])
}



BascetAddSeuratMetadataDF <- function(adata, newmeta, prefix=NULL, default_val="NA") {
  
  colnames(newmeta) <- stringr::str_replace_all(colnames(newmeta), stringr::fixed(" "),"_")
#  colnames(newmeta)
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


if(TRUE){
  
  ################### MOCK COMMUNITY
  
  
  #one cluster: bacillus
  if(FALSE){
    DimPlot(adata, group.by = "abricate_BcII", reduction = "kraken_umap")
    DimPlot(adata, group.by = "abricate_bla1", reduction = "kraken_umap")
    DimPlot(adata, group.by = "abricate_fosB_gen", reduction = "kraken_umap")
    DimPlot(adata, group.by = "abricate_lnu(D)", reduction = "kraken_umap") 
    DimPlot(adata, group.by = "abricate_satA_Ba", reduction = "kraken_umap")
    DimPlot(adata, group.by = "abricate_vanZ-F", reduction = "kraken_umap")
  }
  DimPlot(adata, group.by = c(
    "abricate_BcII",
    "abricate_bla1",
    "abricate_fosB_gen", 
    "abricate_lnu(D)", 
    "abricate_satA_Ba",
    "abricate_vanZ-F"
  ), reduction = "kraken_umap")
  ggsave(file.path(plotDir,"abricate_bacillus.pdf"), width = 15, height = 15)
  
  
  
  #upright cluster: enterococcus
  if(FALSE){
    DimPlot(adata, group.by = "abricate_dfrE", reduction = "kraken_umap") 
    DimPlot(adata, group.by = "abricate_lsa(A)", reduction = "kraken_umap") 
  }
  DimPlot(adata, group.by = c(
    "abricate_dfrE",
    "abricate_lsa(A)"
  ), reduction = "kraken_umap")
  ggsave(file.path(plotDir,"abricate_enterococcus.pdf"), width = 5, height = 5)
  
  
  #human cluster? homo!! (really!)
  DimPlot(adata, group.by = "abricate_blaEC", reduction = "kraken_umap")  
  ggsave(file.path(plotDir,"abricate_homo.pdf"), width = 5, height = 5)
  
  
  #cereibacter -- very few
  if(FALSE){
    DimPlot(adata, group.by = "abricate_tet(K)", reduction = "kraken_umap") 
    ggsave(file.path(plotDir,"abricate_cerei.pdf"), width = 5, height = 5)
  }
  
  
  ########## To correlates with 
  DimPlot(adata, group.by = "genus", reduction = "kraken_umap", label = TRUE) 
  
  
  ### Some other correlation table display. Heatmap?
  table(adata@meta.data[,c("genus","abricate_tet(K)")])
  table(adata@meta.data[,c("genus","abricate_blaEC")])
  
  
  
} else {

  ################### saliva
  
  DimPlot(adata, group.by = paste0("abricate_",colnames(abricate)), reduction = "kraken_umap") 
  ggsave(file.path(plotDir,"abricate_saliva.pdf"), width = 30, height = 30)
  
  ###### barely anything
  
  #table(abricate)
  
}


#this one stinks
table(adata@meta.data[,c(
  paste0("abricate_", colnames(cellid_abricate)),
  "genus")]
)
  

###############
############### FASTQC
###############

aggr_fastqc <- readRDS(file.path(bascetRoot,"cache_fastqc.RDS"))


fastqc_stat_passfail_r1 <- GetFASTQCpassfailStats(aggr_fastqc, "1")  #saliva got _ here; need to remove

rownames(fastqc_stat_passfail_r1) <- stringr::str_sub(rownames(fastqc_stat_passfail_r1),2)  #hack for saliva
colnames(adata)

adata <- BascetAddSeuratMetadataDF(adata, fastqc_stat_passfail_r1, prefix = "fastqc")

if(FALSE){
  FeaturePlot(adata, features = "fastqc_Per sequence GC content", reduction = "kraken_umap", label = TRUE) 
  FeaturePlot(adata, features = "fastqc_Sequence Duplication Levels", reduction = "kraken_umap", label = TRUE) 
  FeaturePlot(adata, features = "fastqc_Sequence Length Distribution", reduction = "kraken_umap", label = TRUE) 
  FeaturePlot(adata, features = "fastqc_Sequence Duplication Levels", reduction = "kraken_umap", label = TRUE) 
  
  adata@meta.data$fastqc_Per_sequence_GC_content
  DimPlot(adata, group.by = "fastqc_Per_sequence_GC_content", reduction = "kraken_umap", label = TRUE) 
}

#FeaturePlot(adata, features = "fastqc_Adapter Content", reduction = "kraken_umap", label = TRUE) 
DimPlot(adata, group.by = "fastqc_Adapter_Content", reduction = "kraken_umap", label = TRUE) 
ggsave(file.path(plotDir,"fastqc_fastqc_pf_adaptercontent.pdf"), width = 5, height = 5)























################################################################################
################## Host % DNA ##################################################
################################################################################



dat <- read.table(file.path(bascetRoot,"alignment_stats.csv"),sep="\t")
colnames(dat) <- c("chrom","seqlen","count","count2")
dat$organism <- "Human"
dat$organism[stringr::str_starts(dat$chrom,stringr::fixed("NZ_"))] <- "Xanthomonas"
dat$organism[dat$chrom=="*"] <- "Unmapped"
#dat$organism[dat$organism=="*"] <- "Unmapped"

dat$count[dat$organism=="Unmapped"] <- dat$count2[dat$organism=="Unmapped"]

dat_cnt <- sqldf::sqldf("select sum(count) as count, organism from dat group by organism")
dat_cnt$frac <- dat_cnt$count / sum(dat_cnt$count) *100


ggplot(dat_cnt, aes(frac, organism)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() +
  xlab("% reads") + 
  ylab("")

ggsave(file.path(plotDir,"host_fraction.pdf"), width = 3, height = 4)


