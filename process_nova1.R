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

################################################################################
################## Postprocessing with Bascet/Zorn #############################
################################################################################

bascet_runner <- LocalRunner(direct = TRUE)
bascetRoot = "/husky/henriksson/atrandi/v4_wgs_novaseq1/"
bascetRoot = "/husky/henriksson/atrandi/v4_wgs_novaseq3/"
bascetRoot = "/husky/henriksson/atrandi/v4_wgs_saliva1//"
#bascetRoot = "/husky/henriksson/atrandi/v2_wgs_miseq2/"  #for development


################################################################################
################## Kraken-based analysis #######################################  
################################################################################


mat <- ReadBascetKrakenMatrix(bascetRoot, "kraken")


### Compress the representation to avoid trouble with some tools
compressed_mat <- SetTaxonomyNamesFeatures(mat)

taxid_ob <- CreateAssayObject(compressed_mat)
adata <- CreateSeuratObject(counts = taxid_ob, project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!


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
sum(adata$nCount_RNA>100)
adata <- adata[,adata$nCount_RNA>100] ## Reduce to sensible number
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


DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap")
#ggsave("/husky/henriksson/atrandi/wgs_saliva1/umap_genus.pdf", width=15, height=5)

DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap")
#ggsave("/husky/henriksson/atrandi/wgs_saliva1/umap_species.pdf", width=40, height=10, limitsize=FALSE)

df <- adata@meta.data
df <- sqldf("select species, count(*) as cnt from df group by species")
df <- df[order(df$cnt),]
df$species <- factor(df$species, levels = df$species)
ggplot(df[df$cnt>5,], aes(species, cnt)) + geom_bar(stat="identity") + coord_flip() + theme_bw()
ggsave("/husky/henriksson/atrandi/wgs_saliva1/hist_species.pdf", width=7, height=10, limitsize=FALSE)




#saveRDS(adata, "/husky/henriksson/atrandi/wgs_saliva1/kraken.RDS")


#Compare with depth. "human" cells got few counts and end up in the middle

FeaturePlot(adata, features = "log_cnt")
FeaturePlot(adata, features = "Xanthomonas campestris")
FeaturePlot(adata, features = "perc_human")

## Picks the wrong ones
KneeplotPerSpecies(adata, max_species = 10)

#v2_wgs_novaseq1  xanthomonas also here??? in mock??


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
sum(adata$nCount_RNA>100)
adata <- adata[,adata$nCount_RNA>100] ## Reduce to sensible number
#sum(adata$nCount_RNA>1000)
adata <- adata[,adata$nCount_RNA>1000] ## Reduce to sensible number
adata <- adata[,adata$nCount_RNA>5000] ## Reduce to sensible number
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
################## informative KMER-analysis ###################################
################################################################################


kmer_hist <- BascetReadMinhashHistogram(bascetRoot)


### Figure 4b
kmer_hist$rank <- 1:nrow(kmer_hist)
ggplot(kmer_hist[kmer_hist$cnt>10,], aes(rank, cnt)) +   # not>2!!
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("Count")+
  xlab("Rank")
ggsave(file.path(bascetRoot,"info_kmer_hist.pdf"), width = 5, height = 5)
ggsave(file.path(bascetRoot,"info_kmer_hist.png"), width = 5, height = 5)


picked_kmers <- ChooseInformativeKMERs(kmer_hist)

writeLines(picked_kmers, file.path(bascetRoot,"use_kmers.txt"))


#todo mapcell: check if script can do threads. 
#1. if so, run single-threaded and give all threads to the script
#2. skip merge at end if a single thread

#todo mapcell script: 
# reader responsible for extracting the data needed!!








################################################################################
################## Analysis of assembled genomes ###############################
################################################################################


#BascetCacheComputation(bascetRoot,"saved_quast_agggr",f2(1,2))

#/husky/henriksson/atrandi/bar

######### Aggregate data from quast mapcell call 
quast_aggr.df <- BascetCacheComputation(bascetRoot,"saved_quast_agggr",MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  "quast",
  aggr.quast,
  verbose=TRUE
  #  include_cells = c(cellname_coord$cell[1:500])
)))
#quast_aggr.df <- (quast_aggr)


quast_aggr.df <- as.data.frame(quast_aggr.df)
quast_aggr.df$tot_length <- as.integer(quast_aggr.df$`Total length`)
quast_aggr.df$largest_contig <- as.integer(quast_aggr.df$`Largest contig`)
quast_aggr.df$num_contigs <- as.integer(quast_aggr.df$`Number of contigs`)
quast_aggr.df$N50 <- as.integer(quast_aggr.df$N50)


ggplot(quast_aggr.df, aes(tot_length, largest_contig)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()

ggplot(quast_aggr.df, aes(tot_length, num_contigs)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()


plot(log10(sort(as.integer(quast_aggr.df$tot_length))))
plot((sort(as.integer(quast_aggr.df$N50))))
plot((sort(as.integer(quast_aggr.df$num_contigs))))






###############
############### Abricate
###############

aggr_abricate <- BascetAggregateAbricate(
    "/data/henriksson/github/zorn/test_aggr/abricate/2", 
    #cacheFile=NULL, #option
)

#aggr_abricate
#colnames(aggr_abricate)
#rownames(aggr_abricate)

##include_cells = c("_G1_H6_B7_A10","_H1_D4_H8_H12","_C2_G4_B9_E12")


###############
############### FASTQC
###############

bascetRoot <- "/home/mahogny/github/zorn/test_aggr/fastqc"

BascetCacheComputation(
  bascetRoot,
  "cache_fastqc",
  BascetAggregateFASTQC(
    bascetRoot, verbose=TRUE
  )
)


