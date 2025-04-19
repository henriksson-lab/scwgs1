 


################################################################################
################## Postprocessing with Bascet/Zorn #############################
################################################################################

bascet_runner <- LocalRunner(direct = TRUE)
bascetRoot = "/husky/henriksson/atrandi/v2_wgs_novaseq1/"
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

saveRDS(adata@meta.data, "/husky/henriksson/atrandi/v2_wgs_novaseq1/metadata.RDS")


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


cnt <- t(ReadBascetCountMatrix(bascetRoot,"chromcount"))  #should not need t!

adata <- CreateSeuratObject(counts = CreateAssayObject(cnt), project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!


sort(colSums(cnt))
sort(rowSums(cnt))

#def need a kneeplot!
#test loading one only


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
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:16, reduction.name = "kraken_umap")  ## depth seems to be less of a problem here
} else {
  
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 17)
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:16, reduction.name = "kraken_umap")
}
DimPlot(object = adata, label = TRUE) + NoLegend()




### TODO label cell by maximum species.


### TODO kneeplot

### TODO doublet removal

### TODO barnyard plot



# v2_wgs_novaseq1 chromcount ready!





################################################################################
################## de novo KMER-analysis #######################################  
################################################################################


kmer_hist <- BascetReadMinhashHistogram(bascetRoot)


### Figure 4b
kmer_hist$rank <- 1:nrow(kmer_hist)
ggplot(kmer_hist[kmer_hist$cnt>2,], aes(rank, cnt)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("Count")+
  xlab("Rank")

picked_kmers <- ChooseInformativeKMERs(kmer_hist)

writeLines(picked_kmers, file.path(bascetRoot,"use_kmers.txt"))


#todo mapcell: check if script can do threads. 
#1. if so, run single-threaded and give all threads to the script
#2. skip merge at end if a single thread

#todo mapcell script: 
# reader responsible for extracting the data needed!!



################################################################################
################## SKESA - post analysis ####################################### 
################################################################################


bascetRoot <- "/husky/henriksson/atrandi/bar"
bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_novaseq1/"

### Assemble all genomes
#system("echo START skesa >> time.txt; echo `date +%s` >> time.txt")
my_job <- BascetMapCell(
  bascetRoot,
  withfunction = "_quast",
  inputName = "skesa",
  outputName = "quast",
  runner=bascet_runner
#  bascet_runnerance=bascet_runnerance
)
#WaitForJob(my_job)
#system("echo END skesa >> time.txt; echo `date +%s` >> time.txt")


#/husky/henriksson/atrandi/bar

######### Aggregate data from previous Map call

quast_aggr <- BascetAggregateMap(
  bascetRoot,
  "quast",
  aggr.quast,
  verbose=TRUE
  #,
#  include_cells = c(cellname_coord$cell[1:500])
)

quast_aggr.df <- MapListAsDataFrame(quast_aggr)


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



if(FALSE){
  #check cell!
  cellname_coord <- BascetCellNames(bascetRoot, "quast")
  
  #Open the file, prep for reading
  #bascetFile <- OpenBascet(bascetRoot, bascetName)
  
}




