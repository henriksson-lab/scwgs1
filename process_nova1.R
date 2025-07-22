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

library(future)
plan("multicore", workers = 10)


mapSeq2strain <- read.csv("~/github/scwgs/mapSeq2strain.csv")
strain_genomesize <- sqldf::sqldf("select sum(len) as len, strain from mapSeq2strain group by strain")



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

setwd("/home/mahogny/github/zorn") #to put SQL in the right place

################################################################################
################## Kraken-based analysis #######################################  
################################################################################

mat <- ReadBascetCountMatrix(bascetRoot,"kraken", verbose=FALSE)
#mat <- ReadBascetCountMatrix("/home/mahogny/test","kraken_count", verbose=FALSE)
## we really should rename taxids in bascet!

#mat <- ReadBascetKrakenMatrix("/home/mahogny/test", "kraken_count")  #rename to hdf5?

### Read matrix, rename 
#mat <- ReadBascetKrakenMatrix(bascetRoot, "kraken")  #rename to hdf5?
rownames(mat) <- stringr::str_remove(rownames(mat),stringr::fixed("BASCET_"))


taxid_ob <- CreateAssayObject(t(mat))  ## replace any BASCET_? or keep? seurat: ('_'), replacing with dashes ('-')
adata <- CreateSeuratObject(counts = taxid_ob, project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!

# 
# if(TRUE){
#   ### Compress the representation to avoid trouble with some tools; one way of doing it
#   #compressed_mat <- SetTaxonomyNamesFeatures(mat)   #### or add name to the end of taxid_asdasd_XXX
#   
#   head(sort(colSums(compressed_mat), decreasing = TRUE), n=30)
#   
#   ############ This way we can 
#   #note, t() needed from anndata
#   #taxid_ob <- CreateAssayObject(t(compressed_mat))  ## replace any BASCET_? or keep? seurat: ('_'), replacing with dashes ('-')
# } else {
#   ### Compress the representation to avoid trouble with some tools; one way of doing it
#   #compressed_mat <- SetTaxonomyNamesFeatures(mat) 
#   
#   ############ This way we can 
#   #note, t() needed from anndata
#   #taxid_ob <- CreateAssayObject(t(compressed_mat))  ## replace any BASCET_? or keep? seurat: ('_'), replacing with dashes ('-')
#   #adata <- CreateSeuratObject(counts = taxid_ob, project = "proj", min.cells = 0, min.features = 0) ### do we need this? can overload also on assayobject!
#   
# }


# map_cellid_depth <- data.frame(
#   row.names = colnames(adata),
#   cellid=colnames(adata),
#   depth=adata$nCount_RNA
# )


## Add KRAKEN consensus taxonomy to metadata
kraken_taxid <- KrakenFindConsensusTaxonomy(mat)
rownames(kraken_taxid) <- kraken_taxid$cell_id
kraken_taxid <- kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]) #AddMetaData behaves weirdly!!


#### How many species?
##################KrakenSpeciesDistribution(adata)  ### Could use metadata column! TODO   rewrite

saveRDS(kraken_taxid, file.path(bascetRoot,"kraken_metadata.RDS"))


## Dimensional reduction using kraken
#DefaultAssay(adata) <- "kraken"
sum(adata$nCount_RNA>min_nCount_RNA)
adata <- adata[,adata$nCount_RNA>min_nCount_RNA] ## Reduce to sensible number  ; we reduced earlier. I think we should not
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


if(FALSE){
  adata$index_taxid_max <- paste("t",apply(adata@assays$RNA@counts,2,which.max))
  DimPlot(object = adata, label = TRUE, group.by = "index_taxid_max", reduction = "kraken_umap") + NoLegend()
}

############# Barplot of species abundance according to KRAKEN

DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap") + NoLegend()
ggsave(file.path(plotDir, "kraken_umap_genus.pdf"), width=15, height=5)

DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap") + NoLegend()
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
################## Alignment-based analysis ####################################    basics
################################################################################


cnt <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)
adata <- CreateSeuratObject(
  counts = CreateAssayObject(t(cnt)), ## according to anndata standard, there should be a transpose here, unless we abstract this step 
  project = "proj", min.cells = 0, min.features = 0, assay = "chrom_cnt"
) 

## Dimensional reduction using kraken
#DefaultAssay(adata) <- "kraken"
sum(adata$nCount_chrom_cnt>min_nCount_RNA)
adata <- adata[,adata$nCount_chrom_cnt>min_nCount_RNA] ## Reduce to sensible number
adata


########## Produce a count matrix on strain level
adata[["species_cnt"]] <- ChromToSpeciesCount(adata, mapSeq2strain)  #gives warning. coerce ourselves to dgCMatrix

#Figure out which species has most reads in which cell
cnt <- adata@assays$species_cnt$counts
adata$species_aln <- rownames(cnt)[apply(cnt, 2, which.max)]

adata$log_cnt <- log10(1+adata$nCount_chrom_cnt)

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
DimPlot(object = adata, label = TRUE, group.by = "species_aln") + 
  xlab("") + 
  ylab("") #+ NoLegend()
ggsave(file.path(plotDir,"umap_alignment.pdf"), width = 5, height = 5)


### TODO kneeplot
### TODO doublet removal

##########
########## Barnyard plot
##########

#BarnyardPlotMatrix()  this is an option

bp <- data.frame(
  maxc = MatrixGenerics::colMaxs(adata@assays$species_cnt$counts),
  totc = colSums(adata@assays$species_cnt$counts)
)
bp$restc <- bp$totc - bp$maxc

ggplot(bp, aes(maxc,restc)) +   #,color=log10(totc)
  geom_point() +
  theme_bw() +
  xlab("Dominant species count") +
  ylab("Other species count")
ggsave(file.path(plotDir,"barnyard_alignment.pdf"), width = 5, height = 5)



################################################################################
################## informative KMER-analysis ################################### plotting of distribution
################################################################################

p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")
if(file.exists(p_minh)) {
  kmerHist <- BascetReadMinhashHistogram(bascetRoot)
  #kmerHist <- kmerHist[order(kmerHist$cnt, decreasing=TRUE),]
  if(!file.exists(p_use_kmers)) {
    picked_kmers <- ChooseInformativeKMERs( ### 0.002 => 11212 this killed conda!
      kmerHist,
      minFreq = 0.002
      )
    writeLines(picked_kmers, p_use_kmers)
  }
}

#kmerHist <- BascetReadMinhashHistogram(bascetRoot)


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



################################################################################
################## informative KMER-analysis ################################### umap etc
################################################################################



### Read count matrix ///////////////// 
cnt <- ReadBascetCountMatrix(  #  x[.,.] <- val : x being coerced from Tsparse* to CsparseMatrix  ---- do this conversion manually
  bascetRoot,
  "kmer_counts"
)  #### TODO: should read and concatenate multiple matrices; different name?
dim(cnt)
#rownames(cnt) <- paste0("_",rownames(cnt))  #### todo fix naming


if(FALSE){
  colSums(cnt)
  colSums(cnt>0)
  
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
  
  #colnames(cnt) <- paste0("BASCET_",colnames(cnt)) ### compatibilíty with fragments
}



sum(rowSums(cnt)>20000)
adata <- CreateSeuratObject(
  counts = CreateAssayObject(t(cnt[rowSums(cnt)>20000,])),  #cutoff on kmer count
#  counts = CreateAssayObject(cnt),
  assay = "infokmer"
)

if(FALSE){
  sum(rowSums(cnt)>20000)
  adata <- CreateSeuratObject(
    counts = CreateAssayObject(t(cnt[rowSums(cnt)>20000,]>0)),
    #  counts = CreateAssayObject(cnt),
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
    #  counts = CreateAssayObject(cnt),
    assay = "infokmer"
  )
  adata  
}


#Add kraken metadata
kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]) #AddMetaData behaves weirdly!!


##### Perform ATAC-seq style analysis
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata <- RunSVD(adata)

DepthCor(adata)
ggsave(file.path(plotDir,"info_kmer_depthcor.pdf"), width = 5, height = 5)

adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "infokmers_umap")  ## dim 1 affected plenty if binarizing matrix

### Plots
DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "infokmers_umap") #+ NoLegend()
FeaturePlot(adata, "nCount_infokmer", reduction = "infokmers_umap")

DimPlot(object = adata[,adata$genus=="Bacillus"], label = TRUE, group.by = "genus", reduction = "infokmers_umap") + NoLegend()

DimPlot(object = adata, label = TRUE, reduction = "infokmers_umap") + NoLegend()




table(adata$genus) #cereibacter dominates plenty!!




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



#TODO: which cell has the most unclassified reads? id=0. does bascet care?
#bascet discards such reads ... store in some matrix too?












################################################################################
################## Count sketch-based analysis #################################
################################################################################

#Kneeplot not possible due to cutoff
if(FALSE){
  df <- data.frame(
    depth=sort(decreasing = TRUE,adata$celldepth)
  )
  df$index <- 1:nrow(df)
  ggplot(df, aes(index, depth)) + scale_x_log10() + scale_y_log10() + geom_line()
}

#Load data as seurat object 
adata <- BascetLoadCountSketchMatrix(bascetRoot)

#Subset cells
keep_cells <- adata$celldepth>200000
sum(keep_cells)
adata <- adata[,keep_cells]


#Add kraken metadata
kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
adata@meta.data <- cbind(adata@meta.data,kraken_taxid[colnames(adata),c("taxid","phylum","class","order","family","genus","species")]) #AddMetaData behaves weirdly!!
#adata$genus <- kraken_meta[colnames(adata),]$genus

#Non-linear dimensional reduction
adata <- RunUMAP(adata, dims = 1:ncol(adata@reductions$kmersketch@cell.embeddings), reduction = "kmersketch")  #Searching Annoy index using 1 thread, search_k = 3000 ; can do more
DimPlot(object = adata, label = TRUE) + NoLegend()

#Plotting
DimPlot(object = adata, label = TRUE, group.by = "genus") + NoLegend()
DimPlot(object = adata, label = TRUE, group.by = "species") + NoLegend()
FeaturePlot(adata, features = "depth")


if(FALSE){
  
  ## For saliva
  
  table_species <- sort(table(adata@meta.data$species))
  table_genus <- sort(table(adata@meta.data$genus))
  
  keep_genus <- setdiff(names(table_genus)[table_genus>10],"Xanthomonas")
  keep_genus <- setdiff(names(table_genus)[table_genus>100],c("Xanthomonas","Streptococcus"))
  
  
  DimPlot(object = adata[,adata$genus %in% keep_genus], label = TRUE, group.by = "genus") #+ NoLegend()
  
}


