
################################################################################
################# SNP v1
################################################################################

if(FALSE){
  
  adata <- readRDS("/husky/henriksson/atrandi/v4_wgs_saliva1/cache_adata_kraken.RDS")
  #DefaultAssay(adata1) <- "species_cnt"
  colnames(adata)
  sort(unique(adata$phylum))
  sort(unique(adata$genus))
  adata$genus=="Streptococcus"
  
  
  ####################
  ####################
  ####################
  
  # # /pfs/proj/nobackup/fs/projnb10/carroll_hpc2n/hadrien/sc_strep/cellsnp_results/cellSNP.cells.vcf
  # 
  # # https://cellsnp-lite.readthedocs.io/en/latest/main/manual.html#mode-1-pileup-with-given-snps
  # 
  # #cellSNP.base.vcf  one line per feature
  # snp_info <- read.table("/husky/henriksson/atrandi/v4_wgs_saliva1/cellsnp_results/cellSNP.base.vcf")
  # colnames(snp_info) <- c("chrom","pos","id","ref","alt","qual","filter","info")
  # rownames(snp_info) <- sprintf("Feature%s",1:nrow(snp_info))
  # 
  # snp_cellid <- readLines("/husky/henriksson/atrandi/v4_wgs_saliva1/cellsnp_results/cellSNP.samples.tsv")
  # length(snp_cellid)
  # 
  # #depth of ALT alleles
  # snp_ad <- Matrix::readMM("/husky/henriksson/atrandi/v4_wgs_saliva1/cellsnp_results/cellSNP.tag.AD.mtx")
  # colnames(snp_ad) <- snp_cellid
  # dim(snp_ad)
  # 
  # #depth of ALT+REF alleles
  # snp_dp <- Matrix::readMM("/husky/henriksson/atrandi/v4_wgs_saliva1/cellsnp_results/cellSNP.tag.DP.mtx")
  # colnames(snp_dp) <- snp_cellid
  # dim(snp_dp)
  # 
  # #snp_ratio <- snp_ad/snp_dp  #this is huge! not a good way? divisions by 0
  
  
  pbmc <- CreateSeuratObject(counts = snp_ad, project = "pbmc3k", min.cells = 3, min.features = 0)
  pbmc
  
  #only keep high quality cells
  pbmc <- pbmc[,colnames(pbmc) %in% colnames(adata)]#@meta.data
  
  #rownames(pbmc)
  #pbmc$genus
  
  #transfer metadata
  snp_meta <- adata@meta.data[colnames(pbmc),]
  pbmc$genus <- snp_meta$genus
  pbmc$species <- snp_meta$species
  pbmc$phylum <- snp_meta$phylum
  
  pbmc <- pbmc[,which(pbmc$genus=="Streptococcus")] ### only one species
  
  #pbmc$genus=="Streptococcus"
  #unique(pbmc$genus)
  
  ### Normalize data, PCA etc
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes) ########## should normalize by sequencing depth of original cell! not SNP count
  
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  #DimPlot(pbmc, reduction = "pca") + NoLegend()
  #ElbowPlot(pbmc)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  
  DimPlot(pbmc, reduction = "umap")
  
  DimPlot(pbmc, reduction = "umap", group.by = "species")
  DimPlot(pbmc, reduction = "umap", group.by = "genus")
  
  #DimPlot(pbmc[,pbmc$genus=="Streptococcus"], reduction = "umap", group.by = "genus")
  
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers$feature <- pbmc.markers$gene
  pbmc.markers_m <- merge(snp_info, pbmc.markers)
  
  pbmc.markers_m %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 2) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  as.data.frame(top10)
  #head(snp_info)
  
  snp_info$feature <- rownames(snp_info)
  
  snp_info["Feature508",]
  
  FeaturePlot(pbmc, reduction = "umap", features = c("Feature1052","Feature1055"))
  FeaturePlot(pbmc, reduction = "umap", features = c("Feature10092","Feature10095","Feature10096","Feature1002"))
  
  FeaturePlot(pbmc, reduction = "umap", features = c("Feature1052","Feature1055"))
  FeaturePlot(pbmc, reduction = "umap", features = c("Feature10092","Feature10095","Feature10096","Feature1002"))
  
  FeaturePlot(pbmc, reduction = "umap", features = c("Feature1052","Feature1002","Feature10096"))
  #FeaturePlot(pbmc, reduction = "umap", features = c("Feature10092","Feature10095","Feature10096","Feature1002"))
  
  
  FeaturePlot(pbmc, reduction = "umap", features = c("Feature1052")) + 
    ggtitle(paste(snp_info["Feature1052",]$pos,snp_info["Feature1052",]$ref,">",snp_info["Feature1052",]$alt))
  
  plot_one_snp <- function(snp_name){
    FeaturePlot(pbmc, reduction = "umap", features = c(snp_name)) + 
      ggtitle(paste(snp_info[snp_name,]$pos,snp_info[snp_name,]$ref,">",snp_info[snp_name,]$alt))
  }
  
  egg::ggarrange(
    plot_one_snp("Feature1002"),
    plot_one_snp("Feature1052"),
    plot_one_snp("Feature10096"),
    nrow=1
  )
  
  
  
  if(str_detect(bascetRoot,"saliva")){
    
    ##########
    ########## Read the aligned counts
    ##########
    
    cnt <- ReadBascetCountMatrix(bascetRoot,"cnt_mutans", verbose=FALSE)
    adata <- CreateSeuratObject(
      counts = CreateAssayObject(t(cnt)), ## according to anndata standard, there should be a transpose here, unless we abstract this step 
      project = "proj", min.cells = 0, min.features = 0, assay = "chrom_cnt"
    ) 
  }
  
  
}



################################################################################
################# SNP -- which cells have well-aligned reads? ##################
################################################################################

# get aligned counts


cnt_mutans <- ReadBascetCountMatrix(bascetRoot,"cnt_mutans", verbose=FALSE)

cnt_mutans$obs$tot_reads <- as.double(cnt_mutans$obs$`_unmapped` + cnt_mutans$X)
cnt_mutans$obs$frac_mapped <- as.double(cnt_mutans$X/cnt_mutans$obs$tot_reads)
cnt_mutans$obs$`_index` <- stringr::str_remove(cnt_mutans$obs$`_index`,"BASCET_")
rownames(cnt_mutans$obs) <- cnt_mutans$obs$`_index`

#hist(as.double(cnt_mutans$obs$frac_mapped))

df <- data.frame(
  frac_mapped = sort(cnt_mutans$obs$frac_mapped, decreasing = TRUE),
  index = 1:nrow(cnt_mutans$obs)
)
ggplot(df, aes(index, frac_mapped)) + 
  geom_line() + 
  scale_x_log10() + 
  theme_bw()

list_cell_mutans <- cnt_mutans$obs$`_index`[which(cnt_mutans$obs$frac_mapped>0.8)]
#list_cell_mutans <- stringr::str_remove(list_cell_mutans,"BASCET_")


################################################################################
################# SNP v2
################################################################################

#basedir <- "/husky/henriksson/atrandi/v4_wgs_saliva1/cellsnp_results"

snp_ad <- ReadCellSNPmatrix(
  "/husky/henriksson/atrandi/v4_wgs_saliva1",
  "cellsnp"
)
  
snp_ad <- allmat

pbmc <- CreateSeuratObject(counts = t(snp_ad), project = "pbmc3k", min.cells = 3, min.features = 0)
pbmc

#only keep high quality cells
#pbmc <- pbmc[,colnames(pbmc) %in% colnames(adata)]
pbmc <- pbmc[,colnames(pbmc) %in% list_cell_mutans]

pbmc$frac_mapped <-  cnt_mutans$obs[colnames(pbmc),]$frac_mapped
pbmc$tot_reads <-  cnt_mutans$obs[colnames(pbmc),]$tot_reads
#adata$_unmapped <-  cnt_mutans$obs[colnames(adata),]$_unmapped


#rownames(pbmc)
#pbmc$genus

#transfer metadata
snp_meta <- adata@meta.data[colnames(pbmc),]
pbmc$genus <- snp_meta$genus
pbmc$species <- snp_meta$species
pbmc$phylum <- snp_meta$phylum

#pbmc <- pbmc[,which(pbmc$genus=="Streptococcus")] ### only one species


### Normalize data, PCA etc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) ########## should normalize by sequencing depth of original cell! not SNP count

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

pbmc$snp_count <- rowSums(pbmc@assays$RNA@layers$counts)

DimPlot(pbmc, reduction = "umap")

DimPlot(pbmc, reduction = "umap", group.by = "species")
#DimPlot(pbmc, reduction = "umap", group.by = "genus")

FeaturePlot(pbmc, reduction = "umap", features = c("snp_count"))

FeaturePlot(pbmc, reduction = "umap", features = c("frac_mapped"))
FeaturePlot(pbmc, reduction = "umap", features = c("tot_reads"))

#DimPlot(pbmc[,pbmc$genus=="Streptococcus"], reduction = "umap", group.by = "genus")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  slice_head(n = 100) %>%
  ungroup() -> top10

as.data.frame(top10)

sum(pbmc.markers$p_val_adj<0.05)

pbmc.markers[pbmc.markers$p_val_adj<0.05,]

FeaturePlot(pbmc, reduction = "umap", features = c("NZ-CP077404.1-7148-T-to-C"))
FeaturePlot(pbmc, reduction = "umap", features = c("NZ-CP077404.1-1444618-C-to-T"))

FeaturePlot(pbmc, reduction = "umap", features = c("NZ-CP077404.1-1634018-A-to-G"))




