

################################################################################
################# SNP analysis -- find s.mutans based on % mapping #############
################################################################################

# get aligned counts
cnt_snp <- ReadBascetCountMatrix(bascetRoot,"cnt_bestmatch", verbose=FALSE)

dim(cnt_snp$obs)
dim(cnt_snp$X)

cnt_snp$obs$mapped <- rowSums(cnt_snp$X)
cnt_snp$obs$tot_reads <- as.double(cnt_snp$obs$`_unmapped` + cnt_snp$obs$mapped)
cnt_snp$obs$frac_mapped <- as.double(cnt_snp$obs$mapped/cnt_snp$obs$tot_reads)
cnt_snp$obs$`_index` <- stringr::str_remove(cnt_snp$obs$`_index`,"BASCET_")
rownames(cnt_snp$obs) <- cnt_snp$obs$`_index`

### Knee plot of alignment
df <- data.frame(
  frac_mapped = sort(cnt_snp$obs$frac_mapped, decreasing = TRUE),
  index = 1:nrow(cnt_snp$obs)
)
ggplot(df, aes(index, frac_mapped)) + 
  geom_line() + 
  scale_x_log10() + 
  theme_bw()

list_cell_bestmatch <- cnt_snp$obs$`_index`[which(cnt_snp$obs$frac_mapped>0.75)]


################################################################################
################# SNP umap - bestmatch #########################################
################################################################################

snp_ad <- ReadCellSNPmatrix(
  "/husky/henriksson/atrandi/v4_wgs_saliva1",
  "cellsnp_bestmatch",
  listCell = list_cell_bestmatch
)
#snp_ad <- allmat

snp_adata_bestmatch <- CreateSeuratObject(counts = t(snp_ad), project = "snp_adata_bestmatch3k", min.cells = 3, min.features = 0)
snp_adata_bestmatch

#only keep high quality cells
snp_adata_bestmatch <- snp_adata_bestmatch[,colnames(snp_adata_bestmatch) %in% list_cell_bestmatch]

snp_adata_bestmatch$frac_mapped <-  cnt_snp$obs[colnames(snp_adata_bestmatch),]$frac_mapped
snp_adata_bestmatch$tot_reads <-  cnt_snp$obs[colnames(snp_adata_bestmatch),]$tot_reads

#transfer metadata
snp_meta <- adata@meta.data[colnames(snp_adata_bestmatch),]
snp_adata_bestmatch$genus <- snp_meta$genus
snp_adata_bestmatch$species <- snp_meta$species
snp_adata_bestmatch$phylum <- snp_meta$phylum

### Normalize data, PCA etc
snp_adata_bestmatch <- NormalizeData(snp_adata_bestmatch, normalization.method = "LogNormalize", scale.factor = 10000)
snp_adata_bestmatch <- NormalizeData(snp_adata_bestmatch)
snp_adata_bestmatch <- FindVariableFeatures(snp_adata_bestmatch, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(snp_adata_bestmatch)
snp_adata_bestmatch <- ScaleData(snp_adata_bestmatch, features = all.genes) ########## should normalize by sequencing depth of original cell! not SNP count

snp_adata_bestmatch <- RunPCA(snp_adata_bestmatch, features = VariableFeatures(object = snp_adata_bestmatch))

snp_adata_bestmatch <- FindNeighbors(snp_adata_bestmatch, dims = 1:10)
snp_adata_bestmatch <- FindClusters(snp_adata_bestmatch, resolution = 0.2)

snp_adata_bestmatch <- RunUMAP(snp_adata_bestmatch, dims = 1:10)

snp_adata_bestmatch$snp_count <- rowSums(snp_adata_bestmatch@assays$RNA@layers$counts)

DimPlot(snp_adata_bestmatch, reduction = "umap")

DimPlot(snp_adata_bestmatch, reduction = "umap", group.by = "species")
#DimPlot(snp_adata_bestmatch, reduction = "umap", group.by = "genus")

FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("snp_count"))

FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("frac_mapped"))
FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("tot_reads"))

snp_adata_bestmatch.markers <- FindAllMarkers(snp_adata_bestmatch, only.pos = TRUE)

snp_adata_bestmatch.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  slice_head(n = 100) %>%
  ungroup() -> top10

as.data.frame(top10)

sum(snp_adata_bestmatch.markers$p_val_adj<0.05)

snp_adata_bestmatch.markers[snp_adata_bestmatch.markers$p_val_adj<0.05,]

FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("Contig-1649-82.271-323-A-to-G"))
#FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("NZ-CP077404.1-1444618-C-to-T"))
#FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("NZ-CP077404.1-1634018-A-to-G"))


################################################################################
#### SNP - bestmatch - Look at differences in bases to find true strain differences ######
################################################################################

long_basepos <- data.frame(
  pos=paste(str_split_i(rownames(snp_adata_bestmatch),"-",2), str_split_i(rownames(snp_adata_bestmatch),"-",3)),
  #pos=paste(str_split_i(rownames(snp_adata_bestmatch),"-",1), str_split_i(rownames(snp_adata_bestmatch),"-",2), str_split_i(rownames(snp_adata_bestmatch),"-",3)),
  base=str_split_i(rownames(snp_adata_bestmatch),"-",7),
  cnt=rowSums(snp_adata_bestmatch@assays$RNA$counts>0)
)
long_basepos

long_basepos <- sqldf::sqldf("select pos, base, sum(cnt) as cnt from long_basepos group by pos, base")
basepos <- reshape2::acast(long_basepos, pos~base, value.var = "cnt", fill=0)

basepos[rowSums(basepos>1)>1,,drop=FALSE]

#rownames(snp_adata_bestmatch)[str_detect(rownames(snp_adata_bestmatch),"1825618")]
#FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("NZ-CP077404.1-1825618-C-to-T"))
#FeaturePlot(snp_adata_bestmatch, reduction = "umap", features = c("NZ-CP077404.1-1825618-T-to-C"))




