

################################################################################
################# SNP analysis -- find s.mutans based on % mapping #############
################################################################################

# get aligned counts
cnt_mutans <- ReadBascetCountMatrix(bascetRoot,"cnt_mutans", verbose=FALSE)

cnt_mutans$obs$tot_reads <- as.double(cnt_mutans$obs$`_unmapped` + cnt_mutans$X)
cnt_mutans$obs$frac_mapped <- as.double(cnt_mutans$X/cnt_mutans$obs$tot_reads)
cnt_mutans$obs$`_index` <- stringr::str_remove(cnt_mutans$obs$`_index`,"BASCET_")
rownames(cnt_mutans$obs) <- cnt_mutans$obs$`_index`

### Knee plot of alignment
df <- data.frame(
  frac_mapped = sort(cnt_mutans$obs$frac_mapped, decreasing = TRUE),
  index = 1:nrow(cnt_mutans$obs)
)
ggplot(df, aes(index, frac_mapped)) + 
  geom_line() + 
  scale_x_log10() + 
  theme_bw()

list_cell_mutans <- cnt_mutans$obs$`_index`[which(cnt_mutans$obs$frac_mapped>0.8)]


################################################################################
################# SNP umap - s. mutans #########################################
################################################################################

snp_ad <- ReadCellSNPmatrix(
  "/husky/henriksson/atrandi/v4_wgs_saliva1",
  "cellsnp"
)
#snp_ad <- allmat

snp_adata <- CreateSeuratObject(counts = t(snp_ad), project = "snp_adata3k", min.cells = 3, min.features = 0)
snp_adata

#only keep high quality cells
snp_adata <- snp_adata[,colnames(snp_adata) %in% list_cell_mutans]

snp_adata$frac_mapped <-  cnt_mutans$obs[colnames(snp_adata),]$frac_mapped
snp_adata$tot_reads <-  cnt_mutans$obs[colnames(snp_adata),]$tot_reads

#transfer metadata
snp_meta <- adata@meta.data[colnames(snp_adata),]
snp_adata$genus <- snp_meta$genus
snp_adata$species <- snp_meta$species
snp_adata$phylum <- snp_meta$phylum

### Normalize data, PCA etc
snp_adata <- NormalizeData(snp_adata, normalization.method = "LogNormalize", scale.factor = 10000)
snp_adata <- NormalizeData(snp_adata)
snp_adata <- FindVariableFeatures(snp_adata, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(snp_adata)
snp_adata <- ScaleData(snp_adata, features = all.genes) ########## should normalize by sequencing depth of original cell! not SNP count

snp_adata <- RunPCA(snp_adata, features = VariableFeatures(object = snp_adata))

snp_adata <- FindNeighbors(snp_adata, dims = 1:10)
snp_adata <- FindClusters(snp_adata, resolution = 0.5)

snp_adata <- RunUMAP(snp_adata, dims = 1:10)

snp_adata$snp_count <- rowSums(snp_adata@assays$RNA@layers$counts)

DimPlot(snp_adata, reduction = "umap")

DimPlot(snp_adata, reduction = "umap", group.by = "species")
#DimPlot(snp_adata, reduction = "umap", group.by = "genus")

FeaturePlot(snp_adata, reduction = "umap", features = c("snp_count"))

FeaturePlot(snp_adata, reduction = "umap", features = c("frac_mapped"))
FeaturePlot(snp_adata, reduction = "umap", features = c("tot_reads"))

snp_adata.markers <- FindAllMarkers(snp_adata, only.pos = TRUE)

snp_adata.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  slice_head(n = 100) %>%
  ungroup() -> top10

as.data.frame(top10)

sum(snp_adata.markers$p_val_adj<0.05)

snp_adata.markers[snp_adata.markers$p_val_adj<0.05,]

FeaturePlot(snp_adata, reduction = "umap", features = c("NZ-CP077404.1-7148-T-to-C"))
FeaturePlot(snp_adata, reduction = "umap", features = c("NZ-CP077404.1-1444618-C-to-T"))

FeaturePlot(snp_adata, reduction = "umap", features = c("NZ-CP077404.1-1634018-A-to-G"))


################################################################################
#### SNP - s. mutans - Look at differences in bases to find true strain differences ######
################################################################################

long_basepos <- data.frame(
  pos=paste(str_split_i(rownames(snp_adata),"-",2), str_split_i(rownames(snp_adata),"-",3)),
  #pos=paste(str_split_i(rownames(snp_adata),"-",1), str_split_i(rownames(snp_adata),"-",2), str_split_i(rownames(snp_adata),"-",3)),
  base=str_split_i(rownames(snp_adata),"-",6),
  cnt=rowSums(snp_adata@assays$RNA$counts>0)
)
long_basepos

long_basepos <- sqldf::sqldf("select pos, base, sum(cnt) as cnt from long_basepos group by pos, base")

basepos <- reshape2::acast(long_basepos, pos~base, value.var = "cnt", fill=0)

basepos[rowSums(basepos>1)>1,,drop=FALSE]
#                   A C G N T
#CP077404.1 1825618 0 2 0 0 3

rownames(snp_adata)[str_detect(rownames(snp_adata),"1825618")]

FeaturePlot(snp_adata, reduction = "umap", features = c("NZ-CP077404.1-1825618-C-to-T"))
FeaturePlot(snp_adata, reduction = "umap", features = c("NZ-CP077404.1-1825618-T-to-C"))




