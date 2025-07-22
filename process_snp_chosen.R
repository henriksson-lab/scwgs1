

################################################################################
################# SNP analysis -- find chosen based on % mapping ###############
################################################################################

chosen_cell <- "H1_E5_F9_E11"

# get aligned counts
cnt_snp <- ReadBascetCountMatrix(bascetRoot,"cnt_chosen", verbose=FALSE)

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

list_cell_chosen <- cnt_snp$obs$`_index`[which(cnt_snp$obs$frac_mapped>0.75)]


################################################################################
################# SNP umap - chosen ############################################
################################################################################

snp_ad <- ReadCellSNPmatrix(
  "/husky/henriksson/atrandi/v4_wgs_saliva1",
  "cellsnp_chosen",
  listCell = list_cell_chosen
)
#snp_ad <- allmat

snp_adata_chosen <- CreateSeuratObject(counts = t(snp_ad), project = "snp_adata_chosen3k", min.cells = 3, min.features = 0)
snp_adata_chosen

snp_adata_chosen$is_chosen <- colnames(snp_adata_chosen)=="H1_E5_F9_E11"

#only keep high quality cells
snp_adata_chosen <- snp_adata_chosen[,colnames(snp_adata_chosen) %in% list_cell_chosen]

snp_adata_chosen$frac_mapped <-  cnt_snp$obs[colnames(snp_adata_chosen),]$frac_mapped
snp_adata_chosen$tot_reads <-  cnt_snp$obs[colnames(snp_adata_chosen),]$tot_reads

#transfer metadata from kraken adata
snp_meta <- adata@meta.data[colnames(snp_adata_chosen),]
snp_adata_chosen$genus <- snp_meta$genus
snp_adata_chosen$species <- snp_meta$species
snp_adata_chosen$phylum <- snp_meta$phylum


#transfer metadata from aligned counts
snp_adata_chosen$frac_mapped <- cnt_snp$obs[colnames(snp_adata_chosen),]$frac_mapped
snp_adata_chosen$tot_reads <- cnt_snp$obs[colnames(snp_adata_chosen),]$tot_reads


### Normalize data, PCA etc
snp_adata_chosen <- NormalizeData(snp_adata_chosen, normalization.method = "LogNormalize", scale.factor = 10000)
snp_adata_chosen <- NormalizeData(snp_adata_chosen)
snp_adata_chosen <- FindVariableFeatures(snp_adata_chosen, selection.method = "vst", nfeatures = 200)

snp_adata_chosen$snp_count <- as.double(colSums(snp_adata_chosen@assays$RNA@layers$counts))
snp_adata_chosen$log_snp_count <- log10(snp_adata_chosen$snp_count)
snp_adata_chosen$log_tot_reads <- log10(snp_adata_chosen$tot_reads)

#plot1 <- VariableFeaturePlot(snp_adata_chosen)
#plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(snp_adata_chosen), 10), repel = TRUE)
#plot1 + plot2

all.genes <- rownames(snp_adata_chosen)
#snp_adata_chosen <- ScaleData(snp_adata_chosen, features = all.genes) ########## should normalize by sequencing depth of original cell! not SNP count
snp_adata_chosen <- ScaleData(snp_adata_chosen, features = all.genes, vars.to.regress="tot_reads") ########## should normalize by sequencing depth of original cell! not SNP count

snp_adata_chosen <- RunPCA(snp_adata_chosen, features = all.genes)
#snp_adata_chosen <- RunPCA(snp_adata_chosen, features = VariableFeatures(object = snp_adata_chosen))

snp_adata_chosen <- FindNeighbors(snp_adata_chosen, dims = 1:10)
snp_adata_chosen <- FindClusters(snp_adata_chosen, resolution = 1)

snp_adata_chosen <- RunUMAP(snp_adata_chosen, dims = 1:10)

saveRDS(adata, file.path("/husky/henriksson/atrandi/v4_wgs_saliva1","cache_adata_snp.RDS"))
adata <- readRDS(file.path("/husky/henriksson/atrandi/v4_wgs_saliva1","cache_adata_snp.RDS"))


#hist(snp_adata_chosen$snp_count)
#hist(rowSums(snp_adata_chosen@assays$RNA$counts))


DimPlot(snp_adata_chosen, reduction = "umap")

DimPlot(snp_adata_chosen, reduction = "umap", group.by = "species")
#DimPlot(snp_adata_chosen, reduction = "umap", group.by = "genus")

FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("log_snp_count"))
FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("snp_count"))
FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("frac_mapped"))
FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("tot_reads"))
FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("log_tot_reads"))
DimPlot(snp_adata_chosen, reduction = "umap", group.by = "is_chosen")




snp_adata_chosen.markers <- FindAllMarkers(snp_adata_chosen, only.pos = FALSE)
snp_adata_chosen.markers[snp_adata_chosen.markers$p_val_adj<0.05,]

snp_adata_chosen.markers <- snp_adata_chosen.markers[
  snp_adata_chosen@assays$RNA$counts[snp_adata_chosen.markers$gene,"H1_E5_F9_E11"]==0,
]


snp_adata_chosen.markers[snp_adata_chosen.markers$p_val_adj<0.5,] %>%
  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 0.5) %>%
  #slice_head(n = 100) %>%
  ungroup() -> top10

as.data.frame(top10)



#sum(snp_adata_chosen.markers$p_val_adj<0.05)
#snp_adata_chosen.markers[snp_adata_chosen.markers$p_val_adj<0.05,]

FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("Contig-2211-617.016-344-A-to-G"))
#FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("NZ-CP077404.1-1444618-C-to-T"))
#FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("NZ-CP077404.1-1634018-A-to-G"))

FeaturePlot(snp_adata_chosen, reduction = "umap", features = top10$gene)

FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("Contig-2869-663.065-338-T-to-A"))
FeaturePlot(snp_adata_chosen, reduction = "umap", features = c("Contig-2869-663.065-338-T-to-A"))




# in glucose transferase
#3 4.043781e-08   9.393761 0.392 0.000 0.007429276       1 Contig-2589-882.912-2480-C-to-T       no_features; blastn (default megablast on 15 June 2025)
#4 1.007050e-07   9.092323 0.373 0.000 0.018501618       1 Contig-2589-882.912-2515-A-to-G       no_features
#5 1.072092e-07   4.470328 0.412 0.015 0.019696576       1 Contig-2589-882.912-2398-C-to-T       no_features


FeaturePlot(snp_adata_chosen, reduction = "umap", features = c(
  "Contig-2589-882.912-2480-C-to-T",
  "Contig-2589-882.912-2515-A-to-G",
  "Contig-2589-882.912-2398-C-to-T"
))


snp_ad["H1_E5_F9_E11","Contig-2589-882.912-2480-C-to-T"]
snp_adata_chosen@assays$RNA$counts["Contig-2589-882.912-2480-C-to-T",]


############### verify umap
toplot <- data.frame(
  cell=colnames(snp_adata_chosen),
  cnt_alt=snp_adata_chosen@assays$RNA$counts["Contig-432-477.703-51-A-to-G",]>0,
  x=snp_adata_chosen@reductions$umap@cell.embeddings[,1],
  y=snp_adata_chosen@reductions$umap@cell.embeddings[,2],
  tot_reads=snp_adata_chosen$tot_reads
)
toplot
x_ref <- snp_adata_chosen@reductions$umap@cell.embeddings[chosen_cell,1]
y_ref <- snp_adata_chosen@reductions$umap@cell.embeddings[chosen_cell,2]
ggplot(toplot, aes(x,y, color=cnt_alt)) +
  geom_point() + 
  geom_point(x=x_ref, y=y_ref, color="#BB5566",size=3) +
  theme_bw() +
  xlab("SNP1 Contig-432-477.703:51 A > G") + 
  ylab("SNP2") +
  scale_color_manual(values = c("#004488","#ddcc77","#BB5566")) +
  theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank())
ggsave(file.path(plotDirAll, "saliva_selected_snp.svg"), width = 4, height = 3)




if(FALSE){
  ggplot(toplot, aes(x,y, color=log(tot_reads))) +
    geom_point() + 
    geom_point(x=x_ref, y=y_ref, color="red",size=3)+
    #  ggrepel::geom_text_repel(x=x_ref, y=y_ref, label=chosen_cell, max.overlaps = 1000)+
    theme_bw()
}


#toplot[toplot$cell=="H1_E5_F9_E11",]

################### use for paper! focus on 3 SNPs

#> as.data.frame(top10)
#p_val avg_log2FC pct.1 pct.2  p_val_adj cluster                            gene
#1 4.128794e-07   2.250739 0.538 0.036 0.02522404       1    Contig-432-477.703-51-A-to-G
#2 7.029042e-07   2.641722 0.462 0.018 0.04294253       1    Contig-432-477.703-57-A-to-T
#3 9.583142e-07   9.396747 0.385 0.000 0.05854629       1 Contig-2221-507.388-3012-A-to-G

