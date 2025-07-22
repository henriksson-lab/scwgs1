jonas.pal <- c("#15889C", "#ED496F", "#8E1382", 
               "#FFB300", "#43A047", "#FF6F00", 
               "#C62828", "#2D62A3", "#1F7530", 
               "#573794", "#28AAE1", "#870F0F", 
               "#9E9E9E", "#795548", "#CD6155", 
               "#F7DC6F", "#7DCEA0", "#85C1E9", 
               "#EB984E", "#EAD9D5", "#03CD4A", 
               "#CDDC39", "#E0594B", "#C76CDE", 
               "#24B177", "#8D6E63", "#486FF7", 
               "#6300B5", "#88E200", "#012824", 
               "#0D3290", "#A347FB", "#54FC7A", 
               "#EB1388", "#B0978D", "#FE52CF", 
               "#83F1F6", "#F1F847", "#2B1DFC", 
               "#6C6F15", "#6CA05C", "#7788CD", 
               "#F502F3", "#0DC290", "#FA0E03", 
               "#3CAA0A", "#BEFC8D", "#08F8EB", 
               "#B1CD3F", "#D6A5FA", "#CE606C", 
               "#AB1EBA", "#6ECC9F", "#054DDC", 
               "#486FF7", "#854F49", "#F22B21", 
               "#3A0E43", "#225805", "#37D160", 
               "#E4B974", "#A8BADE", "#47EDD1", 
               "#F47A92", "#C76CDE", "#9106EB", 
               "#81AA20", "#D7FDFD", "#5DEB2E", 
               "#F82745", "#6435E0", "#027FFE", 
               "#8E3101", "#16F648", "#1C15BC",
               "#8BE46E", "#8D6FA0", "#E68FC6", 
               "#058CA9", "#9E018A", "#BDFD0B", 
               "#B22760", "#2BF49F", "#CB9348", 
               "#9D8303", "#C251A1", "#46ADAF", 
               "#A3E3AF", "#22BB34", "#6EA3FA", 
               "#260374", "#1C3854", "#405D37", 
               "#C21DF3", "#FCEA92", "#537F88", 
               "#FD4C18", "#F2D71E", "#FD4C7A")



dataset_name <- "v4_wgs_saliva1"
bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)


adata1 <- readRDS("/husky/henriksson/atrandi/v4_wgs_saliva1/cache_adata_kraken.RDS")
adata2 <- readRDS("/husky/henriksson/atrandi/v4_wgs_saliva1/cache_adata_infokmer.RDS")
adata3 <- readRDS("/husky/henriksson/atrandi/v4_wgs_saliva1/cache_adata_cs.RDS")

adata1 <- adata1[,adata1$phylum!="NA"]
adata2 <- adata2[,adata2$phylum!="NA"]
adata3 <- adata3[,adata3$phylum!="NA"]

# stats
adata1
#10k cells were kept
sum(adata1$genus=="Streptococcus", na.rm = TRUE) #6034
sum(adata1$species=="Streptococcus salivarius", na.rm = TRUE)


################################################################################
########################## Fig 4 umap matrix ###################################
################################################################################

#plotlevel <- "genus"
plotlevel <- "phylum"

b_p1 <- DimPlot(object = adata1, group.by = plotlevel, reduction = "kraken_umap") + 
  xlab("KRAKEN1") + ylab("KRAKEN2") + 
  scale_color_manual(values = jonas.pal)
b_p2 <- DimPlot(object = adata2, group.by = plotlevel, reduction = "infokmers_umap") + 
  xlab("IK1") + ylab("IK2") + 
  scale_color_manual(values = jonas.pal)
b_p3 <- DimPlot(object = adata3, group.by = plotlevel) + 
  xlab("CS1") + ylab("CS2") + 
  scale_color_manual(values = jonas.pal)

ptot <- egg::ggarrange(
  b_p1+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  b_p2+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  b_p3+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),

  nrow=1)
ptot
ggsave(plot = ptot, file.path(plotDirAll, "fig4_umaps.svg"), width = 7, height = 2.5)
ggsave(plot = b_p1, file.path(plotDirAll, "fig4_one_umap.svg"), width = 6, height = 4)





################################################################################
########################## Fig 4 kneeplot all ###################################
################################################################################

### filtered too hard?

cnt <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)

df <- data.frame(
  cnt=sort(cnt$obs$`_unmapped` + rowSums(cnt$X), decreasing = TRUE)
#  cnt=sort(colSums(adata@assays$RNA$counts), decreasing = TRUE)
)
df$index <- 1:nrow(df) 
ggplot(df, aes(index, cnt)) + c
  geom_line() + 
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + 
  ylab("Read count") + 
  xlab("Cell index")

#ggsave(file.path(plotDir,"alignment_kneeplot_all.pdf"), width = 4, height = 3)
ggsave(file.path(plotDirAll, "fig4_all_kneeplots.svg"), width=3, height=3, limitsize=FALSE)


################################################################################
########################## Fig 4 kneeplot per phylum ###########################
################################################################################


plotlevel <- "phylum"


adata <- readRDS("/husky/henriksson/atrandi/v4_wgs_saliva1/cache_adata_kraken.RDS")

show_num_spec <- 100
KrakenKneePlot(adata, groupby = "phylum", show_num_spec=show_num_spec, sortbyname=TRUE) + 
  scale_color_manual(values = jonas.pal)
ggsave(file.path(plotDirAll, "fig4_strain_kneeplots.svg"), width=5, height=3, limitsize=FALSE)


#ggsave(plot = ptot, file.path(plotDirAll, "fig4_strain_kneeplots.svg"), width = 5, height = 4)








################################################################################
######################### AmrFinder for saliva #################################
################################################################################


bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_saliva1/"
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))

amrfinder <- read.csv("/home/mahogny/github/scwgs/saliva_contigs_amrfinder/amrfinder_matrix.tsv",sep="\t")
rownames(amrfinder) <- amrfinder$cell

amrfinder <- amrfinder[colnames(adata),]
amrfinder[is.na(amrfinder)] <- 0
amrfinder <- amrfinder[,-1]

adata@meta.data <- cbind(adata@meta.data, amrfinder)

FeaturePlot(adata, features = colnames(amrfinder), ncol=3)
#ggsave(file.path(plotDirAll, "saliva_amrfinder.png"), width = 10, height = 12)



############ Nicer looking UMAPs

list_plots <- list()
for(cur_feature in colnames(amrfinder)){
  
  toplot <- data.frame(
    cnt=adata@meta.data[,cur_feature]>0,
    x=adata@reductions$kraken_umap@cell.embeddings[,1],
    y=adata@reductions$kraken_umap@cell.embeddings[,2]
  )
  toplot <- toplot[order(toplot$cnt, decreasing = FALSE),]
  
  list_plots[[cur_feature]] <- ggplot(toplot, aes(x,y, color=cnt)) +
    geom_point() + scale_color_manual(values=c("gray","red"), guide="none")+    
    theme(legend.position = "none")+
    theme_bw() + 
    xlab(cur_feature) + ylab("")
  
}
ptot <- egg::ggarrange(plots=list_plots, ncol=3)
ggsave(plot = ptot, file.path(plotDirAll, "saliva_amrfinder.png"), width = 10, height = 12)











################################################################################
######################### GECCO for saliva #####################################
################################################################################


bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_saliva1/"
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))

gecco <- read.csv("/home/mahogny/github/scwgs/saliva_contigs_bgc/bgc_matrix.tsv",sep="\t")
rownames(gecco) <- gecco$cell
gecco <- gecco[,-1]

gecco <- gecco[colnames(adata),]
gecco[is.na(gecco)] <- 0

gecco <- gecco[,colSums(gecco)>1]
adata@meta.data <- cbind(adata@meta.data, gecco)

FeaturePlot(adata, features = colnames(gecco), ncol=4)
#ggsave(file.path(plotDirAll, "saliva_amrfinder.png"), width = 10, height = 12)



############ Nicer looking UMAPs

list_plots <- list()
for(cur_feature in colnames(gecco)){
  
  toplot <- data.frame(
    cnt=adata@meta.data[,cur_feature]>0,
    x=adata@reductions$kraken_umap@cell.embeddings[,1],
    y=adata@reductions$kraken_umap@cell.embeddings[,2]
  )
  toplot <- toplot[order(toplot$cnt, decreasing = FALSE),]
  
  list_plots[[cur_feature]] <- ggplot(toplot, aes(x,y, color=cnt)) +
    geom_point() + scale_color_manual(values=c("gray","red"), guide="none")+    
    theme(legend.position = "none")+
    theme_bw() + 
    xlab(cur_feature) + ylab("")
  
}
ptot <- egg::ggarrange(plots=list_plots, ncol=4)
ggsave(plot = ptot, file.path(plotDirAll, "saliva_gecco.png"), width = 10, height = 12)



############ Selected BGC UMAPs

toplot <- data.frame(
  x=adata@reductions$kraken_umap@cell.embeddings[,1],
  y=adata@reductions$kraken_umap@cell.embeddings[,2]
)
toplot$category <- "other"
toplot$category[adata@meta.data[,"GCF0000004"]>0] <- "GCF0000004"
toplot$category[adata@meta.data[,"GCF0000005"]>0] <- "GCF0000005"

toplot <- toplot[order(toplot$category, decreasing = TRUE),]

ggplot(toplot, aes(x,y, color=category)) +
  geom_point() + scale_color_manual(values=c("red","blue","gray"), guide="none")+    
  #  geom_point() + scale_color_manual(values=c("red","blue","gray"))+    
  theme(legend.position = "none")+
  theme_bw() + 
  xlab("KRAKEN1") + ylab("KRAKEN2")

ggsave(file.path(plotDirAll, "saliva_gecco_selected.png"), width = 5, height = 5)







################################################################################
################################## GTDB for saliva #############################
################################################################################

## Load GTDB metadata
gtdb_tab <- read.csv("/home/mahogny/github/scwgs/saliva_contigs_gtdbtk/gtdbtk.bac120.summary.tsv", sep="\t")
gtdb_classification <- as.data.frame(str_split_fixed(gtdb_tab$classification,";",7))
colnames(gtdb_classification) <- c("domain","phylum","class","order","family","genus","species")
rownames(gtdb_classification) <- gtdb_tab$user_genome

adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))

#Add GTDB data to kraken object
adata <- adata[,colnames(adata) %in% rownames(gtdb_classification)]
gtdb_classification <- gtdb_classification[colnames(adata),]
gtdb_classification[is.na(gtdb_classification)] <- "NA"
gtdb_classification[gtdb_classification==""] <- "NA"

adata$gtdb_domain <- gtdb_classification$domain
adata$gtdb_phylum <- gtdb_classification$phylum
adata$gtdb_class <- gtdb_classification$class
adata$gtdb_order <- gtdb_classification$order
adata$gtdb_family <- gtdb_classification$family
adata$gtdb_genus <- gtdb_classification$genus
adata$gtdb_species <- gtdb_classification$species

table(gtdb_classification$species)

DimPlot(adata, group.by = c("gtdb_domain")) + xlab("KRAKEN1") + ylab("KRAKEN2")

DimPlot(adata, group.by = c("gtdb_family"), label = TRUE) + 
  xlab("KRAKEN1") + ylab("KRAKEN2")

DimPlot(adata, group.by = c("gtdb_genus"), label = TRUE) + 
  xlab("KRAKEN1") + ylab("KRAKEN2") #+
#scale_color_manual(values = c(khroma::color("muted")(n = 11),"#DDDDDD"))

DimPlot(adata, group.by = c("gtdb_species"), label = TRUE, label.size = 3) + 
  xlab("KRAKEN1") + ylab("KRAKEN2")+theme(legend.position = "none")

DimPlot(adata, group.by = c("gtdb_species"), label = TRUE, label.size = 3) + 
  xlab("KRAKEN1") + ylab("KRAKEN2")#+theme(legend.position = "none")

DimPlot(adata, group.by = c("gtdb_species"), label = FALSE, label.size = 3) + 
  xlab("KRAKEN1") + ylab("KRAKEN2")#+theme(legend.position = "none")

DimPlot(adata, group.by = c("gtdb_species"), label = FALSE) + 
  xlab("KRAKEN1") + ylab("KRAKEN2")+theme(legend.position = "none")










################################################################################
######### Analysis of reads: Host% DNA ######################################### saliva only. overestimates!
################################################################################





cnt <- ReadBascetCountMatrix(bascetRoot,"chromcount", verbose=FALSE)

cnt_cs <- colSums(cnt$X)
cnt_xanthan <- sum(cnt_cs[str_starts(names(cnt_cs),"NZ_")])
cnt_human <- sum(cnt_cs[!str_starts(names(cnt_cs),"NZ_")])
cnt_other <- sum(cnt$obs$`_unmapped`)


df_src <- rbind(
  data.frame(src="Human", cnt=cnt_human),
  data.frame(src="Xanthan", cnt=cnt_xanthan),
  data.frame(src="Other", cnt=cnt_other)
)
df_src$frac <- df_src$cnt/sum(df_src$cnt)*100

df_src

ggplot(df_src, aes(frac, src)) + 
  geom_bar(stat = "identity", fill="#85C1E9", color="#85C1E9") + 
#  geom_bar(stat = "identity", fill="#8E1382", color="#8E1382") + 
  coord_flip() + 
  theme_bw() +
  xlab("% reads") + 
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path(plotDirAll,"saliva_host_fraction.svg"), width = 1.2, height = 3)







