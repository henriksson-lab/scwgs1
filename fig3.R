


################################################################################
########################## Fig 3 umap matrix ###################################
################################################################################


#adata1 <- readRDS(file.path("/husky/henriksson/atrandi//v4_wgs_novaseq1","cache_adata_align.RDS"))

######### 
######### simulated mock data
######### 
bdata1 <- readRDS("/husky/henriksson/atrandi/simulated4/cache_adata_align.RDS")
bdata2 <- readRDS("/husky/henriksson/atrandi/simulated4/cache_adata_kraken.RDS")
bdata3 <- readRDS("/husky/henriksson/atrandi/simulated4/cache_adata_infokmer.RDS")
bdata4 <- readRDS("/husky/henriksson/atrandi/simulated4/cache_adata_cs.RDS")

colnames(bdata1) <- str_remove(colnames(bdata1),"BASCET_")

#take cell names from bdata1
bdata2$species_aln_short <- bdata1@meta.data[colnames(bdata2),]$species_aln_short
bdata3$species_aln_short <- bdata1@meta.data[colnames(bdata3),]$species_aln_short
bdata4$species_aln_short <- bdata1@meta.data[colnames(bdata4),]$species_aln_short

######### 
######### real mock data
######### 

adata1 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/cache_adata_align.RDS")
adata2 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/cache_adata_kraken.RDS")
adata3 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/cache_adata_infokmer.RDS")
adata4 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/cache_adata_cs.RDS")

adata1 <- adata1[,adata1$nCount_chrom_cnt>5000] #removes a lot of crap

colnames(adata1) <- str_remove(colnames(adata1),"BASCET_")


#take cell names from adata1
adata2$species_aln_short <- adata1@meta.data[colnames(adata2),]$species_aln_short
adata3$species_aln_short <- adata1@meta.data[colnames(adata3),]$species_aln_short
adata4$species_aln_short <- adata1@meta.data[colnames(adata4),]$species_aln_short



do_label <- FALSE
a_p1 <- DimPlot(object = adata1, label = do_label, group.by = "species_aln_short") + 
  xlab("BWA1") + ylab("BWA2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)
a_p1
#ggsave(plot = a_p1, file.path(plotDirAll, "fig3_one_umap.svg"), width = 4, height = 4)

a_p2 <- DimPlot(object = adata2[,!is.na(adata2$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("KRAKEN1") + ylab("KRAKEN2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

a_p3 <- DimPlot(object = adata3[,!is.na(adata3$species_aln_short)], label = do_label, group.by = "species_aln_short") + #one cluster is a bit meah
  xlab("INFOKMER1") + ylab("INFOKMER2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

a_p4 <- DimPlot(object = adata4[,!is.na(adata4$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("CS1") + ylab("CS2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)





b_p1 <- DimPlot(object = bdata1, label = do_label, group.by = "species_aln_short") + 
  xlab("BWA1") + ylab("BWA2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

b_p2 <- DimPlot(object = bdata2[,!is.na(bdata2$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("KRAKEN1") + ylab("KRAKEN2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

b_p3 <- DimPlot(object = bdata3[,!is.na(bdata3$species_aln_short)], label = do_label, group.by = "species_aln_short") + #one cluster is a bit meah
  xlab("INFOKMER1") + ylab("INFOKMER2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

b_p4 <- DimPlot(object = bdata4[,!is.na(bdata4$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("CS1") + ylab("CS2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)


#ptot <- egg::ggarrange(p1,p2,p3,p4, nrow=1)

ptot <- egg::ggarrange(
  b_p1+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  b_p2+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  b_p3+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  b_p4+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  
  a_p1+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  a_p2+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  a_p3+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  a_p4+theme(legend.position = "none", axis.title=element_blank(),axis.text=element_blank(), axis.ticks=element_blank()),
  
  nrow=2)
ptot
#ggsave(plot = ptot, file.path(plotDirAll, "fig3_umaps.svg"), width = 16, height = 8)
ggsave(plot = ptot, file.path(plotDirAll, "fig3_umaps.png"), width = 16*0.5, height = 8*0.5)






ptot <- egg::ggarrange(
  b_p1+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  b_p2+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  b_p3+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  b_p4+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  
  a_p1+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  a_p2+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  a_p3+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  a_p4+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+xlab("")+ylab(""),
  
  nrow=2)
ptot
ggsave(plot = ptot, file.path(plotDirAll, "fig3_umaps.png"), width = 16, height = 8)
ggsave(plot = ptot, file.path(plotDirAll, "fig3_umaps.png"), width = 16*0.5, height = 8*0.5)


FeaturePlot(object = adata1, features = "nCount_species_cnt") + 
  xlab("BWA1") + ylab("BWA2") + 
  ggtitle("") 






################################################################################
########################## Fig 3 JL graph ######################################
################################################################################

PlotJohnsonLindenstraussMinDim(list_eps = c(0.1,0.2,0.3))
ggsave(file.path(plotDirAll, "fig3_JL.svg"), width = 4, height = 2)





################################################################################
########################## Fig 3 kmer hisogram #################################
################################################################################

dataset_name <- "v4_wgs_novaseq3"
bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)


p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")

kmerHist <- BascetReadMinhashHistogram(bascetRoot)

### KMER count histogram
kmerHist$rank <- 1:nrow(kmerHist)
ggplot(kmerHist[sample(1:nrow(kmerHist),min(30000, nrow(kmerHist))),], aes(rank, cnt)) +   
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("Count")+
  xlab("Rank")
ggsave(file.path(plotDirAll,"fig3_info_kmerHist.svg"), width = 3, height = 3)
#ggsave(file.path(plotDir,"info_kmerHist.png"), width = 5, height = 5)




################################################################################
########################## Fig 3 kmer on umap ##################################
################################################################################

dataset_name <- "v4_wgs_novaseq3"
bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)

adata <- readRDS(file.path(bascetRoot,"cache_adata_infokmer.RDS"))

FeaturePlot(adata, features = "AGATCGGAAGAGCACACGTCTGAACTCCAGT") + 
  xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDirAll,"fig3_info_kmer_ex.pdf"), width = 4, height = 4)


