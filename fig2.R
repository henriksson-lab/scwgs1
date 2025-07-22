

################################################################################
################# Fig 2 kneeplots for different lysis ########################## per species
################################################################################


adata1 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq1/cache_adata_align.RDS")
DefaultAssay(adata1) <- "species_cnt"

adata2 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq2/cache_adata_align.RDS")
DefaultAssay(adata2) <- "species_cnt"

adata3 <- readRDS("/husky/henriksson/atrandi/v4_wgs_novaseq3/cache_adata_align.RDS")
DefaultAssay(adata3) <- "species_cnt"

p1 <- KneeplotPerSpecies(adata1) + scale_color_manual(values = list_color_mock)
p2 <- KneeplotPerSpecies(adata2) + scale_color_manual(values = list_color_mock)
p3 <- KneeplotPerSpecies(adata3) + scale_color_manual(values = list_color_mock)

ptot <- egg::ggarrange(
  p1,p2,p3,
  nrow=1)
ptot

ggsave(plot = ptot, file.path(plotDirAll, "fig2_strain_kneeplots.svg"), width = 20, height = 4)



################################################################################
################# Fig 2 kneeplots for different lysis ########################## all species as a whole
################################################################################

p1 <- plot_one_kneeplot_all(adata1)
p2 <- plot_one_kneeplot_all(adata2)
p3 <- plot_one_kneeplot_all(adata3)

ptot <- egg::ggarrange(
  p1,p2,p3,
  nrow=1)
ptot

ggsave(plot = ptot, file.path(plotDirAll, "fig2_allread_kneeplots.svg"), width = 8, height = 3)



################################################################################
################# Fig 2 count per species ########################## 
################################################################################

#10k reads make sense?

df <- rbind(
  data.frame(
    species=adata1@meta.data$species_aln[adata1$nCount_species_cnt>10000],
    rnd=1
  ),
  data.frame(
    species=adata2@meta.data$species_aln[adata2$nCount_species_cnt>10000],
    rnd=2
  ),
  data.frame(
    species=adata3@meta.data$species_aln[adata3$nCount_species_cnt>10000],
    rnd=3
  )
)

#df <- data.frame(species=adata3@meta.data$species_aln[adata3$nCount_species_cnt>10000])

lysiseff <- merge(
  sqldf::sqldf("select species, rnd, count(*) as cnt from df group by species, rnd"),
  sqldf::sqldf("select rnd, count(*) as cnt_tot from df group by rnd")
)
lysiseff$frac <- lysiseff$cnt / lysiseff$cnt_tot * 100
ggplot(lysiseff, aes(frac, species, fill=paste("",rnd))) + 
  geom_bar(stat="identity", position="dodge", position = position_fill(reverse = TRUE))


plot_one_speciesfrac <- function(adata1){
  df <- data.frame(species=adata1@meta.data$species_aln[adata1$nCount_species_cnt>10000])
  df <- sqldf::sqldf("select species, count(*) as cnt from df group by species")
  df$frac <- df$cnt / sum(df$cnt) * 100
  df$species <- stringr::str_split_i(df$species," \\(",1)
  #df
  
  df$gram <- "pos"
  df$gram[stringr::str_detect(df$species,"coli")] <- "neg"
  df$gram[stringr::str_detect(df$species,"sphaer")] <- "neg"
  
  df$species <- fct_rev(df$species)
  #levels(df$species) <- stringr::str_split_i(levels(df$species)," \\(",1)
  
  ggplot(df, aes(species, frac, fill=gram)) + 
    geom_bar(stat="identity") +   #, position = position_stack(reverse = TRUE)
    coord_flip() + 
    ylab("Fraction %") +
    xlab("Species") +
    theme_bw() #position = position_fill(reverse = TRUE)
}


ptot <- egg::ggarrange(
  plot_one_speciesfrac(adata1),
  plot_one_speciesfrac(adata2),
  plot_one_speciesfrac(adata3),
  nrow=1)
ptot
ggsave(plot = ptot, file.path(plotDirAll, "fig2_allread_kneeplots.svg"), width = 16, height = 3)

#sqldf::sqldf("select rnd, count(*) as cnt_tot from df group by rnd")




##########
########## Barnyard plot for mock R3
##########

bascetRoot <- "/husky/henriksson/atrandi/v4_wgs_novaseq3/"
adata <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))


bp <- data.frame(
  maxc = MatrixGenerics::colMaxs(adata@assays$species_cnt$counts),
  totc = colSums(adata@assays$species_cnt$counts)
)
bp$restc <- bp$totc - bp$maxc

ggplot(bp, aes(maxc/1e6,restc/1e6)) + 
  geom_point() +
  theme_bw() +
  xlab("Dominant species count (M reads)") +
  ylab("Other species count (M reads)") +
  xlim(0,6) + 
  ylim(0,6)
ggsave(file.path(plotDirAll,"fig2_alignment_barnyard.pdf"), width = 4, height = 4)


####### Statistics on mixing

sum(bp$maxc)/sum(bp$totc)
sum(bp$restc)/sum(bp$totc)
nrow(bp)









################################################################################
######### Coverage of cell #####################################################
################################################################################

# Take top cell, show typical coverage
# Do subsampling analysis, coverage vs depth

#adata from alignment analysis must be loaded. novaseq1



plot(sort(as.double(adata$nCount_chrom_cnt[adata$species_aln_short=="Escherichia coli"])))

onemat <- ReadBascetCountMatrix_one(file.path(bascetRoot,"chromcount.1.h5"), verbose=FALSE)
rownames(onemat)

### Pick cells to analyze
submeta <- adata@meta.data[colnames(adata) %in% rownames(onemat),]
set.seed(666)
ind_samples <- sample(1:nrow(submeta), prob = submeta$nCount_chrom_cnt, size = 20)
submeta <- submeta[ind_samples,]
submeta <- submeta[order(submeta$nCount_chrom_cnt),]
submeta
table(submeta$species_aln_short)


### For each cell
list_alldf <- list()
for(cur_cell_num in 1:nrow(submeta)){
  
  ## Get the aligned reads for this cell
  cellid <- rownames(submeta)[cur_cell_num]
  lines_for_cell <- system(paste("grep",cellid,"/husky/henriksson/atrandi/v4_wgs_novaseq1/unsorted_aligned.1.bed"), intern = TRUE)
  zz <- textConnection(lines_for_cell)
  cellbedfile <- read.table(zz)[,1:3]
  close(zz)
  colnames(cellbedfile) <- c("seqname","start","end")
  
  ## For simplicity, we only consider the major chromosome
  cellbedfile <- cellbedfile[cellbedfile$seqname == names(sort(table(cellbedfile$seqname), decreasing = TRUE))[1],]
  
  ## Figure out length of chromosome
  len_of_chrom <- mapSeq2strain$len[stringr::str_replace_all(mapSeq2strain$id,"-","_")==cellbedfile$seqname[1]]
  cellid_strain <- mapSeq2strain$strain[stringr::str_replace_all(mapSeq2strain$id,"-","_")==cellbedfile$seqname[1]]
  
  ## Subsample reads
  for(frac in pracma::logseq(0.001, 1, n = 30)) {
    print(paste(cur_cell_num,cellid,frac))  #future: can set max 20x coverage
    
    ###### Actual coverage
    
    ## Compute overlaps
    sampled_cellbedfile <- cellbedfile[sample(1:nrow(cellbedfile),size=round(nrow(cellbedfile)*frac)),]
    gr <- makeGRangesFromDataFrame(sampled_cellbedfile) #cellbedfile[sample(1:nrow(cellbedfile),size=round(nrow(cellbedfile)*frac)),])
    grcov <- coverage(gr)[[1]]
    unique_read_cov <- sum(grcov@lengths[grcov@values>0])
    tot_read_cov <- sum(gr@ranges@width)
    
    ###### Theoretically best coverage
    
    #Bootstrap coordinates
    #sampled_cellbedfile <- cellbedfile[sample(1:nrow(cellbedfile),size=round(nrow(cellbedfile)*frac)),]
    sampled_cellbedfile$length <- sampled_cellbedfile$end - sampled_cellbedfile$start
    sampled_cellbedfile$start <- round(runif(nrow(sampled_cellbedfile),min = 1, max = len_of_chrom - sampled_cellbedfile$length))
    sampled_cellbedfile$end <- sampled_cellbedfile$start + sampled_cellbedfile$length
    
    ## Compute overlaps
    gr <- makeGRangesFromDataFrame(sampled_cellbedfile)
    grcov <- coverage(gr)[[1]]
    sim_unique_read_cov <- sum(grcov@lengths[grcov@values>0])
    sim_tot_read_cov <- sum(gr@ranges@width)
    
    onedf <- data.frame(
      cellid=cellid,
      strain=cellid_strain,
      num_read=length(gr),
      len_of_chrom=len_of_chrom,
      
      unique_read_cov=unique_read_cov/len_of_chrom,
      tot_read_cov=tot_read_cov/len_of_chrom,
      
      sim_unique_read_cov=sim_unique_read_cov/len_of_chrom,
      sim_tot_read_cov=sim_tot_read_cov/len_of_chrom  #should be the same? almost.
    )
    
    list_alldf[[paste(cellid,frac)]] <- onedf
  }  
}
alldf <- do.call(rbind, list_alldf)

if(FALSE){
  ggplot(alldf, aes(tot_read_cov, unique_read_cov,group=cellid, color=strain)) + 
    geom_line() + 
    theme_bw() +
    scale_x_log10() +
    xlab("x Coverage") +
    ylab("Fraction of overlapped genome")
  
  
  ggplot(alldf, aes(num_read, unique_read_cov,group=cellid, color=strain)) + 
    geom_line() + 
    theme_bw() +
    xlab("Number of reads") +
    ylab("Fraction of overlapped genome")
  #strain
  
  
  ###
  ggplot(alldf, aes(sim_tot_read_cov, sim_unique_read_cov,group=cellid, color=strain)) + 
    geom_line() + 
    theme_bw() +
    scale_x_log10() +
    xlab("Bootstrapped, x Coverage") +
    ylab("Bootstrapped, fraction of overlapped genome")
  
}

### Make a plot with combined efficiencies, and optimal efficiency
#alldf[which.max(alldf$tot_read_cov),]
cellid_max_cov <- alldf[which.max(alldf$tot_read_cov),]$cellid
df_with_sim <- rbind(
  alldf[,c("tot_read_cov","unique_read_cov","cellid","strain")],
  data.frame(
    tot_read_cov=alldf$sim_tot_read_cov[alldf$cellid==cellid_max_cov],
    unique_read_cov=alldf$sim_unique_read_cov[alldf$cellid==cellid_max_cov],
    cellid="Optimal sample",
    strain="Optimal sample"
  )
)
ggplot(df_with_sim, aes(tot_read_cov, unique_read_cov,group=cellid, color=strain)) + 
  geom_line() + 
  theme_bw() +
  #  scale_x_log10() +
  xlim(0,50) +
  xlab("x Coverage") +
  ylab("Fraction of overlapped genome")
ggsave(file.path(plotDir,"coverage_saturation.pdf"), width = 15, height = 3)





