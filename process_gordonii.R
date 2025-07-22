
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
KrakenKneePlot(adata, groupby = "species", show_num_spec=10)


#Most common species
#Streptococcus gordonii

adata$log_cnt
adata$species=="Streptococcus gordonii"

met <- adata@meta.data[adata@meta.data$species=="Streptococcus gordonii",]
tail(met[order(met$log_cnt),])

#G1_E5_H7_E11  has high count, 7.540792e-05 human

met <- adata@meta.data[order(adata@meta.data$log_cnt, decreasing = TRUE),]
list_cells <- rownames(met)[1:10000]


################################################################################
##################### Output assembly for best gordonii cell ###################
################################################################################


#bascet_instance <- getBascetDockerImage("/Users/mahogny/Desktop/rust/bascet/docker_image")  #return ok on linux. wtf?
bascet_instance <- getBascetSingularityImage("/home/mahogny/github/bascet/singularity")
TestBascetInstance(bascet_instance)

skesa_file <- OpenBascet(bascetRoot, "skesa", bascet_instance)

BascetListFilesForCell(skesa_file, "G1_E5_H7_E11", bascet_instance = bascet_instance)

contig_file <- BascetReadFile(
  skesa_file,cellID = "G1_E5_H7_E11", filename = "contigs.fa",
  as="text"#, verbose = TRUE
)
str_length(contig_file)

CloseBascet(skesa_file)


writeLines(contig_file, file.path(plotDirAll,"gordonii.fa"))


gordonii_fa <- seqinr::read.fasta(file.path(plotDirAll,"gordonii.fa"))
#length(gordonii_fa)
#length(gordonii_fa$Contig_1_81.1742)

gordonii_fa_red <- gordonii_fa[sapply(gordonii_fa,length)>100]

sum(sapply(gordonii_fa_red,length))
sort(sapply(gordonii_fa_red,length),decreasing = TRUE)








################################################################################
##################### Output assembly for all cells ############################
################################################################################

#only 10k cells with most reads, list of cells from above

#bascet_instance <- getBascetDockerImage("/Users/mahogny/Desktop/rust/bascet/docker_image")  #return ok on linux. wtf?
bascet_instance <- getBascetSingularityImage("/home/mahogny/github/bascet/singularity")
TestBascetInstance(bascet_instance)


BascetDumpContigs(
  bascetRoot,
  inputName="skesa",
  listCells=list_cells,
  bascet_instance=bascet_instance,
  outputDir="/home/mahogny/github/scwgs/saliva_contigs"
)

