################################################################################ 
################## Prepare info about mock community chroms ####################
################################################################################ 

map_strain2gram <- read.csv("/husky/fromsequencer/240701_wgs_atcc1/all/straintype.csv")

get_mapSeq2strain <- function(){
  mapSeq2strain <- NULL
  refdir <- "/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/separate"
  for(f in list.files(refdir, pattern = "*.fasta")){
    print(f)
    thel <- readLines(file.path(refdir,f))
    thel <- thel[str_starts(thel,">")]
    onedf <- data.frame(line=thel)
    onedf$strain <- f
    mapSeq2strain <- rbind(mapSeq2strain, onedf)
  }
  mapSeq2strain$id <- str_split_fixed(str_sub(mapSeq2strain$line,2)," ",2)[,1]
  mapSeq2strain$name <- str_split_fixed(str_sub(mapSeq2strain$line,2)," ",2)[,2]
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain,".fasta")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," chr1")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," chr2")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," pRSPH01")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," pMP1")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," pCP1")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," NC_003909.8")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," NC_005707.1")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," NC_005707.1")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," plasmid")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," CP086328.1")
  mapSeq2strain$strain <- str_remove_all(mapSeq2strain$strain," CP086329.1")
  
  mapSeq2strain <- merge(mapSeq2strain,map_strain2gram)
  mapSeq2strain <- mapSeq2strain[,colnames(mapSeq2strain)!="line"]
  
  idx <- read.table(pipe("samtools idxstats /husky/fromsequencer/241206_novaseq_wgs3/trimmed/sorted.bam"),sep="\t")
  colnames(idx) <- c("id","len","mapped","unmapped")
  mapSeq2strain <- merge(idx[,c("id","len")], mapSeq2strain)  #note, only using id and len
  
  mapSeq2strain$id <- stringr::str_replace_all(mapSeq2strain$id,stringr::fixed("_"), "-") #### for compatibility with signac
  
  mapSeq2strain
}

mapSeq2strain <- get_mapSeq2strain()[,c("id","len","strain")]
write.csv(mapSeq2strain, "~/github/scwgs/mapSeq2strain.csv", row.names = FALSE)

