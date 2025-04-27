################################################################################ 
################## Prepare info about mock community chroms ####################
################################################################################ 

map_strain2gram <- read.csv("/husky/fromsequencer/240701_wgs_atcc1/all/straintype.csv")

get_map_seq2strain <- function(){
  map_seq2strain <- NULL
  refdir <- "/husky/fromsequencer/240809_novaseq_wgs1/trimmed/ref10/separate"
  for(f in list.files(refdir, pattern = "*.fasta")){
    print(f)
    thel <- readLines(file.path(refdir,f))
    thel <- thel[str_starts(thel,">")]
    onedf <- data.frame(line=thel)
    onedf$strain <- f
    map_seq2strain <- rbind(map_seq2strain, onedf)
  }
  map_seq2strain$id <- str_split_fixed(str_sub(map_seq2strain$line,2)," ",2)[,1]
  map_seq2strain$name <- str_split_fixed(str_sub(map_seq2strain$line,2)," ",2)[,2]
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain,".fasta")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," chr1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," chr2")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," pRSPH01")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," pMP1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," pCP1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," NC_003909.8")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," NC_005707.1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," NC_005707.1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," plasmid")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," CP086328.1")
  map_seq2strain$strain <- str_remove_all(map_seq2strain$strain," CP086329.1")
  
  map_seq2strain <- merge(map_seq2strain,map_strain2gram)
  map_seq2strain <- map_seq2strain[,colnames(map_seq2strain)!="line"]
  
  idx <- read.table(pipe("samtools idxstats /husky/fromsequencer/241206_novaseq_wgs3/trimmed/sorted.bam"),sep="\t")
  colnames(idx) <- c("id","len","mapped","unmapped")
  map_seq2strain <- merge(idx[,c("id","len")], map_seq2strain)  #note, only using id and len
  
  map_seq2strain$id <- stringr::str_replace_all(map_seq2strain$id,stringr::fixed("_"), "-") #### for compatibility with signac
  
  map_seq2strain
}

map_seq2strain <- get_map_seq2strain()[,c("id","len","strain")]
write.csv(map_seq2strain, "~/github/scwgs/map_seq2strain.csv", row.names = FALSE)

