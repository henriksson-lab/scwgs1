# bedtools bamtobed -i sorted.9.bam | cut -f1-3 > conv.bed

#/husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/conv.bed

#sort -R conv.bed -o random.bed

#sort -R

# /husky/henriksson/atrandi/bwa_ref/atcc/ref10_xp.fa

# bedtools bamtobed -i /husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/sorted.9.bam | cut -f1-3

# grep ATCC_10987_contig_1

############
# bedtools bamtobed -i /husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/unsorted.9.bam | cut -f1-3 | grep ATCC_10987_contig_1 > /husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/conv.bed
# bedtools bamtobed -i /husky/henriksson/atrandi/v6_260116_mock_pta/unsorted.9.bam | cut -f1-3 | grep ATCC_10987_contig_1 > /husky/henriksson/atrandi/v6_260116_mock_pta/conv.bed
# sort -R /husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/conv.bed > /husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/random.bed
# sort -R /husky/henriksson/atrandi/v6_260116_mock_pta/conv.bed > /husky/henriksson/atrandi/v6_260116_mock_pta/random.bed


library(GenomicRanges)
library(ggplot2)


subsamp_one <- function(fname) {
  all_len <- NULL
  for(numread in c(100,500,1000,5000,10000,20000,30000,40000,50000,100000,200000,300000,400000)) {
    #numread <- 100
    print(numread)
    subs <- readLines(fname, n=numread)
    subs <- stringr::str_split_fixed(subs,"\t",3)
    
    df <- data.frame(
      start=as.integer(subs[,2]),
      end=as.integer(subs[,3])
    )
    df$seqnames <- "b"
    
    read_bp <- sum(df$end- df$start)
    
    gr <- GRanges(seqnames = Rle(df$seqnames),
                  ranges = IRanges(start = df$start,
                                   end = df$end))
    #print(gr)
    gr <- reduce(gr)
    totlen <- sum(end(gr)-start(gr))
    
    all_len <- rbind(
      all_len,
      data.frame(numread=numread, read_bp=read_bp, coverage=totlen)
    )  
  }  
  all_len
}

subsamp_mda <- subsamp_one("/husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/random.bed")
#subsamp_mda <- subsamp_one("/husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/conv.bed")
subsamp_mda$chem <- "mda"

subsamp_pta <- subsamp_one("/husky/henriksson/atrandi/v6_260116_mock_pta/random.bed")
#subsamp_pta <- subsamp_one("/husky/henriksson/atrandi/v6_260116_mock_pta/conv.bed")
subsamp_pta$chem <- "pta"

#fname <- "/husky/henriksson/atrandi/v6_251128_jyoti_mock_bulk/conv.bed"


#ggplot(subsamp_mda, aes(read_bp, coverage)) + geom_line()
#ggplot(all_len, aes(numread, coverage)) + geom_line()

ggplot(rbind(subsamp_mda, subsamp_pta), aes(read_bp, coverage, color=chem)) + geom_line()


#how big is this thing??
