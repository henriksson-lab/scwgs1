library(seqinr)


get_gc_frac <- function(fname) {
  onefa <- seqinr::read.fasta(fname)
  tot_cg <- 0
  tot_len <- 0
  for(n in names(onefa)) {
    onecontig <- onefa[[n]]
    this_cg <- sum(onecontig %in% c("g","c"))
    this_len <- length(onecontig)
    tot_cg <- tot_cg + this_cg
    tot_len <- tot_len + this_len
  }
  tot_cg/tot_len
}


get_gc_frac("/husky/henriksson/atrandi/bwa_ref/atcc/ATCC-AGP-24087822e79a496f/genome_ATCC_10987_687931d9b06b4cb4/assembly/Bacillus_pacificus_ATCC_10987.fasta")

rootdir <- "/husky/henriksson/atrandi/bwa_ref/atcc/copied_out_assembly"
all_fname <- list.files(rootdir)
all_df <- NULL
for(fname in all_fname) {
  print(fname)

  df <- data.frame(
    fname=fname,
    gc=get_gc_frac(file.path(rootdir, fname))
  )  
  all_df <- rbind(all_df, df)
}
all_df[order(all_df$gc),]

# The proportion of mapped reads decreased for DNA genomes with higher GC content. 
# https://academic.oup.com/ismecommun/article/4/1/ycae024/7606642


# contig 6 AATTATGGTGCGACAAGAAGTCGCATTCCCACA..... GCTTACCCGACTAAATAAGCTG
# contig 7 TGAACTGGACCCGA .... TAACAAATTGGGTGCGGTTCA


