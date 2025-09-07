

library(ggplot2)

################################################################################
######################## Fig 1 xxx #############################################
################################################################################

df <- merge(
  data.frame(genome_size=seq(from=1000000, to=5000000, by=10000)),
  data.frame(num_cell=seq(from=1000, to=1000000, by=1000)),
)

totreads <- 25e9
readlen <- 150*2

df$depth <- totreads*readlen/(df$genome_size*df$num_cell)

#ggplot(df, aes(genome_size, log10(numcell), z=depth)) +  geom_contour_filled()
#ggplot(df, aes(genome_size, log10(numcell), z=log10(depth))) +  geom_contour_filled() + theme_bw()
#ggplot(df, aes(genome_size, num_cell, z=log10(depth))) +  geom_contour_filled() + theme_bw() + scale_y_log10()
ggplot(df, aes(genome_size/1e6, num_cell/1000, z=log10(depth))) +  geom_contour_filled() + theme_bw() + scale_y_log10()
#ggplot(df, aes(genome_size, num_cell, z=(depth))) +  geom_contour_filled() + theme_bw() + scale_y_log10()
ggsave(file.path(plotDirAll, "fig1_coverage.svg"), width = 4, height = 2)


#sort(df$depth)
#tail(sort(df$depth))

ggplot(df, aes(genome_size, num_cell, z=(depth))) +  geom_contour_filled(breaks=c(1,2, 5,10,30,100,300,1000,3000,5000,10000)) + theme_bw() + scale_y_log10()


df[df$num_cell==1e6 & df$genome_size==3e6,]
df[df$num_cell==1e6 & df$genome_size==5e6,]

################################################################################
######################## multiplet rate 10k cells ##############################
################################################################################


bcspace <- 24**4

table(table(sample(1:bcspace, 10000, replace = TRUE)))
#9674  163 
163/(9674 + 163 )
