dat <- read.table("/home/mahogny/readdis.tsv")




setwd("/husky/julian")


alldat <- list()
for(i in 1:20) {
  dat <- read.table("/husky/julian/debarcoded.1.tirp.gz.hist")
  alldat[[i]] <- dat
}
alldat <- do.call(rbind, alldat)
colnames(alldat) <- c("cell","cnt")

sumdat <- sqldf::sqldf("select cell, sum(cnt) as cnt from alldat group by cell")

sum(sumdat$cnt)/1e9  #1B reads

sumdat <- sumdat[order(sumdat$cnt, decreasing = TRUE), ]
head(sumdat)
sumdat$index <- 1:nrow(sumdat)

library(ggplot2)
ggplot(sumdat, aes(log10(index), log10(cnt))) + geom_line()

10**5.5 * (150+100)/1e6
