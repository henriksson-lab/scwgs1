qpcr_capsule <- read.csv("/home/mahogny/github/scwgs/qpcr/capsule.csv",sep="\t", row.names = "X")
rownames(qpcr_capsule)[3] <- "control"
qpcr_capsule <- melt(as.matrix(qpcr_capsule))
colnames(qpcr_capsule) <- c("condition","bacteria","ct")

qpcr_capsule$bacteria <- factor(qpcr_capsule$bacteria, levels=sort(unique(qpcr_capsule$bacteria), decreasing = TRUE))

p1 <- ggplot(qpcr_capsule, aes(bacteria, ct, fill=condition)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() + 
  theme_bw() + xlab("") +
  scale_fill_manual(values=khroma::color("muted")(n = 3))


########

qpcr_outside <- read.csv("/home/mahogny/github/scwgs/qpcr/outside.csv",sep="\t", row.names = "X")
rownames(qpcr_outside)[3] <- "control"
qpcr_outside <- melt(as.matrix(qpcr_outside))
colnames(qpcr_outside) <- c("condition","bacteria","ct")

qpcr_outside$bacteria <- factor(qpcr_outside$bacteria, levels=sort(unique(qpcr_outside$bacteria), decreasing = TRUE))

p2 <- ggplot(qpcr_outside, aes(bacteria, ct, fill=condition)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() + 
  theme_bw() + xlab("") +
  scale_fill_manual(values=khroma::color("muted")(n = 3))

ptot <- egg::ggarrange(p1,p2, nrow = 1)
ptot

ggsave(plot = ptot, file.path(plotDirAll, "qpcr.svg"), width = 7, height = 3)

