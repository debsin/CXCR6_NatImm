

df = as.data.frame(t(assay(vsd)))
df$sample = ntd$population
df$norm = "vst"

df  = melt(df)
head(df)

df2= as.data.frame(t(assay(ntd)))
df2$sample = ntd$population
df2$norm = "ntd"

df2  = melt(df2)

df1 = rbind(df, df2)
head(df1)

ggplot(df2, aes(x=value,  color = sample)) + 
  geom_density()


library(vsn)
meanSdPlot(assay(ntd), ranks=FALSE)
x <- assay(ntd)[,1]
y <- assay(ntd)[,2]
plot(.5*(x + y), y - x,
     cex=.5, col=rgb(0,0,0,.1), pch=20)
abline(h=0, col="red", lwd=3)
meanSdPlot(assay(vsd), ranks=FALSE)
x <- assay(vsd)[,1]
y <- assay(vsd)[,2]
plot(.5*(x + y), y - x,
     cex=.5, col=rgb(0,0,0,.1), pch=20)
abline(h=0, col="red", lwd=3)

sq.diffs <- (assay(ntd)[,2] - assay(ntd)[,1])^2
ave <- .5 * (assay(ntd)[,2] + assay(ntd)[,1])
plot(ave, sq.diffs)

sum(sq.diffs[ave < 5])/sum(sq.diffs[ave > 5])
sum(sq.diffs[ave < 3])/sum(sq.diffs[ave > 3])

library(matrixStats)
rv <- rowVars(assay(vsd))
o <- order(rv,decreasing=TRUE)
dists <- dist(t(assay(vsd)[head(o,500),]))


hc <- hclust(dists)
plot(hc, labels=vsd$population)


library(magrittr)
table(vsd$population)

idx <- vsd$population == "CD8posCD103pos"
vsd2 <- vsd
vsd2$patient %<>% factor


dists <- dist(t(assay(vsd)[head(o,100),]))
hc <- hclust(dists)
dend <- as.dendrogram(hc)


suppressPackageStartupMessages(library(dendextend))
library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
o.dend <- order.dendrogram(dend)
labels(dend) <- vsd2$population[o.dend]
labels_colors(dend) <- as.integer(vsd2$population[o.dend])
plot(dend)

plotPCA(vsd2, intgroup=c("population"))


vsd2$cluster <- cutree(hc, k=5)
head(vsd2$cluster)

table(vsd2$cluster)
labels_colors(dend) <- as.integer(vsd2$cluster[o.dend])
plot(dend)
