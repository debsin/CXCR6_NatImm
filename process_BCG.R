library(ggplot2)
library(stringr)
library(ggpubr)
library("GEOquery")

gse=getGEO(filename="BCG/GSE124218_series_matrix.txt.gz")

colnames(pData(phenoData(gse)))
# title, source_name_ch1  characteristics_ch1 characteristics_ch1.1 characteristics_ch1.2 characteristics_ch1.3
count = read.csv("BCG/GSE124218_BCG_BM_cirovic_counts.txt", row.names = 1)

meta = data.frame(full = colnames(count))

a = meta$full[1]
a
a = str_split(meta$full, "_")

d = do.call(rbind, a)
head(d)
colnames(d) = c("condition", "donor", "cell", "tissue", "timepoint")
meta = cbind(meta, d)
head(meta)

plot_df = meta
plot_df$condition = factor(plot_df$condition, levels = c("Ctrl", "BCG"))
plot_df$CXCR6 =  t(count["CXCR6", ])
head(plot_df)
g1 = ggplot(plot_df, aes(x = condition, y = CXCR6, fill = condition))+
  facet_wrap(~timepoint)+
  scale_fill_hue(direction = -1)+
  geom_boxplot(outlier.shape = NA, show.legend = F)+
  geom_jitter(width = 0.1, show.legend = F)+
  stat_compare_means( method = "anova", label = "p.format")+
  theme_classic()
  
g1

am = list()
am$mice = readRDS("TRM/Mice Lung/top_mice_am.rds")
am$human_homolog = mus2humanGENE(sort(am$mice))

am_module = scale(apply(count[intersect(am$human_homolog,rownames(count)), ], 2, function(x) median(x)))
plot_df$am = am_module[,1]

g2 = ggplot(plot_df, aes(x = condition, y = am, fill = condition))+
  facet_wrap(~timepoint)+
  geom_boxplot(outlier.shape = NA, show.legend = F)+
  scale_fill_hue(direction = -1)+
  geom_jitter(width = 0.1, show.legend = F)+
  stat_compare_means( method = "anova", label = "p.format")+
  ylab("AM signature from mice")+
  theme_classic()

gridExtra::grid.arrange(g2, ncol = 1, top=textGrob("HSPC BM Dataset, Cirovic B et al."))

library(gridExtra)
library(grid)
pdf("BM_data_am_signature.pdf", width = 5, height = 5)
gridExtra::grid.arrange(g2, ncol = 1, top=textGrob("HSPC BM Dataset, Cirovic B et al."))
dev.off()


# am = readRDS("mice_am_sig.rds")
# am_module = scale(apply(count[am$human_homolog, ], 2, function(x) mean(x)))
# plot_df$am = am_module
# 
# g2 = ggplot(plot_df, aes(x = condition, y = am, fill = condition))+
#   facet_wrap(~timepoint)+
#   geom_boxplot(outlier.shape = NA, show.legend = F)+
#   scale_fill_hue(direction = -1)+
#   geom_jitter(width = 0.1, show.legend = F)+
#   stat_compare_means( method = "t.test", label = "p.format")+
#   ylab("am signature from mice")+
#   theme_classic()
# 
# library(gridExtra)
# library(grid)
# pdf("BM_data.pdf", width = 5, height = 5)
# gridExtra::grid.arrange(g2, ncol = 1, top=textGrob("HSPC BM Dataset, Cirovic B et al."))
# dev.off()
