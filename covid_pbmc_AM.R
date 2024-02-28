library(Seurat)
library(org.Mm.eg.db)
library(dittoSeq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(reshape2)


cov_pbmc = readRDS("C:/Projects/PBMC_covid_sc_/covid_meta_AM_sig.rds")

colnames(cov_pbmc)

table(cov_pbmc$Status_on_day_collection)
table(cov_pbmc$Status_on_day_collection_summary)

cov_pbmc$phenotype = cov_pbmc$Status_on_day_collection_summary

plot_df = cov_pbmc[intersect(grep("CD14|CD16", cov_pbmc$initial_clustering),
                                     grep("Healthy|^Moder|Severe", cov_pbmc$phenotype)), ]

plot_df$condition  = factor(plot_df$phenotype , levels = c("Healthy","Moderate","Severe"))

gcovid = ggplot(plot_df, aes(x = condition, y = AM, fill = condition))+
  # facet_wrap(~predicted.celltype.l1, ncol = 6)+
  geom_violin(show.legend = F)+
  geom_jitter(width = 0.2, size = 0.1, show.legend = F,)+
  geom_boxplot(outlier.shape = NA, show.legend = F, width =0.5, fill = "#FFFFFFAA")+
  scale_fill_hue(direction = -1)+
  # geom_jitter(width = 0.1, show.legend = F)+
  stat_compare_means(ref.group = "Healthy",   method = "wilcox.test", label = "p.format")+
  ylab("AM signature from mice")+
  xlab("Monocytes")+
  theme_classic()

pdf("cov_pbmc_AM.pdf", width = 9, height = 7)
gridExtra::grid.arrange(gcovid, ncol = 1, top=textGrob("Covid PBMC Dataset,  Emily, et al. "))
dev.off()
