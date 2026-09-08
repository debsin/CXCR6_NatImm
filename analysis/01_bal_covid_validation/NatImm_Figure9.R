
library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
# map_signif_level <- c(`****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = Inf)

gene.path = "C:/Projects/external/jem_20190249_tables1/"

lungTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T)
lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[lungTRM_noTRM$`log2 fold change`>2 & lungTRM_noTRM$`P adj`<0.01]

immune.combined.sct =  readRDS("covid_nat_med/T_all.rds")
DefaultAssay(immune.combined.sct) = "integrated"

immune.combined.sct$condition = immune.combined.sct$group
immune.combined.sct$condition = gsub("HC", "Healthy", immune.combined.sct$condition)
immune.combined.sct$condition = gsub("^O", "Moderate", immune.combined.sct$condition)
immune.combined.sct$condition = gsub("S/C", "Severe", immune.combined.sct$condition)

immune.combined.sct = AddModuleScore(immune.combined.sct, features = list(lungTRM_noTRM.genes), name = "TRM_", ctrl = 20)
immune.combined.sct = AddModuleScore(immune.combined.sct, features = "CXCR6",assay = "integrated",  name = "CXCR6_", ctrl = 50)
immune.combined.sct = AddModuleScore_UCell(immune.combined.sct, features = list("TRM" = lungTRM_noTRM.genes,"CXCR6" = c("CXCR6")))


hist(immune.combined.sct@assays$integrated@data["CXCR6", ], breaks = 100, main = "CXCR6", xlab = "CXCR6")
abline(v = 1 , col = "red")

CXCR6hi_th = 1
CXCR6hi_cells = which(immune.combined.sct@assays$integrated@data["CXCR6", ] > 1)#0
immune.combined.sct$CXCR6hi = "CXCR6 low"
immune.combined.sct$CXCR6hi[CXCR6hi_cells] = "CXCR6 high"

hist(immune.combined.sct$TRM_1, breaks = 100, main = "TRM", xlab = "TRM")
abline(v = 0.05 , col = "red")

TRMhi_th = 0.05
TRMhi_cells = which(immune.combined.sct$TRM_1 > 0.05)#0.1
immune.combined.sct$TRMhi = "TRM low"
immune.combined.sct$TRMhi[TRMhi_cells] = "TRM high"
table(immune.combined.sct$TRMhi)

my_comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe"),c("Moderate","Severe"))

count_tcells = as.data.frame(table(immune.combined.sct@meta.data[CXCR6hi_cells, "condition"]))
count_tcells2 = as.data.frame(table(immune.combined.sct@meta.data[intersect(TRMhi_cells , CXCR6hi_cells),  "condition"]))

count_t_indv = as.data.frame(table(immune.combined.sct@meta.data[CXCR6hi_cells, "condition"], immune.combined.sct@meta.data[CXCR6hi_cells, "sample_new"])) %>% 
  dplyr::group_by(Var1)%>% dplyr::filter(Freq>0) %>% dplyr::count(Var1)%>% data.frame()

count_t_indv2 = as.data.frame(table(immune.combined.sct@meta.data[intersect(TRMhi_cells , CXCR6hi_cells), "condition"], immune.combined.sct@meta.data[intersect(TRMhi_cells , CXCR6hi_cells), "sample_new"]))%>% 
  dplyr::group_by(Var1)%>% dplyr::filter(Freq>0) %>% dplyr::count(Var1)%>% data.frame()


pdf("manuscript_ready/Fig9_TCells_CXCR6.pdf", width = 7, height = 4)
dittoPlot(immune.combined.sct, var = "CXCR6", assay = "integrated" ,slot = "data", 
          cells.use = CXCR6hi_cells,
          boxplot.lineweight = 0.7,
          vlnplot.lineweight = 0.7,
          plots = c("jitter", "vlnplot", "boxplot"),max = 4.5,
          jitter.size = 0.1, jitter.color = "#00000055",jitter.width = 0.3,
          group.by = "condition")+xlab("")+
  stat_compare_means(label = "p.format", tip.length = 0.03,size = 2.5,method = "wilcox.test",comparisons = my_comparisons)+
  ggtitle("T cells compartment")+scale_fill_brewer()+NoLegend()+
  annotate("text",
           x = count_tcells$Var1,
           y = 3, hjust = -0.1, vjust = -0.5,
           label = paste(count_t_indv$n, "individuals\n", count_tcells$Freq, "cells", sep = " "),
           size = 2, fontface =2)+
  scale_x_discrete(labels=c("Healthy" = "Healthy", "Moderate" = "Moderate\n(infection)",
                            "Severe" = "Severe\n(sepsis)"))+
  
  dittoPlot(immune.combined.sct, var = "CXCR6", assay = "integrated" ,slot = "data", 
            cells.use = intersect(TRMhi_cells , CXCR6hi_cells),
            boxplot.lineweight = 0.7,
            vlnplot.lineweight = 0.7,
            plots = c("jitter", "vlnplot", "boxplot"),max =4.5, 
            jitter.size = 0.1, jitter.color = "#00000055",jitter.width = 0.3,
            group.by = "condition")+xlab("")+
  stat_compare_means(label = "p.format", tip.length = 0.03,size = 2.5,method = "wilcox.test", comparisons = my_comparisons)+
  ggtitle("TRM cells compartment")+scale_fill_brewer()+NoLegend()+
  annotate("text",
           x = count_tcells2$Var1,
           y = 3, hjust = -0.1, vjust = -0.5,
           label = paste(count_t_indv2$n, "individuals\n", count_tcells2$Freq, "cells", sep = " "),
           size = 2, fontface =2)+
  scale_x_discrete(labels=c("Healthy" = "Healthy", "Moderate" = "Moderate\n(infection)",
                            "Severe" = "Severe\n(sepsis)"))

dev.off()


####

cov = readRDS("C:/Projects/CXCR6/covid_nat_med/nCoV.rds")
mac = subset(cov,idents = c('0','1','2','3','4','5','7','8','10','11','12','18','21','22','23','26'))
table(mac$group)
Idents(mac) = "group"
DefaultAssay(mac) = "integrated"

mac$condition = mac$group
mac$condition = gsub("HC", "Healthy", mac$condition)
mac$condition = gsub("^O", "Moderate", mac$condition)
mac$condition = gsub("S/C", "Severe", mac$condition)


my_comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe"),c("Moderate","Severe"))

count_maccells = as.data.frame(table(mac@meta.data[which(mac@assays$integrated@data["CXCL16", ]>0), "condition"]))
count_maccells2 = as.data.frame(table(mac@meta.data[which(mac@assays$integrated@data["CCL7", ] > 0), "condition"]))
count_maccells3 = as.data.frame(table(mac@meta.data[which(mac@assays$integrated@data["CCL2", ] >0), "condition"]))

count_macindv = as.data.frame(table(mac@meta.data[which(mac@assays$integrated@data["CXCL16", ]>0), "condition"], 
                                    mac@meta.data[which(mac@assays$integrated@data["CXCL16", ]>0), "sample_new"])) %>% 
  dplyr::group_by(Var1)%>% dplyr::filter(Freq>0) %>% dplyr::count(Var1)%>% data.frame()

count_macindv2 = as.data.frame(table(mac@meta.data[which(mac@assays$integrated@data["CCL7", ] > 0), "condition"], 
                                     mac@meta.data[which(mac@assays$integrated@data["CCL7", ] > 0), "sample_new"]))%>% 
  dplyr::group_by(Var1)%>% dplyr::filter(Freq>0) %>% dplyr::count(Var1)%>% data.frame()

count_macindv3 = as.data.frame(table(mac@meta.data[which(mac@assays$integrated@data["CCL2", ] >0), "condition"], 
                                     mac@meta.data[which(mac@assays$integrated@data["CCL2", ] >0), "sample_new"]))%>% 
  dplyr::group_by(Var1)%>% dplyr::filter(Freq>0) %>% dplyr::count(Var1)%>% data.frame()


VlnPlot(mac, features = "CCL7", group.by = "condition", assay = "RNA", pt.size = 0)+
  VlnPlot(mac, features = "CCL7", group.by = "condition", assay = "integrated", pt.size = 0)

d1 = dittoPlot(mac, var  = c("CXCL16"), 
               cells.use = which(mac@assays$integrated@data["CXCL16", ]>0),
               assay = "integrated",
               plots = c("jitter", "vlnplot", "boxplot"),max = 5.5,
               jitter.size = 0.1, jitter.color = "#00000044",jitter.width = 0.3,
               group.by = "condition")+
  stat_compare_means(label = "p.format", 
                     tip.length = 0.03,method = "wilcox.test",size = 2,
                     comparisons = my_comparisons)+
  annotate("text",
           x = count_macindv$Var1,
           y = 3, hjust = -0.1, vjust = -1.5,
           label = paste(count_macindv$n, "individuals\n", count_maccells$Freq, "cells", sep = " "),
           size = 2, fontface =2)+
  scale_x_discrete(labels=c("Healthy" = "Healthy", "Moderate" = "Moderate\n(infection)",
                            "Severe" = "Severe\n(sepsis)"))+
  ggtitle("Macrophages")+scale_fill_brewer()+NoLegend()

d2 = dittoPlot(mac, var  = c("CCL7"), 
               cells.use = which(mac@assays$integrated@data["CCL7", ] > 0),
               assay = "integrated", 
               plots = c("jitter", "vlnplot", "boxplot"),max = 8,
               jitter.size = 0.1, jitter.color = "#00000044",jitter.width = 0.3,
               group.by = "condition")+
  stat_compare_means(label = "p.format", tip.length = 0.03,method = "wilcox.test",size = 2,
                     comparisons = my_comparisons)+
  annotate("text",
           x = count_macindv2$Var1,
           y = 4.5, hjust = -0.1, vjust = -1.5,
           label = paste(count_macindv2$n, "individuals\n", count_maccells2$Freq, "cells", sep = " "),
           size = 2, fontface =2)+
  scale_x_discrete(labels=c("Healthy" = "Healthy", "Moderate" = "Moderate\n(infection)",
                            "Severe" = "Severe\n(sepsis)"))+
  ggtitle("Macrophages")+scale_fill_brewer()+NoLegend()

d3 = dittoPlot(mac, var  = c("CCL2"), 
               cells.use = which(mac@assays$integrated@data["CCL2", ] >0),
               assay = "integrated", 
               plots = c("jitter", "vlnplot", "boxplot"),max = 11,
               jitter.size = 0.1, jitter.color = "#00000044",jitter.width = 0.3,
               group.by = "condition")+
  stat_compare_means(label = "p.format", tip.length = 0.03,method = "wilcox.test",size = 2,
                     comparisons = my_comparisons)+
  annotate("text",
           x = count_macindv3$Var1,
           y = 6, hjust = -0.1, vjust = -1.5,
           label = paste(count_macindv3$n, "individuals\n", count_maccells3$Freq, "cells", sep = " "),
           size = 2, fontface =2)+
  scale_x_discrete(labels=c("Healthy" = "Healthy", "Moderate" = "Moderate\n(infection)",
                            "Severe" = "Severe\n(sepsis)"))+
  ggtitle("Macrophages")+scale_fill_brewer()+NoLegend()

pdf("manuscript_ready/Fig_9_mac_genes_vln.pdf", width = 10, height = 4)
gridExtra::grid.arrange(d1,d2,d3, ncol = 3)
dev.off()


##### Plot C
#  C:/Projects/CXCR6/BCG/process_BCG.R


#### Plot E
# C:/Projects/CXCR6/plots_TRM.R line 76
