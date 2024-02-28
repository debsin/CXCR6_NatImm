library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(cowplot)
library(dplyr)

cov = readRDS("C:/Projects/CXCR6/covid_nat_med/nCoV.rds")
anno = read.csv("C:/Projects/CXCR6/covid_nat_med/all.cell.annotation.meta.txt", sep = "\t", header = T)

cov@assays$integrated@data[1:5, 1:5]
cov = subset(cov, cells = match(anno$ID, cov$ID))
cov$anno = anno$celltype
table(cov$sample, cov$disease)
table(cov$seurat_clusters)

DimPlot(cov, group.by = "anno", label = T, label.box = T)
table(cov$cluster, cov$anno)

t.anno = read.csv("C:/Projects/CXCR6/covid_nat_med/NKT.cell.annotation.meta.txt", sep = "\t", header = T)

cov$anno.T_l2 = cov$anno
cov$anno.T_l2[t.anno$ID] = t.anno$celltype
table(cov$anno.T_l2, cov$anno)

table(cov$sample, cov$group)

cov$condition = cov$group
cov$condition = gsub("HC", "Healthy", cov$condition)
cov$condition = gsub("^O", "Moderate", cov$condition)
cov$condition = gsub("S/C", "Severe", cov$condition)

table(cov$anno.T_l2, cov$condition)


genes = c("MR1","HLA-DRA", "HLA-A", "CD3E")


prop_plot = list()
cond_plot= list()
for(gene in genes){
  covid_MR1 = cov@meta.data
  covid_MR1$MR1 = cov@assays$RNA@data[gene,  ]
  covid_MR1$MR1_raw = cov@assays$RNA@counts[gene, ]
  
  colnames(covid_MR1) = gsub("condition", "Condition", colnames(covid_MR1))
  
  table(covid_MR1$sample, covid_MR1$sample_new)
  
  covid_MR1_clean = covid_MR1
  covid_MR1_clean$cell_type = covid_MR1_clean$anno
  
  covid_means <- covid_MR1_clean[covid_MR1_clean$Condition=="Healthy",] %>% 
    group_by( cell_type) %>% 
    summarise( count = n(), mean = mean(MR1), prop = sum(MR1_raw>0)/count)
  head(covid_means)
  
  prop_plot[[gene]] = ggplot(covid_means, aes(x = cell_type, y = prop, fill = cell_type, label = count)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = F, width = 0.7)+
    geom_label(show.legend = F)+
    ylab(paste0("Proportion of ",  gene, ">0 expressing cells"))+
    xlab("Count of annotated cells within Healthy individuals")+
    # ggtitle(label = "", subtitle = "")+
    theme_classic2()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

prop_plot[[4]]

library(gridExtra)
pdf("C:/Projects/CXCR6/covid_nat_med/hamish/prop_4_genes.pdf", width = 22, height = 6)
do.call("grid.arrange", c(prop_plot, ncol = 4)) 
dev.off()


cond_plot= list()
for(gene in genes[2:4]){
  covid_MR1 = cov@meta.data
  covid_MR1$MR1 = cov@assays$integrated@data[gene,  ]
  covid_MR1$MR1_raw = cov@assays$RNA@counts[gene, ]
  
  colnames(covid_MR1) = gsub("condition", "Condition", colnames(covid_MR1))
  

  
  cond_plot[[gene]] = ggplot(covid_MR1[covid_MR1$Condition %in% c("Healthy", "Moderate", "Severe") & 
                                         # substr(covid_MR1$patient_id, start = 1, stop = 2)=="MH" &
                                         covid_MR1$anno=="Macrophages" & 
                                         covid_MR1$MR1_raw>0,] , 
                             aes(y =  MR1, x = Condition, fill = Condition))+
    geom_violin(scale = "count", trim = F,  show.legend = F)+
    geom_jitter(size = 0.5,  width = 0.1,show.legend = F)+
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
    stat_compare_means(comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe") ), label = "p.format")+
    theme_classic2()+
    xlab("Macrophages Only")+
    scale_y_log10()+ylab(paste(gene, "(log scale)"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}



library(gridExtra)


pdf("C:/Projects/CXCR6/covid_nat_med/mac/MAC_genes_full_intg.pdf", width = 15, height = 6)
do.call("grid.arrange", c(cond_plot, ncol = 4)) 
dev.off()

## sTART maC ----
mac.Integrated = readRDS("C:/Projects/CXCR6/covid_nat_med/mac/2-Macrophage.rds")
mac.Integrated$condition = mac.Integrated$group
mac.Integrated$condition = gsub("HC", "Healthy", mac.Integrated$condition)
mac.Integrated$condition = gsub("^O", "Moderate", mac.Integrated$condition)
mac.Integrated$condition = gsub("S/C", "Severe", mac.Integrated$condition)
###first generate data and scaledata in RNA assay
# DefaultAssay(mac.Integrated) <- "RNA"
# mac.Integrated[['percent.mito']] <- PercentageFeatureSet(mac.Integrated, pattern = "^MT-")
# mac.Integrated <- NormalizeData(object = mac.Integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
# mac.Integrated <- FindVariableFeatures(object = mac.Integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
# mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
#
DefaultAssay(mac.Integrated) <- "integrated"
mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE)
# mac.Integrated <- RunPCA(mac.Integrated, npcs = 30, verbose = FALSE)
# mac.Integrated <- RunUMAP(mac.Integrated, reduction = "pca", dims = 1:30)
# 
genes = c("MR1","HLA-DRA", "HLA-A", "CD3E")

cond_plot= list()
for(gene in genes){
  covid_MR1 = mac.Integrated@meta.data
  covid_MR1$MR1 = mac.Integrated@assays$integrated@scale.data[gene,  ]
  covid_MR1$MR1_raw = mac.Integrated@assays$RNA@counts[gene, ]
  
  colnames(covid_MR1) = gsub("condition", "Condition", colnames(covid_MR1))
  
  table(covid_MR1$sample, covid_MR1$sample_new)
  
  cond_plot[[gene]] = ggplot(covid_MR1[covid_MR1$Condition %in% c("Healthy", "Moderate", "Severe"),] , 
                             aes(y =  MR1, x = Condition, fill = Condition))+
    geom_violin(scale = "count", trim = F,  show.legend = F, size = 0.1)+
    geom_jitter(size = 0.5,  width = 0.1,show.legend = F)+
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.size = 0.5,  show.legend = F)+
    stat_compare_means(comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe") ), label = "p.format")+
    theme_classic2()+
    xlab("Macrophages Only")+
    # scale_y_log10()+ylab(paste(gene, "(log scale)"))+
    ylab(paste(gene))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


library(gridExtra)
pdf("MAC_genes_intg.pdf", width = 15, height = 6)
do.call("grid.arrange", c(cond_plot, ncol = 4)) 
dev.off()


for(gene in genes){
  covid_MR1 = mac.Integrated@meta.data
  covid_MR1[[gene]] =  mac.Integrated@assays$integrated@scale.data[gene,  ]
  colnames(covid_MR1) = gsub("condition", "Condition", colnames(covid_MR1))
  meta_tab = covid_MR1[covid_MR1$Condition %in% c("Healthy", "Moderate", "Severe"),c("Condition", gene, "sample")]
  meta_tab = cbind(index = seq(nrow(meta_tab))-1, meta_tab)
  rownames(meta_tab) = NULL
  write.csv(x = meta_tab,row.names = F,col.names = T,
                   file = file.path("C:/Projects/CXCR6/covid_nat_med/mac/sv", 
                                    paste0(gene, ".csv", collapse = "")))
}



