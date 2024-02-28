library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(cowplot)
library(dplyr)

cov = readRDS("C:/Projects/CXCR6/nCoV.rds")
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


covid_MR1 = cov@meta.data
covid_MR1$MR1 = cov@assays$RNA@data["MR1", ]
covid_MR1$MR1_raw = cov@assays$RNA@counts["MR1", ]

colnames(covid_MR1) = gsub("condition", "Condition", colnames(covid_MR1))

table(covid_MR1$sample, covid_MR1$sample_new)


DefaultAssay(cov)
cov = NormalizeData(cov)
cov = FindVariableFeatures(cov)
cov@assays$RNA@var.features = c(cov@assays$RNA@var.features , "MR1")
cov = ScaleData(cov, vars.to.regress = c("sample", "nCount_RNA", "percent.mito"))
cov = RunPCA(cov)
cov = RunUMAP(cov, dim = 1:15)
DimPlot(cov, group.by = "sample")+
DimPlot(cov, group.by = "anno")


g1 = ggplot(covid_MR1, aes(y =  MR1, x = Condition, fill = Condition))+
  geom_violin(scale = "count", trim = F,  show.legend = F)+
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
  theme_classic2()+
  scale_fill_brewer(type = "qual", palette = "Set1")+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

head(covid_MR1)

ggplot(covid_MR1, aes(y =  MR1, x = Condition, fill = sample_new))+
  geom_violin(scale = "count", trim = F,  show.legend = F)+
  geom_boxplot(width = 0.8, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
  facet_wrap(~anno, ncol = 5)+
  theme_classic2()+
  scale_fill_brewer(type = "qual", palette = "Set1")+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

covid_MR1_clean = covid_MR1
covid_MR1_clean$cell_type = covid_MR1_clean$anno


# making the plot
g2 = 
  ggplot(covid_MR1_clean[covid_MR1_clean$Condition=="Healthy",] , 
         aes(y =  MR1, x = cell_type, fill = cell_type, col = cell_type))+
  geom_violin(scale = "count", trim = F,  show.legend = F)+
  geom_boxplot(width = 0.3, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
  geom_jitter(width = 0.1, size = 0.2, show.legend = F)+
  theme_classic2()+
  xlab("Cells from Healthy individuals")+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


library(dplyr)
covid_means <- covid_MR1_clean[covid_MR1_clean$Condition=="Healthy",] %>% 
  group_by( cell_type) %>% 
  summarise( count = n(), mean = mean(MR1), prop = sum(MR1_raw>0)/count)
head(covid_means)

g3 = ggplot(covid_means, aes(x = cell_type, y = prop, fill = cell_type, label = count)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = F)+
  geom_label(show.legend = F)+
  ylab("Proportion of MR1>0 expressing cells")+
  xlab("Count of annotated cells within Healthy individuals")+
  # ggtitle(label = "", subtitle = "")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

table(covid_MR1$Condition)

g4 = ggplot(covid_MR1[covid_MR1$Condition %in% c("Healthy", "Moderate", "Severe") & 
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
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


png("C:/Projects/CXCR6/covid_nat_med/MR1_mac_liao.png", width = 1600, height = 1400, res = 150, units = "px")
gridExtra::grid.arrange(g1,g2,g3,g4, ncol = 2)
dev.off()


#### 
 
############## mac Integration -----
inp = readRDS("covid_nat_med/mac.rds")
mac = CreateSeuratObject(counts = inp[["RNA"]]@counts, meta.data = inp@meta.data, project = inp@project.name)
rm(inp)
gc()
DefaultAssay(mac) <- "RNA"
Subsemacs.list <- SplitObject(mac, split.by = "sample")
anchor.features = list()
for (i in 1:length(Subsemacs.list)) {
  Subsemacs.list[[i]] <- NormalizeData(Subsemacs.list[[i]], verbose = FALSE)
  Subsemacs.list[[i]] <- FindVariableFeatures(Subsemacs.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  Subsemacs.list[[i]]@assays$RNA@var.features <-  c(Subsemacs.list[[i]]@assays$RNA@var.features, "MR1")
  
  anchor.features[[i]] = Subsemacs.list[[i]]@assays$RNA@var.features
  Subsemacs.list[[i]] <- ScaleData(Subsemacs.list[[i]])
  Subsemacs.list[[i]] <- RunPCA(Subsemacs.list[[i]], verbose = F)
}
# samples_name = c('C51','C52','C100','GSM3660650','C141','C142','C144','C143','C145','C146','C148','C149','C152')
# reference.list <- Subsemacs.list[samples_name]
mac.temp <- FindIntegrationAnchors(object.list = Subsemacs.list,reduction = "rpca", 
                                   anchor.features = unique(unlist(anchor.features)),
                                     dims = 1:50,k.filter = 115)


# saveRDS(mac.temp, "covid_nat_med/mac_temp.rds")
# mac.temp = readRDS("mac_temp.rds")
# mac.Integrated <- IntegrateData(anchorset = mac.temp, dims = 1:50)
# saveRDS(mac.temp, "covid_nat_med/mac_temp.rds")
mac.Integrated = readRDS("mac_int_MR1.rds")
###first generate data and scaledata in RNA assay
# DefaultAssay(mac.Integrated) <- "RNA"
# mac.Integrated[['percent.mito']] <- PercentageFeatureSet(mac.Integrated, pattern = "^MT-")
# mac.Integrated <- NormalizeData(object = mac.Integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
# mac.Integrated <- FindVariableFeatures(object = mac.Integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
# mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
#
DefaultAssay(mac.Integrated) <- "integrated"
mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE)
mac.Integrated <- RunPCA(mac.Integrated, npcs = 30, verbose = FALSE)
mac.Integrated <- RunUMAP(mac.Integrated, reduction = "pca", dims = 1:30)