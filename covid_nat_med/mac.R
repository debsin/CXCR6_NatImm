library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(cowplot)

cov = readRDS("C:/Projects/CXCR6/nCoV.rds")
mac = subset(cov,idents = c('0','1','2','3','4','5','7','8','10','11','12','18','21','22','23','26'))

rm(cov)
gc()
table(mac$group)
Idents(mac) = "group"
DefaultAssay(mac) = "integrated"
dittoPlot(mac, var  = c("CXCL16"), 
          cells.use = which(
            # mac$group %in% c("O", "S/C") &
            mac@assays$RNA@data["CXCL16", ]>0),
          assay = "integrated",
          # jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
          jitter.color = "#00000044",
          group.by = "group", jitter.size = 0.1,boxplot.fill = T,
          plots = c("vlnplot", "jitter", "boxplot"))+ 
  stat_compare_means(label = "p.signif", ref = "HC")+ggtitle("Macrophages", subtitle  = "Uncorrected")


dittoPlot(mac, var  = c("CXCL16"), 
          cells.use = which(
            # mac$group %in% c("O", "S/C") &
            mac@assays$RNA@data["CXCL16", ]>0),
          assay = "integrated",
          # jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
          jitter.color = "#00000044",max = 6,
          group.by = "group", jitter.size = 0.1,boxplot.fill = T,
          plots = c("vlnplot", "jitter", "boxplot"))+ 
  stat_compare_means(comparisons = list(c("HC", "O"),c("O", "S/C"), c("HC", "S/C")),
                     label = "p.signif")+
  ggtitle("Macrophages", subtitle  = "full integration")+ylab(c(0,6))

d2 = dittoPlot(mac, var  = c("CCL7"), 
               cells.use = which(mac@assays$integrated@data["CCL7", ] > 0),
               assay = "integrated", 
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  stat_compare_means(label = "p.signif", ref = "O")+ggtitle("Macrophages")

d3 = dittoPlot(mac, var  = c("CCL2"), 
               cells.use = which(mac@assays$integrated@data["CCL2", ] >0),
               assay = "integrated", 
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  stat_compare_means(label = "p.signif", ref = "O")+ggtitle("Macrophages")

gridExtra::grid.arrange(d1,d2,d3, ncol = 3)

# 
# Idents(mac) = "group"
# mac = subset(mac, cells = grep("HC", mac$group, invert = T))
# 
# 
# saveRDS(mac, "covid_nat_med/mac.rds")
# 
##############mac###############
DefaultAssay(mac) <- "RNA"
Subsemacs.list <- SplitObject(mac, split.by = "sample")
for (i in 1:length(Subsemacs.list)) {
  Subsemacs.list[[i]] <- NormalizeData(Subsemacs.list[[i]], verbose = FALSE)
  Subsemacs.list[[i]] <- FindVariableFeatures(Subsemacs.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}
samples_name = c('C51','C52','C100','GSM3660650','C141','C142','C144','C143','C145','C146','C148','C149','C152')
reference.list <- Subsemacs.list[samples_name]
mac.temp <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50,k.filter = 115)
mac.Integrated <- IntegrateData(anchorset = mac.temp, dims = 1:50)
saveRDS(mac.Integrated, "mac_Integrated.rds")

mac.Integrated = readRDS("covid_nat_med/mac_Integrated.rds")
###first generate data and scaledata in RNA assay
DefaultAssay(mac.Integrated) <- "RNA"
mac.Integrated[['percent.mito']] <- PercentageFeatureSet(mac.Integrated, pattern = "^MT-")
mac.Integrated <- NormalizeData(object = mac.Integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
mac.Integrated <- FindVariableFeatures(object = mac.Integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(mac.Integrated) <- "integrated"

VlnPlot(object = mac.Integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(object = mac.Integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Run the standard workflow for visualization and clustering
mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
mac.Integrated <- RunPCA(mac.Integrated, verbose = FALSE)
#visulaization pca result
mac.Integrated <- ProjectDim(object = mac.Integrated)
ElbowPlot(object = mac.Integrated,ndims = 50)

###cluster
mac.Integrated <- FindNeighbors(object = mac.Integrated, dims = 1:50)
mac.Integrated <- FindClusters(object = mac.Integrated, resolution = 0.8)

###tsne and umap
# mac.Integrated <- RunTSNE(object = mac.Integrated, dims = 1:50)
mac.Integrated <- RunUMAP(mac.Integrated, reduction = "pca", dims = 1:50)
DimPlot(object = mac.Integrated, reduction = 'umap',label = TRUE, group.by = "sample")

mac.Integrated$condition = mac.Integrated$group
mac.Integrated$condition = gsub("HC", "Healthy", mac.Integrated$condition)
mac.Integrated$condition = gsub("^O", "Moderate", mac.Integrated$condition)
mac.Integrated$condition = gsub("S/C", "Severe", mac.Integrated$condition)

saveRDS(mac.Integrated, "covid_nat_med/mac_Integrated.rds")


my_comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe"),c("Moderate","Severe"))
d1 = dittoPlot(mac.Integrated, var  = c("CCL2"),slot = "data",max = 11,
               cells.use = which(mac.Integrated[["integrated"]]@data["CCL2", ]>1),
               plots = c("jitter", "vlnplot", "boxplot"),
               jitter.size = 0.1, jitter.color = "#00000044",
               group.by = "condition", assay = "integrated")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("Macrophages")+scale_fill_hue()+NoLegend()

d2 = dittoPlot(mac.Integrated, var  = c("CCL7"),slot = "data",max = 8,
               cells.use = which(mac.Integrated[["integrated"]]@data["CCL7", ]>0.3),
               plots = c("jitter", "vlnplot", "boxplot"),
               jitter.size = 0.1, jitter.color = "#00000044",
               group.by = "condition", assay = "integrated")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("Macrophages")+scale_fill_hue()+NoLegend()

d3 = dittoPlot(mac.Integrated, var  = c("CXCL16"),slot = "data",max = 5,
               cells.use = which(mac.Integrated[["integrated"]]@data["CXCL16", ]>0.3),
               plots = c("jitter", "vlnplot", "boxplot"),
               jitter.size = 0.1, jitter.color = "#00000044",
               group.by = "condition", assay = "integrated")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("Macrophages")+scale_fill_hue()+NoLegend()

pdf("covid_nat_med/mac_cutoff.pdf", width = 10, height = 6)
plot(density(mac.Integrated[["integrated"]]@data["CCL2", ]), main = "CCL2")
abline(v = 1, col = "red")
plot(density(mac.Integrated[["integrated"]]@data["CCL7", ]), main = "CCL7")
abline(v = 0.3, col = "red")
plot(density(mac.Integrated[["integrated"]]@data["CXCL16", ]), main = "CXCL16")
abline(v = 0.3, col = "red")
dev.off()

pdf("covid_nat_med/mac.pdf", width = 10, height = 6)
gridExtra::grid.arrange(d1,d2,d3, ncol = 3)
dev.off()
DefaultAssay(mac.Integrated) =  "integrated"
FeatureScatter(mac.Integrated, feature1  = "CCL2", feature2 = "CXCL16",slot = "scale.data")
