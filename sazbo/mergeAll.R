
setwd("C:/Projects/CXCR6/covid_szabo/sims-farber/")

library(Seurat)
aww =readRDS("aww.rds")
table(aww$Annotation)
mac = subset(aww, cells = grep("myeloid", aww$Annotation))
rownames(mac[["RNA"]]@counts) = gsub("GRCh38---","", rownames(mac[["RNA"]]@counts))
rownames(mac[["RNA"]]@data) = gsub("GRCh38---","", rownames(mac[["RNA"]]@data))
rownames(mac[["RNA"]]@scale.data) = gsub("GRCh38---","", rownames(mac[["RNA"]]@scale.data))

mac[['percent.mito']] <- PercentageFeatureSet(mac, pattern = "^MT-")
mac = NormalizeData(mac, normalization.method = "LogNormalize", scale.factor = 1e4)
mac = FindVariableFeatures(mac, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
mac <- ScaleData(mac, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
mac = RunPCA(mac, verbose = FALSE)
# ElbowPlot(object = mac,ndims = 50)
mac <- RunUMAP(mac, reduction = "pca", dims = 1:50, seed.use = 100, n.neighbors = 20)
# DimPlot(mac, group.by = "Patient", label = T, label.box = T)
grep("CCL12", rownames(mac@assays$RNA@counts))
mac_genes = c("CXCL16", "CCL7","CCL2" )
DefaultAssay(mac) <- "RNA"
SubseTs.list <- SplitObject(mac, split.by = "Patient")
for (i in 1:length(SubseTs.list)) {
  SubseTs.list[[i]] <- NormalizeData(SubseTs.list[[i]], verbose = FALSE)
  SubseTs.list[[i]] <- FindVariableFeatures(SubseTs.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  SubseTs.list[[i]]@assays$RNA@var.features = union(SubseTs.list[[i]]@assays$RNA@var.features, mac_genes)
  
}

patient_name = c('COV022','COV026','COV027','COV028')
reference.list <- SubseTs.list[patient_name]
mac.temp <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50,k.filter = 140)
mac.int <- IntegrateData(anchorset = mac.temp, dims = 1:50)

saveRDS(mac.int, "mac_int.rds")

mac.int = readRDS("mac_int.rds")
grep("CXCL16", rownames(mac.int@assays$integrated@data))

###first generate data and scaledata in RNA assay
DefaultAssay(mac.int) <- "RNA"
# mac.int[['percent.mito']] <- PercentageFeatureSet(mac.int, pattern = "^MT-")
mac.int <- NormalizeData(object = mac.int, normalization.method = "LogNormalize", scale.factor = 1e4)
mac.int <- FindVariableFeatures(object = mac.int, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
mac.int <- ScaleData(mac.int, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

VlnPlot(object = mac.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 2)
FeatureScatter(object = mac.int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

##change to integrated assay
DefaultAssay(mac.int) <- "integrated"

# Run the standard workflow for visualization and clustering
mac.int <- FindVariableFeatures(object = mac.int, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
mac.int <- ScaleData(mac.int, verbose = TRUE, vars.to.regress = c("nCount_RNA", "percent.mito"),
                     features = union(mac.int@assays$integrated@var.features, mac_genes))
mac.int <- RunPCA(mac.int, verbose = FALSE)
#visulaization pca result
mac.int <- ProjectDim(object = mac.int)
ElbowPlot(object = mac.int,ndims = 50)

mac.int <- RunUMAP(mac.int, reduction = "pca", dims = 1:50, seed.use = 100, n.neighbors = 20)
DimPlot(mac.int,  group.by = "Patient", label = T, label.box = T, order = T, shuffle = T, pt.size = 0.1)
mac.int@assays$integrated@counts = mac.int@assays$integrated@data

mac.int = AddModuleScore(mac.int, assay = "integrated",ctrl = 50,
                           features = list(c("CXCL16")), name = "CXCL16_")
mac.int = AddModuleScore(mac.int, assay = "integrated",ctrl = 50,
                         features = list(c("CCL7")), name = "CCL7_")
mac.int = AddModuleScore(mac.int, assay = "integrated",ctrl = 50,
                         features = list(c("CCL2")), name = "CCL2_")

library(ggpubr)
library(dittoSeq)
DefaultAssay(mac.int) = "integrated"
d1 = dittoPlot(mac.int, var =  "CXCL16_1",
               assay = "integrated", 
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "Timepoint", jitter.size = 0.1,plots = c("boxplot", "jitter"))+ylab("CXCL16")+
  stat_compare_means(label = "p.signif", ref = "IntubationDay1")+ggtitle("Myeloid")+NoLegend()
d2 = dittoPlot(mac.int, var =  "CCL7_1",
               assay = "integrated", 
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "Timepoint", jitter.size = 0.1,plots = c("boxplot", "jitter"))+ylab("CCL7")+
  stat_compare_means(label = "p.signif", ref = "IntubationDay1")+ggtitle("Myeloid")+NoLegend()

d3 = dittoPlot(mac.int, var =  "CCL2_1",
          assay = "integrated", 
          jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
          group.by = "Timepoint", jitter.size = 0.1,plots = c("boxplot", "jitter"))+ylab("CCL2")+
  stat_compare_means(label = "p.signif", ref = "IntubationDay1")+ggtitle("Myeloid")+NoLegend()

gridExtra::grid.arrange(d1,d2,d3, ncol = 2)
