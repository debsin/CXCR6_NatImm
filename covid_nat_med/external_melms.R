library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)

gene.path = "C:/Projects/external/jem_20190249_tables1/"

lungTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T)
lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[lungTRM_noTRM$`log2 fold change`>2 & lungTRM_noTRM$`P adj`<0.01]
# lungTRM_noTRM.genes = setdiff(lungTRM_noTRM.genes, "CXCR6")

counts = Read10X("C:/Projects/external/melms/", gene.column = 1, cell.column = 1)

meta = read.csv("C:/Projects/external/melms/lung_metaData.txt", header = T, row.names = 1, sep = "\t")
meta = meta[-1,]

counts = counts[,rownames(meta)]
cov = CreateSeuratObject(counts = counts,assay = "RNA", meta.data = meta,min.cells = 0,min.features = 0,project = "melms")
cov[["RNA"]]@data = cov[["RNA"]]@counts

table(cov$group)

cov = FindVariableFeatures(cov)
VariableFeatures(cov) = union(VariableFeatures(cov), lungTRM_noTRM.genes)
cov = ScaleData(cov)
cov = RunPCA(cov)
cov = RunUMAP(cov, dims = 1:50)


# moderate covid patients (=infection) and severe COVID patients (=sepsis)
# Can you compare the CXCR6 signature in the T cell compartment between these 2 groups ? 
# and of CXCR6 TRM signature ? 
# Ideally, we want to show that sepsis patients have higher CXCR6 expression in T cells (TRM cells) than infection patients .
# and if you can look at CXCL16, Ccl2, Ccl7 and ccl12 in the macrophage compartment ++

table(cov$cell_type_intermediate, cov$cell_type_main)

table(cov$donor_id)
DimPlot(cov, group.by = "cell_type_intermediate", label = T, label.box = T)



cov.T = subset(cov, cells = grep("T cells", cov$cell_type_main))


DimPlot(cov.T, group.by = "cell_type_intermediate", label = T, label.box = T)
cov.T= FindVariableFeatures(cov.T)
cov.T = ScaleData(cov.T, features = union(VariableFeatures(cov.T), lungTRM_noTRM.genes))
cov.T = RunPCA(cov.T)
ElbowPlot(cov.T)
cov.T = RunUMAP(cov.T, dims = 1:15)
DimPlot(cov.T, group.by = "cell_type_intermediate", label = T, label.box = T)
DimPlot(cov.T, group.by = "donor_id", label = T, label.box = T)

NKT_T = cov.T



## Here are the plots
NKT_T = AddModuleScore(NKT_T, features = list("TRM" = lungTRM_noTRM.genes), name = "TRM_")

NKT_T = AddModuleScore(NKT_T, features = list("CXCR6" = c("CXCR6")), name = "CXCR6_")



# plot(density(NKT_T$TRM_1), main = "TRM Signature")
hist(NKT_T$TRM_1,breaks = 500, main = "TRM Signature", xlab = "signature score")
abline(v = quantile(NKT_T$TRM_1, 0.9), col = "red", lwd = 2)
text(0.07, 700, "quantile @ 0.9")

# plot(density(NKT_T$CXCR6_1), main = "CXCR6 Signature")
hist(NKT_T$CXCR6_1,breaks = 1000, main = "CXCR6 Signature",  xlab = "signature score")
abline(v = quantile(NKT_T$CXCR6_1, 0.5), col = "red", lwd = 2)
text(0.6, 200, "quantile @ 0.75")

FeaturePlot(NKT_T, features = "TRM_1", 
            min.cutoff = quantile(NKT_T$TRM_1, 0.9), order = T, ncol = 3)+ggtitle("TRM Signature")+
  FeaturePlot(NKT_T, features = "CXCR6_1", 
              min.cutoff = quantile(NKT_T$CXCR6_1,0.5), order = T)+ggtitle("CXCR6 Signature")+
  FeatureScatter(object = NKT_T, feature1 = "TRM_1", feature2 = "CXCR6_1", group.by = "group", shuffle = T)+
  xlab("TRM Signature Score")+ylab("CXCR6")

Trm_rich = subset(NKT_T, 
                  cells = which(NKT_T$TRM_1 > quantile(NKT_T$TRM_1, 0.5)))
dittoPlot(Trm_rich, var = "CXCR6_1", assay = "RNA" ,split.ncol = 5,slot = "data",
          cells.use = which(Trm_rich$CXCR6_1 > quantile(Trm_rich$CXCR6_1 ,0.75) ),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("TRM Rich cells")+
  stat_compare_means(label = "p.signif", ref = "O")


dittoPlot(NKT_T, var = "CXCR6", 
          # cells.use = which(NKT_T@assays$RNA@scale.data["CXCR6",] > 0),
          split.ncol = 5,slot = "data",
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("T cells")+
  stat_compare_means(label = "p.signif", , ref = "Control")


dittoPlot(NKT_T, var = "CXCR6", 
          cells.use = which(NKT_T@assays$RNA@data["CXCR6",] > 0),
          split.ncol = 5,slot = "data",
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("T cells")+
  stat_compare_means(label = "p.signif", , ref = "Control")

dittoPlot(NKT_T, var = "CXCR6_1", split.ncol = 5,slot = "data",
          cells.use = intersect(which(NKT_T$CXCR6_1 > quantile(NKT_T$CXCR6_1 ,0)),
            which(NKT_T$TRM_1 > quantile(NKT_T$TRM_1 ,0.5) )),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("CXCR6 TRM cells")+
  stat_compare_means(label = "p.signif", ref = "Control")

table(NKT_T@assays$RNA@scale.data["CXCR6",] > 0)


### Macrophages

mac = subset(cov, cells = grep("Macrophages", cov$cell_type_intermediate))
mac[['percent.mito']] <- PercentageFeatureSet(mac, pattern = "^MT-")
mac = FindVariableFeatures(mac, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
mac <- ScaleData(mac, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
mac = RunPCA(mac, verbose = FALSE)
# ElbowPlot(object = mac,ndims = 50)
mac <- RunUMAP(mac, reduction = "pca", dims = 1:50, seed.use = 100)
DimPlot(mac, group.by = "cell_type_fine", label = T, label.box = T)
grep("CCL12", rownames(mac@assays$RNA@counts))
mac_genes = c("CXCL16", "CCL7","CCL2" )

grep("CXCL16", rownames(mac@assays$RNA@data))

mac = AddModuleScore(mac, assay = "RNA",ctrl = 50,
                         features = list(c("CXCL16")), name = "CXCL16_")
mac = AddModuleScore(mac, assay = "RNA",ctrl = 50,
                         features = list(c("CCL7")), name = "CCL7_")
mac = AddModuleScore(mac, assay = "RNA",ctrl = 50,
                         features = list(c("CCL2")), name = "CCL2_")

library(ggpubr)
library(dittoSeq)
DefaultAssay(mac) = "RNA"
d1 = dittoPlot(mac, var =  "CXCL16_1",
               assay = "RNA", 
               # cells.use = mac$CXCL16_1 >0,
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+ylab("CXCL16")+
  stat_compare_means(label = "p.signif", ref = "Control")+ggtitle("Myeloid")+NoLegend()

d2 = dittoPlot(mac, var =  "CCL7_1",
               assay = "RNA", 
               # cells.use = mac$CCL7_1 >0,
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+ylab("CCL7")+
  stat_compare_means(label = "p.signif", ref = "Control")+ggtitle("Myeloid")+NoLegend()

d3 = dittoPlot(mac, var =  "CCL2_1",
               assay = "RNA", 
               # cells.use = mac$CCL2_1 >0,
               jitter.width = 0.3,jitter.color = "#00000055",boxplot.width = 0.9,
               group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+ylab("CCL2")+
  stat_compare_means(label = "p.signif", ref = "Control")+ggtitle("Myeloid")+NoLegend()

gridExtra::grid.arrange(d1,d2,d3, ncol = 3)

