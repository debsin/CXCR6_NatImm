library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)

gene.path = "C:/Projects/external/jem_20190249_tables1/"

lungTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T)
lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[lungTRM_noTRM$`log2 fold change`>2 & lungTRM_noTRM$`P adj`<0.01]
# lungTRM_noTRM.genes = setdiff(lungTRM_noTRM.genes, "CXCR6")


cov = readRDS("nCoV.rds")
anno = read.csv("covid_nat_med/all.cell.annotation.meta.txt", sep = "\t", header = T)

cov = subset(cov, cells = match(anno$ID, cov$ID))
cov$anno = anno$celltype
table(cov$sample, cov$disease)
# moderate covid patients (=infection) and severe COVID patients (=sepsis)
# Can you compare the CXCR6 signature in the T cell compartment between these 2 groups ? 
# and of CXCR6 TRM signature ? 
# Ideally, we want to show that sepsis patients have higher CXCR6 expression in T cells (TRM cells) than infection patients .
# and if you can look at CXCL16, Ccl2, Ccl7 and ccl12 in the macrophage compartment ++
table(cov$seurat_clusters)

DimPlot(cov, group.by = "anno", label = T, label.box = T)
table(cov$cluster, cov$anno)

t.anno = read.csv("covid_nat_med/NKT.cell.annotation.meta.txt", sep = "\t", header = T)

cov$anno.T_l2 = "non"
cov$anno.T_l2[t.anno$ID] = t.anno$celltype

table(cov$sample, cov$group)

# keep only clean T cells
cov.T = subset(cov, cells = grep("T", cov$anno.T_l2))

DimPlot(cov.T, group.by = "anno.T_l2", label = T, label.box = T)

DimPlot(cov.T, group.by = "sample", label = T, label.box = T)

table(cov.T$anno.T_l2)
##############T###############
DefaultAssay(cov.T) <- "RNA"
SubseNKTs.list <- SplitObject(cov.T, split.by = "sample")
for (i in 1:length(SubseNKTs.list)) {
  SubseNKTs.list[[i]] <- NormalizeData(SubseNKTs.list[[i]], verbose = FALSE)
  SubseNKTs.list[[i]] <- FindVariableFeatures(SubseNKTs.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  SubseNKTs.list[[i]]@assays$RNA@var.features = union(SubseNKTs.list[[i]]@assays$RNA@var.features, lungTRM_noTRM.genes)
}

table(cov.T$sample, cov.T$anno.T_l2)
samples_name = c('C51','C100','GSM3660650','C141','C142','C144','C143','C145','C148','C149','C152')
reference.list <- SubseNKTs.list[samples_name]
NKT.temp <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50,k.filter = 140)
NKT.Integrated_clean <- IntegrateData(anchorset = NKT.temp, dims = 1:50)

###first generate data and scaledata in RNA assay
DefaultAssay(NKT.Integrated_clean) <- "RNA"
NKT.Integrated_clean[['percent.mito']] <- PercentageFeatureSet(NKT.Integrated_clean, pattern = "^MT-")
NKT.Integrated_clean <- NormalizeData(object = NKT.Integrated_clean, normalization.method = "LogNormalize", scale.factor = 1e4)
NKT.Integrated_clean <- FindVariableFeatures(object = NKT.Integrated_clean, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
NKT.Integrated_clean@assays$RNA@var.features = union(NKT.Integrated_clean@assays$RNA@var.features, lungTRM_noTRM.genes)

NKT.Integrated_clean <- ScaleData(NKT.Integrated_clean, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(NKT.Integrated_clean) <- "integrated"
VlnPlot(object = NKT.Integrated_clean, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(object = NKT.Integrated_clean, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

NKT.Integrated_clean = subset(NKT.Integrated_clean, cells = grep("T", NKT.Integrated_clean$anno.T_l2))
# Run the standard workflow for visualization and clustering
NKT.Integrated_clean <- ScaleData(NKT.Integrated_clean, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
NKT.Integrated_clean <- RunPCA(NKT.Integrated_clean, verbose = FALSE)
#visulaization pca result
NKT.Integrated_clean <- ProjectDim(object = NKT.Integrated_clean)
ElbowPlot(object = NKT.Integrated_clean,ndims = 50)

NKT.Integrated_clean <- RunUMAP(NKT.Integrated_clean, reduction = "pca", dims = 1:50, seed.use = 100, n.neighbors = 20)

pdf("covid_nat_med/t_cell_integrated_umap.pdf", width = 12, height = 5)
DimPlot(NKT.Integrated_clean, group.by = "anno.T_l2", label = T, label.box = T)+
  DimPlot(NKT.Integrated_clean, group.by = "sample", label = T, label.box = T)
dev.off()

DefaultAssay(NKT.Integrated_clean) = "integrated"
saveRDS(NKT.Integrated_clean, "covid_nat_med/T_integrated.rds")

# gridExtra::grid.arrange(p1,p2, p3, layout_matrix=rbind(c(1,2),c(3,3)))

NKT.Integrated_clean@assays$integrated@counts = NKT.Integrated_clean@assays$integrated@data


# NKT.Integrated_clean = AddModuleScore(NKT.Integrated_clean, assay = "integrated", name = "TRM",
#                                           features = lungTRM_noTRM.genes)
# NKT.Integrated_clean = AddModuleScore_UCell(NKT.Integrated_clean, assay = "integrated",
#                                             features = list("TRM" = lungTRM_noTRM.genes))
# 
# NKT.Integrated_clean = AddModuleScore_UCell(NKT.Integrated_clean, assay = "integrated",
#                                             features = list("CXCR6" = c("CXCR6")))
NKT.Integrated_clean@meta.data = droplevels(NKT.Integrated_clean@meta.data)


density(NKT.Integrated_clean$TRM_UCell)
plot(density(NKT.Integrated_clean$TRM_UCell))

NKT.Integrated_clean$condition = NKT.Integrated_clean$group
NKT.Integrated_clean$condition = gsub("HC", "Healthy", NKT.Integrated_clean$condition)
NKT.Integrated_clean$condition = gsub("^O", "Moderate", NKT.Integrated_clean$condition)
NKT.Integrated_clean$condition = gsub("S/C", "Severe", NKT.Integrated_clean$condition)


## Here are the plots
NKT.Integrated_clean = AddModuleScore(NKT.Integrated_clean, assay = "integrated",ctrl = 50,
                                      features = list("TRM" = lungTRM_noTRM.genes), name = "TRM_")

NKT.Integrated_clean = AddModuleScore(NKT.Integrated_clean, assay = "integrated",ctrl = 50,
                                      features = list("CXCR6" = c("CXCR6")), name = "CXCR6_")


DimPlot(NKT.Integrated_clean, group.by = "anno.T_l2", label = T, label.box = T)

# plot(density(NKT.Integrated_clean$TRM_1), main = "TRM Signature")
hist(NKT.Integrated_clean$TRM_1,breaks = 100, main = "TRM Signature", xlab = "signature score")
abline(v = quantile(NKT.Integrated_clean$TRM_1, 0.5), col = "red", lwd = 2)
text(0.005, 215, "median")

# plot(density(NKT.Integrated_clean$CXCR6_1), main = "CXCR6 Signature")
hist(NKT.Integrated_clean$CXCR6_1,breaks = 100, main = "CXCR6 Signature",  xlab = "signature score")
abline(v = quantile(NKT.Integrated_clean$CXCR6_1, 0.75), col = "red", lwd = 2)
text(0.6, 200, "quantile @ 0.75")

FeaturePlot(NKT.Integrated_clean, features = "TRM_1", slot = "data",
            min.cutoff = quantile(NKT.Integrated_clean$TRM_1, 0.5), order = T, ncol = 3)+ggtitle("TRM Signature")+
  FeaturePlot(NKT.Integrated_clean, features = "CXCR6_1", 
              min.cutoff = quantile(NKT.Integrated_clean$CXCR6_1, 0.75), order = T)+ggtitle("CXCR6 Signature")

cxcr6_pos = which(NKT.Integrated_clean[["integrated"]]@data["CXCR6", ]>0)
pdf("covid_nat_med/t_cells_trm_CXCR6.pdf", width = 15, height = 4.5)
FeatureScatter(object = NKT.Integrated_clean, slot = "data",
               cells = intersect(cxcr6_pos, grep("Healthy", NKT.Integrated_clean$condition, invert = F)),
               feature1 = "TRM_1", feature2 = "CXCR6", group.by = "condition", shuffle = T)+
  xlab("TRM Signature Score")+ylab("CXCR6")+
  FeatureScatter(object = NKT.Integrated_clean, slot = "data",
                 cells = intersect(cxcr6_pos, grep("Healthy", NKT.Integrated_clean$condition, invert = T)), 
                 feature1 = "TRM_1", feature2 = "CXCR6", group.by = "condition", shuffle = T)+
  xlab("TRM Signature Score")+ylab("CXCR6")+
  FeatureScatter(object = NKT.Integrated_clean,  slot = "data",
                 cells = cxcr6_pos, 
                 feature1 = "TRM_1", feature2 = "CXCR6", group.by = "condition", shuffle = T)+
  xlab("TRM Signature Score")+ylab("CXCR6")
dev.off()
Trm_rich = subset(NKT.Integrated_clean, 
                  cells = which(NKT.Integrated_clean$TRM_1 > quantile(NKT.Integrated_clean$TRM_1, 0.5)))
dittoPlot(Trm_rich, var = "CXCR6_1", assay = "integrated" ,split.ncol = 5,slot = "data",
          cells.use = which(Trm_rich$CXCR6_1 > quantile(Trm_rich$CXCR6_1 ,0.75) ),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("TRM Rich cells")+
  stat_compare_means(label = "p.signif", ref = "O")


dittoPlot(NKT.Integrated_clean, var = "CXCR6_1", assay = "integrated" ,split.ncol = 5,slot = "data",
          cells.use = which(NKT.Integrated_clean$CXCR6_1 > quantile(NKT.Integrated_clean$CXCR6_1 ,0.75) ),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("T cells")+
  stat_compare_means(label = "p.signif", , ref = "O")

dittoPlot(NKT.Integrated_clean, var = "CXCR6_1", assay = "integrated" ,split.ncol = 5,slot = "data",
          cells.use = intersect(which(NKT.Integrated_clean$CXCR6_1 > quantile(NKT.Integrated_clean$CXCR6_1 ,0.75)),
                                which(NKT.Integrated_clean$TRM_1 > quantile(NKT.Integrated_clean$TRM_1 ,0.5) )),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "group", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("CXCR6 TRM cells")+
  stat_compare_means(label = "p.signif", , ref = "O")

table(NKT.Integrated_clean$anno.T_l2)









# t1 = dittoPlot(NKT.Integrated_clean, var = "TRM1", assay = "integrated" ,split.ncol = 5,
#                cells.use = grep("Doublets|NK|Uncertain", NKT.Integrated_clean$anno.T_l2, invert = T),
#                group.by = "group", split.by ="anno.T_l2", jitter.size = 0.1)+
#   geom_boxplot(fill = "white", width = 0.3, outlier.shape = NA)+scale_fill_brewer(palette = "Set1")

t1 = dittoPlot(NKT.Integrated_clean, var = "TRM_UCell", assay = "integrated" ,split.ncol = 5,
               cells.use = grep("Doublets|NK|Uncertain", NKT.Integrated_clean$anno.T_l2, invert = T),
               group.by = "group", split.by ="anno.T_l2", jitter.size = 0.1)+
  geom_boxplot(fill = "white", width = 0.3, outlier.shape = NA)+scale_fill_brewer(palette = "Set1")

t2 = dittoPlot(NKT.Integrated_clean, var = "CXCR6_UCell", assay = "integrated" ,split.ncol = 5,
               cells.use = grep("Doublets|NK|Uncertain", NKT.Integrated_clean$anno.T_l2, invert = T),
               group.by = "group", split.by ="anno.T_l2", jitter.size = 0.1)+
  geom_boxplot(fill = "white", width = 0.3, outlier.shape = NA)+scale_fill_brewer(palette = "Set1")


t3 = dittoPlot(NKT.Integrated_clean, var = "CXCR6", assay = "integrated" ,split.ncol = 5,
               cells.use = grep("Doublets|NK|Uncertain", NKT.Integrated_clean$anno.T_l2, invert = T),
               group.by = "group", split.by ="anno.T_l2", jitter.size = 0.1)+
  geom_boxplot(fill = "white", width = 0.3, outlier.shape = NA)+scale_fill_brewer(palette = "Set1")

gridExtra::grid.arrange(t1,t2, t3, ncol = 1)


DimPlot(NKT.Integrated_clean, group.by = "anno.T_l2", order = T, pt.size = 0.5, label.box = T, label = T)
FeaturePlot(NKT.Integrated_clean, features = "CXCR6" ,
            min.cutoff = quantile(NKT.Integrated_clean$CXCR6_UCell, 0.75), pt.size = 0.6, order = T)
FeaturePlot(NKT.Integrated_clean, features = "TRM_UCell" ,
            min.cutoff = quantile(NKT.Integrated_clean$TRM_UCell, 0.75), pt.size = 0.6, order = T)

dittoPlot(NKT.Integrated_clean, var = "CXCR6", assay = "integrated" ,split.ncol = 5,
          cells.use = grep("Doublets|NK", NKT.Integrated_clean$anno.T_l2, invert = T),
          group.by = "group", split.by ="anno.T_l2", jitter.size = 0.1)+
  geom_boxplot(fill = "white", width = 0.3, outlier.shape = NA)+scale_fill_brewer(palette = "Set1")

th = round(quantile(NKT.Integrated_clean$TRM_UCell, 0.75), 1)
plot(density(NKT.Integrated_clean$TRM_UCell), main = "TRM signature")
abline(v =th, col="red", lwd=3)
FeaturePlot(NKT.Integrated_clean, features = "TRM_UCell" ,min.cutoff = th, pt.size = 0.8, order = T)+ggtitle("TRM rich T Cells")
Trm_rich = subset(NKT.Integrated_clean, 
                  cells = which(NKT.Integrated_clean$TRM_UCell > quantile(NKT.Integrated_clean$TRM_UCell, 0.75)))

th_cxcr6 = round(quantile(NKT.Integrated_clean$CXCR6_UCell, 0.75), 1)
dittoPlot(Trm_rich, var = "CXCR6_1", assay = "integrated" ,split.ncol = 5,slot = "data",
          cells.use = which(Trm_rich@assays$integrated@counts["CXCR6",] > 0 ),
          group.by = "group", jitter.size = 0.1,plots = c("vlnplot", "jitter","boxplot"))+
  ggtitle("TRM Rich cells")+
  stat_compare_means(label = "p.signif")


dittoPlot(NKT.Integrated_clean, var = "CXCR6", assay = "integrated" ,split.ncol = 5,slot = "data",
          cells.use = which(NKT.Integrated_clean@assays$integrated@counts["CXCR6",] > 0 ),
          group.by = "group", jitter.size = 0.1,plots = c("vlnplot", "jitter","boxplot"))+
  ggtitle("T cells")+
  stat_compare_means(label = "p.signif")
# geom_boxplot(fill = "white", width = 0.3, outlier.shape = NA)+scale_fill_brewer(palette = "Set1")


##

