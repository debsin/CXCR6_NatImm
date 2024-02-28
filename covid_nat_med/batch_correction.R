
library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)

gene.path = "C:/Projects/external/jem_20190249_tables1/"

lungTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T)
lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[lungTRM_noTRM$`log2 fold change`>2 & lungTRM_noTRM$`P adj`<0.01]
# lungTRM_noTRM.genes = setdiff(lungTRM_noTRM.genes, "CXCR6")

# 
# cov.T =  readRDS("covid_nat_med/T_all.rds")
# 
# count = GetAssayData(cov.T, "RNA",slot =  "count")
# meta = droplevels(cov.T@meta.data)

T.cells = CreateSeuratObject(count, project = "covid", assay = "RNA", meta.data = meta)

DefaultAssay(T.cells) = "RNA"
table(T.cells$sample, T.cells$anno.T_l2)
T.cells.list <- SplitObject(T.cells, split.by = "sample")
T.cells.list <- lapply(X = T.cells.list, SCTransform, variable.features.n = NULL, variable.features.rv.th = 1.1)
# missing = intersect(setdiff(lungTRM_noTRM.genes, SelectIntegrationFeatures(object.list = T.cells.list, nfeatures = 5000)), rownames(count))
present2k = intersect(lungTRM_noTRM.genes, SelectIntegrationFeatures(object.list = T.cells.list, nfeatures = 3000))
# present6k = intersect(lungTRM_noTRM.genes, SelectIntegrationFeatures(object.list = T.cells.list, nfeatures = 4000))
setdiff(present6k, present2k)

VlnPlot(T.cells, features = setdiff(present6k, present2k), assay = "RNA")

features <- SelectIntegrationFeatures(object.list = T.cells.list, nfeatures = 3000)

# res.features = union(features, lungTRM_noTRM.genes)
T.cells.list <- PrepSCTIntegration(object.list = T.cells.list,  anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = T.cells.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")


immune.combined.sct =  readRDS("covid_nat_med/T_all.rds")
DefaultAssay(immune.combined.sct) = "integrated"

immune.combined.sct$condition = immune.combined.sct$group
immune.combined.sct$condition = gsub("HC", "Healthy", immune.combined.sct$condition)
immune.combined.sct$condition = gsub("^O", "Moderate", immune.combined.sct$condition)
immune.combined.sct$condition = gsub("S/C", "Severe", immune.combined.sct$condition)



immune.combined.sct= RunPCA(immune.combined.sct)
immune.combined.sct = RunUMAP(immune.combined.sct, dims = 1:50)
DimPlot(immune.combined.sct, group.by = "anno.T_l2", label.box = T, label = T)+NoLegend()

immune.combined.sct = AddModuleScore(immune.combined.sct, features = list(lungTRM_noTRM.genes), name = "TRM_", ctrl = 20)
immune.combined.sct = AddModuleScore(immune.combined.sct, features = "CXCR6",assay = "integrated",  name = "CXCR6_", ctrl = 50)

immune.combined.sct = AddModuleScore_UCell(immune.combined.sct, features = list("TRM" = lungTRM_noTRM.genes,
                                                                                "CXCR6" = c("CXCR6")))


VlnPlot(immune.combined.sct, features = c("TRM_UCell", "CXCR6_UCell"), group.by = "condition")


df = data.frame(CXCR6 = immune.combined.sct[["integrated"]]@data["CXCR6", ], 
                condition = immune.combined.sct$condition)


f1 = FeaturePlot(immune.combined.sct, features = "CXCR6", slot = "data", 
                 min.cutoff = 0,order = T, pt.size = 1)
f2 = FeaturePlot(immune.combined.sct, features = "TRM_1", 
                 min.cutoff = 0,order = T, pt.size = 1)



f3 = FeatureScatter(immune.combined.sct, feature1 = "CXCR6", feature2 = "TRM_1", group.by = "condition", jitter = T)

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

r1 = RidgePlot(immune.combined.sct, features = "CXCR6", group.by = "condition")+
  geom_vline(xintercept = 1, col = "cyan", size = 1)+NoLegend()
r2 = RidgePlot(immune.combined.sct, features = "TRM_1", group.by = "condition")+
  geom_vline(xintercept = 0.05, col = "cyan", size = 1)+ggtitle("TRM Signature")+
  NoLegend()

pdf("covid_nat_med/TCells_CXCR6_TRMcutoff.pdf", width = 15, height = 3.5)
gridExtra::grid.arrange(r1,r2, ncol = 2)
dev.off()

dittoBarPlot(immune.combined.sct, var = "TRMhi", group.by = "condition")

VlnPlot(immune.combined.sct, features = "CXCR6", group.by = "CXCR6hi", pt.size = 0.1)
FeatureScatter(immune.combined.sct,feature1 = "CXCR6_UCell", feature2 = "CXCR6", group.by = "CXCR6hi")
FeatureScatter(immune.combined.sct,feature1 = "TRM_UCell", feature2 = "TRM_1", group.by = "TRMhi")


colnames(immune.combined.sct@meta.data)[grep("TRM_1",colnames(immune.combined.sct@meta.data) )] = "TRM Signature"

pdf("covid_nat_med/TCells_umap_TRM.pdf", width = 11, height = 3.5)
# RidgePlot(immune.combined.sct,features = "TRM Signature", group.by = "condition")+
#   geom_vline(xintercept = TRMhi_th, col = "cyan", size = 1)+NoLegend()+
FeaturePlot(immune.combined.sct, features = "TRM Signature", slot = "data", split.by = "condition", 
            min.cutoff = TRMhi_th,order = T, pt.size = 1, keep.scale = "feature")

dev.off()

pdf("covid_nat_med/TCells_umap_CXCR6.pdf", width = 11, height = 3.5)
# RidgePlot(immune.combined.sct,features = "CXCR6", group.by = "condition")+
#   geom_vline(xintercept = CXCR6hi_th, col = "cyan", size = 1)+NoLegend()+
  FeaturePlot(immune.combined.sct, features = "CXCR6", slot = "data", split.by = "condition", 
              min.cutoff = CXCR6hi_th,order = T, pt.size = 1, keep.scale = "feature")

dev.off()

FeaturePlot(immune.combined.sct, features = "CXCR6", slot = "data", split.by = "condition", 
            min.cutoff = 1,order = T, pt.size = 1, keep.scale = "feature")

immune.combined.sct$CXCR6_TRM = "others"
immune.combined.sct$CXCR6_TRM[intersect(CXCR6hi_cells, TRMhi_cells)] = "CXCR6hi TRM"
immune.combined.sct$CXCR6_TRM = factor(immune.combined.sct$CXCR6_TRM, levels = c( "CXCR6hi TRM", "others"))

table(immune.combined.sct$sample, immune.combined.sct$condition)

pdf("covid_nat_med/TCells_umap_subtype.pdf", width = 5, height = 4)
DimPlot(immune.combined.sct, group.by = "anno.T_l2", order = T)
dev.off()

pdf("covid_nat_med/TCells_umap_binary.pdf", width = 12, height = 3)
DimPlot(immune.combined.sct, group.by = "TRMhi", order = F)+
  DimPlot(immune.combined.sct, group.by = "CXCR6hi", order = F)+
  DimPlot(immune.combined.sct,group.by = "CXCR6_TRM", order = F)
dev.off()


Idents(immune.combined.sct)  = "condition"
f1 = FeatureScatter(immune.combined.sct, feature1 = "CXCR6", feature2 = "TRM_1", group.by = "condition", jitter = T, )
f2 = FeatureScatter(immune.combined.sct, feature1 = "CXCR6", feature2 = "TRM_1",shuffle = T,   
                    cells = intersect(CXCR6hi_cells, grep("Healthy", immune.combined.sct$condition)))
f3 = FeatureScatter(immune.combined.sct, feature1 = "CXCR6", feature2 = "TRM_1",shuffle = T,   
                    cells = intersect(CXCR6hi_cells, grep("Moderate", immune.combined.sct$condition)))
f4 = FeatureScatter(immune.combined.sct, feature1 = "CXCR6", feature2 = "TRM_1",shuffle = T,   
                    cells = intersect(CXCR6hi_cells, grep("Severe", immune.combined.sct$condition)))
pdf("covid_nat_med/TCells_correlation.pdf", width = 16, height = 3)
gridExtra::grid.arrange(f1,f2,f3,f4, ncol =4)
dev.off()

my_comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe"),c("Moderate","Severe"))

pdf("covid_nat_med/TCells_CXCR6.pdf", width = 6, height = 3.5)
dittoPlot(immune.combined.sct, var = "CXCR6", assay = "integrated" ,slot = "data", 
          cells.use = CXCR6hi_cells,
          plots = c("jitter", "vlnplot", "boxplot"),max = 4.5,
          jitter.size = 0.1, jitter.color = "#00000044",
          group.by = "condition")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("T cells")+scale_fill_brewer()+NoLegend()+

dittoPlot(immune.combined.sct, var = "CXCR6", assay = "integrated" ,slot = "data", 
          cells.use = intersect(TRMhi_cells , CXCR6hi_cells),
          plots = c("jitter", "vlnplot", "boxplot"),max =4.5,
          jitter.size = 0.1, jitter.color = "#00000044",
          group.by = "condition")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("TRM cells")+scale_fill_brewer()+NoLegend()
dev.off()



