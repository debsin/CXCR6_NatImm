
setwd("C:/Projects/CXCR6/covid_szabo/sims-farber/")
library(Seurat)
library(dittoSeq)
library(ggpubr)


gene.path = "C:/Projects/external/jem_20190249_tables1/"

lungTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T)
lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[lungTRM_noTRM$`log2 fold change`>2 & lungTRM_noTRM$`P adj`<0.01]
# lungTRM_noTRM.genes = setdiff(lungTRM_noTRM.genes, "CXCR6")

# pbmc = readRDS("pbmc.rds")
aww = readRDS("aww.rds")
table(aww$Sample, aww$Timepoint)
tcell = subset(aww, cells = grep("Tcell", aww$Annotation))
rownames(tcell[["RNA"]]@counts) = gsub("GRCh38---","", rownames(tcell[["RNA"]]@counts))
rownames(tcell[["RNA"]]@data) = gsub("GRCh38---","", rownames(tcell[["RNA"]]@data))
rownames(tcell[["RNA"]]@scale.data) = gsub("GRCh38---","", rownames(tcell[["RNA"]]@scale.data))


tcell = NormalizeData(tcell, normalization.method = "LogNormalize", scale.factor = 1e4)
tcell = FindVariableFeatures(tcell, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
tcell@assays$RNA@var.features = union(tcell@assays$RNA@var.features, lungTRM_noTRM.genes)
tcell <- ScaleData(tcell, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
tcell = RunPCA(tcell, verbose = FALSE)
ElbowPlot(object = tcell,ndims = 50)
tcell <- RunUMAP(tcell, reduction = "pca", dims = 1:50, seed.use = 100, n.neighbors = 20)
DimPlot(tcell, group.by = "Patient", label = T, label.box = T)

DefaultAssay(tcell) <- "RNA"
SubseTs.list <- SplitObject(tcell, split.by = "Patient")
for (i in 1:length(SubseTs.list)) {
  SubseTs.list[[i]] <- NormalizeData(SubseTs.list[[i]], verbose = FALSE)
  SubseTs.list[[i]] <- FindVariableFeatures(SubseTs.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  SubseTs.list[[i]]@assays$RNA@var.features = union(SubseTs.list[[i]]@assays$RNA@var.features, lungTRM_noTRM.genes)
}

patient_name = c('COV022','COV026','COV027','COV028')
reference.list <- SubseTs.list[patient_name]
tcell.temp <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50,k.filter = 140)
tcell.int <- IntegrateData(anchorset = tcell.temp, dims = 1:50)


###first generate data and scaledata in RNA assay
DefaultAssay(tcell.int) <- "RNA"
tcell.int[['percent.mito']] <- PercentageFeatureSet(tcell.int, pattern = "^MT-")
tcell.int <- NormalizeData(object = tcell.int, normalization.method = "LogNormalize", scale.factor = 1e4)
tcell.int <- FindVariableFeatures(object = tcell.int, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
tcell.int@assays$integrated@var.features = union(tcell.int@assays$integrated@var.features, lungTRM_noTRM.genes)
tcell.int <- ScaleData(tcell.int, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(tcell.int) <- "integrated"
VlnPlot(object = tcell.int, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(object = tcell.int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Run the standard workflow for visualization and clustering
tcell.int <- FindVariableFeatures(object = tcell.int, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
tcell.int <- ScaleData(tcell.int, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
tcell.int <- RunPCA(tcell.int, verbose = FALSE)
#visulaization pca result
tcell.int <- ProjectDim(object = tcell.int)
ElbowPlot(object = tcell.int,ndims = 50)

tcell.int <- RunUMAP(tcell.int, reduction = "pca", dims = 1:50, seed.use = 100, n.neighbors = 20)
DimPlot(tcell.int,  group.by = "Patient", label = T, label.box = T)

tcell.int@assays$integrated@counts = tcell.int@assays$integrated@data

saveRDS(tcell.int, "tcell_int.rds")

## T Cell Figures

tcell.int = AddModuleScore(tcell.int, assay = "integrated",ctrl = 50,
                           features = list("TRM" = lungTRM_noTRM.genes), name = "TRM_")

tcell.int = AddModuleScore(tcell.int, assay = "integrated",ctrl = 50,
                           features = list("CXCR6" = c("CXCR6")), name = "CXCR6_")

DimPlot(tcell.int, group.by = "Timepoint", label = T, label.box = T)

plot(density(tcell.int$TRM_1,adjust = 0.1), main = "TRM Signature")
hist(tcell.int$TRM_1,breaks = 100, main = "TRM Signature", xlab = "signature score")
abline(v = quantile(tcell.int$TRM_1, 0.77), col = "red", lwd = 2)
text(0.1, 140, "quantile @ 0.8")

# plot(density(tcell.int$CXCR6_1), main = "CXCR6 Signature")
hist(tcell.int$CXCR6_1,breaks = 100, main = "CXCR6 Signature",  xlab = "signature score")
abline(v = quantile(tcell.int$CXCR6_1, 0.8), col = "red", lwd = 2)
text(0.5, 400, "quantile @ 0.8")

trm_th = 0.8
cx_th = 0.75
FeaturePlot(tcell.int, features = "TRM_1", 
            min.cutoff = quantile(tcell.int$TRM_1, trm_th), order = T, ncol = 3)+ggtitle("TRM Signature")+
  FeaturePlot(tcell.int, features = "CXCR6_1", 
              min.cutoff = quantile(tcell.int$CXCR6_1, cx_th), order = T)+ggtitle("CXCR6 Signature")+
  FeatureScatter(object = tcell.int, feature1 = "TRM_1", feature2 = "CXCR6_1", group.by = "Timepoint",
                 # cells = intersect(which(tcell.int$TRM_1 > quantile(tcell.int$TRM_1, trm_th)) ,
                 #         which(tcell.int$CXCR6_1 > quantile(tcell.int$CXCR6_1, cx_th))),
                 shuffle = T)+
  xlab("TRM Signature Score")+ylab("CXCR6")



FeatureScatter(object = tcell.int, feature1 = "TRM_1", feature2 = "CXCR6_1", group.by = "Timepoint",
               cells = intersect(which(tcell.int$TRM_1 > quantile(tcell.int$TRM_1, 0.8)),
                                 which(tcell.int$CXCR6_1 > quantile(tcell.int$CXCR6_1, 0.75))),
               shuffle = T)

Trm_rich = subset(tcell.int, 
                  cells = which(tcell.int$TRM_1 > quantile(tcell.int$TRM_1, 0.5)))
dittoPlot(Trm_rich, var = "CXCR6_1", assay = "integrated" ,split.ncol = 5,slot = "data",
          cells.use = which(Trm_rich$CXCR6_1 > quantile(Trm_rich$CXCR6_1 ,0.75) ),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "Timepoint", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("TRM Rich cells")+
  stat_compare_means(label = "p.signif", ref = "O")

table(tcell.int$Timepoint)
tcell.int$Time = "Early"
tcell.int$Time[grep("4|5", tcell.int$Timepoint)] = "Mid"
tcell.int$Time[grep("6|7", tcell.int$Timepoint)] = "Late"

tcell.int$Time = factor(tcell.int$Time, levels = c("Early", "Mid", "Late"))
dittoPlot(tcell.int, var = "CXCR6", assay = "integrated" ,split.ncol = 5,slot = "data",
          # cells.use = which(tcell.int$CXCR6_1 > quantile(tcell.int$CXCR6_1 ,0.75) ),
          jitter.width = 0.3,jitter.color = "#00000077",boxplot.width = 0.8,
          group.by = "Time", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("T cells")+theme(legend.position = "None")+
  stat_compare_means(label = "p.signif",  ref = "Mid")

dittoPlot(tcell.int, var = "CXCR6_1", assay ="integrated" ,split.ncol = 5,slot = "data",
          cells.use = intersect(which(tcell.int$CXCR6_1 > quantile(tcell.int$CXCR6_1 ,0.75)),
                                which(tcell.int$TRM_1 > quantile(tcell.int$TRM_1 ,0.77) )),
          jitter.width = 0.3,jitter.color = "#00000088",boxplot.width = 0.8,
          group.by = "Timepoint", jitter.size = 0.1,plots = c("boxplot", "jitter"))+
  ggtitle("CXCR6 TRM cells")+theme(legend.position = "None")+
  stat_compare_means(label = "p.signif",  ref = "IntubationDay1")

table(tcell.int$anno.T_l2)

