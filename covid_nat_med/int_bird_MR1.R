##integrate data with seurat v3
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
cov = readRDS(file = "nCoV.rds")
cov_rna = CreateSeuratObject(counts = cov@assays$RNA@counts, min.cells = 3, min.features = 200, project = "liao", meta.data = cov@meta.data)
##############Macrophage###############
DefaultAssay(cov_rna) <- "RNA"
nCoV.list <- SplitObject(cov_rna, split.by = "sample")
for (i in 1:length(nCoV.list)) {
  nCoV.list[[i]] <- NormalizeData(nCoV.list[[i]], verbose = FALSE)
  nCoV.list[[i]] <- FindVariableFeatures(nCoV.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  nCoV.list[[i]]@assays$RNA@var.features <-  union(nCoV.list[[i]]@assays$RNA@var.features, c("MR1", "HLA-DRA", "HLA-A", "CD3E"))
  
}
# nCoV.list = list()
# for(sample_s in unique(cov_rna$sample)){
#   print(sample_s)
#   sample_i = samples %>% dplyr::filter(.,sample == sample_s)
#   datadir = paste("/home/data/results/workspace/COVID_matrix/",sample_s,"/outs/filtered_feature_bc_matrix/",sep="")
#   sample.tmp = Read10X(data.dir = datadir)
#   sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = sample_s)
#   sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
#   sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
#   sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
#   sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
#   sample_i$percent.mito = as.numeric(sample_i$percent.mito)
#   sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > sample_i$nFeature_RNA_low & nFeature_RNA < sample_i$nFeature_RNA_high 
#                               & nCount_RNA > sample_i$nCount_RNA & percent.mito < sample_i$percent.mito)
#   sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
#   sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
#   nCoV.list[sample_s] = sample.tmp.seurat
# }
nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))

###first generate data and scale data in RNA assay
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(nCoV.integrated) <- "integrated"
dpi = 300
png(file="qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = nCoV.integrated,ndims = 100)
dev.off()

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2) 

###tsne and umap
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- RunUMAP(nCoV.integrated, reduction = "pca", dims = 1:50)
png(file="tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

# DefaultAssay(nCoV.integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones
# nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
# write.table(nCoV.integrated@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

# dpi = 300
# png(file="feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
# VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
# dev.off()
saveRDS(nCoV.integrated, file = "nCoV_MR1.rds")
