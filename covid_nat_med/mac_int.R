library(Seurat)
library(UCell)
library(ggplot2)
library(ggpubr)
library(dittoSeq)
library(ggrepel)
library(tidyverse)

# DefaultAssay(mac) = "RNA"
# table(mac$sample, mac$group)
# mac.list <- SplitObject(mac, split.by = "sample")
# mac.list <- lapply(X = mac.list, SCTransform, variable.features.n = NULL, variable.features.rv.th = 1.1)
# 
# 
# features <- SelectIntegrationFeatures(object.list = mac.list, nfeatures = 3000)
# 
# # res.features = union(features, lungTRM_noTRM.genes)
# mac.list <- PrepSCTIntegration(object.list = mac.list,  anchor.features = features)
# mac.anchors <- FindIntegrationAnchors(object.list = mac.list, normalization.method = "SCT",
#                                          anchor.features = features)
# mac.combined.sct <- IntegrateData(anchorset = mac.anchors, normalization.method = "SCT")
# 

# 
# saveRDS(mac.combined.sct, "mac_int_heal_SCT.rds")

# mac.combined.sct = readRDS("covid_nat_med/mac_int_heal_SCT.rds")
mac.combined.sct = readRDS("covid_nat_med/mac.rds")
DimPlot(mac.combined.sct, group.by = "condition")

DefaultAssay(mac.combined.sct) = "integrated" 

mac.combined.sct$condition = mac.combined.sct$group
mac.combined.sct$condition = gsub("HC", "Healthy", mac.combined.sct$condition)
mac.combined.sct$condition = gsub("^O", "Moderate", mac.combined.sct$condition)
mac.combined.sct$condition = gsub("S/C", "Severe", mac.combined.sct$condition)

table(mac.combined.sct$sample, mac.combined.sct$condition)
table(mac.combined.sct$cluster)


# mac.Integrated <- ScaleData(mac.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
mac.combined.sct= RunPCA(mac.combined.sct,reduction.name = "pca_re")
ElbowPlot(mac.combined.sct, reduction = "pca_re")
mac.combined.sct = RunUMAP(mac.combined.sct, dims = 1:10, reduction = "pca",reduction.name = "umap_re")

DimPlot(mac.combined.sct, group.by = "condition", reduction = "umap_re")


CXCL16_th = 1.3
mac.combined.sct$CXCL16pos = "neg"
mac.combined.sct$CXCL16pos[mac.combined.sct[["integrated"]]@data["CXCL16", ] > CXCL16_th] = "pos"
dittoBarPlot(mac.combined.sct, var = "CXCL16pos", group.by = "condition" )
df = data.frame(CXCL16pos = mac.combined.sct$CXCL16pos, condition = mac.combined.sct$condition)
head(df)
plotPie(df, "CXCL16")

pdf("covid_nat_med/mac_CXCL16.pdf", width = 11, height = 3.5)
# RidgePlot(mac.combined.sct,"CXCL16", group.by = "condition")+
#   geom_vline(xintercept = CXCL16_th, col = "cyan", size = 1)+NoLegend()+
  FeaturePlot(mac.combined.sct, "CXCL16", min.cutoff = CXCL16_th, split.by = "condition",
              pt.size = 1, order = T, reduction = "umap")
dev.off()


CCL7_th = 1.5
mac.combined.sct$CCL7pos = "neg"
mac.combined.sct$CCL7pos[mac.combined.sct[["integrated"]]@data["CCL7", ] > CCL7_th] = "pos"
dittoBarPlot(mac.combined.sct, var = "CCL7pos", group.by = "condition" )
df_CCL7 = data.frame(CCL7pos = mac.combined.sct$CCL7pos, condition = mac.combined.sct$condition)
head(df_CCL7)
plotPie(df_CCL7, "CCL7")

pdf("covid_nat_med/mac_CCL7.pdf", width = 11, height = 3.5)
# RidgePlot(mac.combined.sct,"CCL7", group.by = "condition")+
#   geom_vline(xintercept = CCL7_th, col = "cyan", size = 1)+NoLegend()+
  FeaturePlot(mac.combined.sct, "CCL7", min.cutoff = CCL7_th, split.by = "condition",
              pt.size = 1, order = T, reduction = "umap")
dev.off()


CCL2_th = 3.5
mac.combined.sct$CCL2pos = "neg"
mac.combined.sct$CCL2pos[mac.combined.sct[["integrated"]]@data["CCL2", ] > CCL2_th] = "pos"
dittoBarPlot(mac.combined.sct, var = "CCL2pos", group.by = "condition" )
df_CCL2 = data.frame(CCL2pos = mac.combined.sct$CCL2pos, condition = mac.combined.sct$condition)
head(df_CCL2)
plotPie(df_CCL2, "CCL2")

pdf("covid_nat_med/mac_CCL2.pdf", width = 11, height = 3.5)
# RidgePlot(mac.combined.sct,"CCL2", group.by = "condition")+
#   geom_vline(xintercept = CCL2_th, col = "cyan", size = 1)+NoLegend()+
  FeaturePlot(mac.combined.sct, "CCL2", min.cutoff = CCL2_th, split.by = "condition",
              pt.size = 1, order = T, reduction = "umap")
dev.off()



my_comparisons = list(c("Healthy","Moderate"),c("Healthy","Severe"),c("Moderate","Severe"))

d1 = dittoPlot(mac.combined.sct, var  = c("CXCL16"), 
               cells.use = which(mac.combined.sct@assays$integrated@data["CXCL16", ]>0),
               assay = "integrated",
               plots = c("jitter", "vlnplot", "boxplot"),max = 5.5,
               jitter.size = 0.1, jitter.color = "#00000044",
               group.by = "condition")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("Macrophages")+scale_fill_brewer()+NoLegend()

d2 = dittoPlot(mac.combined.sct, var  = c("CCL7"), 
               cells.use = which(mac.combined.sct@assays$integrated@data["CCL7", ] > 0),
               assay = "integrated", 
               plots = c("jitter", "vlnplot", "boxplot"),max = 8,
               jitter.size = 0.1, jitter.color = "#00000044",
               group.by = "condition")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("Macrophages")+scale_fill_brewer()+NoLegend()

d3 = dittoPlot(mac.combined.sct, var  = c("CCL2"), 
               cells.use = which(mac.combined.sct@assays$integrated@data["CCL2", ] >0),
               assay = "integrated", 
               plots = c("jitter", "vlnplot", "boxplot"),max = 11,
               jitter.size = 0.1, jitter.color = "#00000044",
               group.by = "condition")+
  stat_compare_means(label = "p.signif", tip.length = 0.03,
                     comparisons = my_comparisons)+
  ggtitle("Macrophages")+scale_fill_brewer()+NoLegend()

pdf("covid_nat_med/mac_genes_vln.pdf", width = 9, height = 3.5)
gridExtra::grid.arrange(d1,d2,d3, ncol = 3)
dev.off()



plotPie<-function(df, gene_name = "CXCL16"){
  
  colnames(df) = c(gene_name, "condition")
  # df[, gene_name] = factor(df[, gene_name], levels = c("pos", "neg")) 
  
  df1 = df %>% group_by(condition, UQ(sym(gene_name))) %>% # Variable to be transformed
    count() %>% 
    ungroup( UQ(sym(gene_name))) %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc)) %>%
    arrange(condition)
  
  
  df2 <- df1 %>% 
    mutate(csum = rev(cumsum(rev(perc))), 
           pos = perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), perc/2, pos)) %>% arrange(condition, desc(UQ(sym(gene_name))))
  df2 = as.data.frame(df2)
  pdf(paste0("covid_nat_med/mac_",gene_name, "_prop.pdf"), width = 10, height = 3)
  par(mfrow = c(1, 3))
  pie(labels = df2[1:2, 2],x = df2[1:2, 4], main = paste( gene_name,df2[1, 1], sep ="\n"), 
      col = RColorBrewer::brewer.pal(3, "Set1") )
  pie(labels = df2[3:4, 2],x = df2[3:4,4], main = paste( gene_name,df2[ 3, 1], sep ="\n") ,
      col = RColorBrewer::brewer.pal(3, "Set1"))
  pie(labels = df2[5:6, 2],x = df2[5:6,4], main = paste( gene_name,df2[5,1], sep ="\n") ,
      col = RColorBrewer::brewer.pal(3, "Set1"))
  dev.off()
  
}


# gg1 = base2grob(~pie(labels = df2[1:2, 2],x = df2[1:2, 4], main = paste( gene_name,df2[1, 1], sep ="\n"), 
#                      col = RColorBrewer::brewer.pal(3, "Set1") ))
# gg2 = base2grob(~pie(labels = df2[3:4, 2],x = df2[3:4,4], main = paste( gene_name,df2[ 3, 1], sep ="\n") ,
#                      col = RColorBrewer::brewer.pal(3, "Set1")))
# gg3  =  base2grob(~pie(labels = df2[5:6, 2],x = df2[5:6,4], main = paste( gene_name,df2[5,1], sep ="\n") ,
#                        col = RColorBrewer::brewer.pal(3, "Set1")))


# ggplot(df1[1:2,], aes(x = "" , y = perc, fill = fct_inorder(CXCL16pos))) +
#   geom_col(width = 1, color = 1) +
#   coord_polar(theta = "y") +
#   scale_fill_brewer(palette = "Set1") +
#   geom_text_repel(data = df2,
#                    aes(y = pos, label = CXCL16pos),
#                    size = 4.5, nudge_x = 0.7, show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Group")) +
#   theme_void()


