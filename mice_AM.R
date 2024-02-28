library(Seurat)
library(org.Mm.eg.db)
library(dittoSeq)
library(ggplot2)
library(ggpubr)
library(reshape2)


mice = readRDS("C:/Projects/CXCR6/TRM/Mice Lung/LungPN1PN2_202203_Deb.rds")
mice$time = gsub("(T[0|1|3|7])-(HTO)([2|3|6|9])", "\\1", mice$HTO_classification)


DimPlot(mice, group.by = "integrated_snn_res.0.5", label = T)+
DimPlot(mice, group.by = "isCT", label = T)
FeaturePlot(mice,features = c("Siglecf" , "Pparg"), label = T, order = T)

Idents(mice) = "integrated_snn_res.0.5"
AM_sig = FindMarkers(mice, ident.1 = 13, assay = "integrated", test.use = "MAST")

table
DimPlot(mice, group.by = "celltype", label = T)

mice.alexis  = subset(mice, cells = grep("Alexis", mice$celltype))
table(mice$)

t0.mice.am = FindMarkers(subset(mice, cells = grep("T0", mice$time)), ident.1 = 13, assay = "integrated", test.use = "MAST")
t1.mice.am = FindMarkers(subset(mice, cells = grep("T1", mice$time)), ident.1 = 13,assay = "integrated", test.use = "MAST")
t3.mice.am = FindMarkers(subset(mice, cells = grep("T3", mice$time)), ident.1 = 13,assay = "integrated", test.use = "MAST" )
t7.mice.am = FindMarkers(subset(mice, cells = grep("T7", mice$time)), ident.1 = 13,assay = "integrated", test.use = "MAST" )
top_mice.am = Reduce(intersect, 
                      list(rownames(t0.mice.am)[t0.mice.am$avg_log2FC>0&t0.mice.am$p_val_adj<0.05],
                           rownames(t1.mice.am)[t1.mice.am$avg_log2FC>0&t1.mice.am$p_val_adj<0.05],
                           rownames(t3.mice.am)[t3.mice.am$avg_log2FC>0&t3.mice.am$p_val_adj<0.05],
                           rownames(t7.mice.am)[t7.mice.am$avg_log2FC>0&t7.mice.am$p_val_adj<0.05]))
top_mice.am
cat(sort(top_mice.am), sep = ", ")

full.am_sig = rownames(AM_sig)[AM_sig$avg_log2FC>0&AM_sig$p_val_adj<0.05]

length(intersect(top_mice.am,full.am_sig ))

saveRDS(top_mice.am , "TRM/Mice Lung/top_mice_am.rds")
saveRDS(full.am_sig , "TRM/Mice Lung/full_am_sig.rds")
