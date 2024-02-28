library(Seurat)
library(org.Mm.eg.db)
library(dittoSeq)
library(ggplot2)
library(ggpubr)
library(reshape2)


mice = readRDS("C:/Projects/CXCR6/TRM/Mice Lung/LungPN1PN2_202203_Deb.rds")
human_t = read.csv("C:/Projects/CXCR6/TRM/HUMAN TRM TUMOR/P26_MALE_Expression_table.tsv", sep = "\t", header = T, row.names = 1)
human_t.meta = as.data.frame(readxl::read_xlsx("C:/Projects/CXCR6/TRM/HUMAN TRM TUMOR/Sample_ID_RNAseq TRM.xlsx"))

rownames(human_t.meta) = human_t.meta$`Sample ID`
human_t.meta = human_t.meta[, -1]

head(rownames(mice))

gene.path = "C:/Projects/external/jem_20190249_tables1/"

trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS3.xlsx"),skip = 2,col_names = T,sheet = 2)
trm_table = trm_table[!is.na(trm_table$`Log2 fold change in human lung`),]
table(trm_table$`Log2 fold change in human lung`>0)
TRM.markers = trm_table$`Cheuk, et al. 2017 human skin TRM signature`[trm_table$`Log2 fold change in human lung`>0]



## tumor TRM versus. tumor non-TRM
trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS4.xlsx"),skip = 2,col_names = T,sheet = 1)
head(trm_table)
TRMtumor = trm_table[,c("Gene ID","log2 fold change", "P adj")]
TRMtumor = TRMtumor[TRMtumor$`log2 fold change`>0 &
                      TRMtumor$`P adj`<0.05, ]$`Gene ID`

## Lung TRM vs non TRM
trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T,sheet = 1)
TRMlung = trm_table[,c("Gene ID","log2 fold change", "P adj")]
TRMlung = TRMlung[TRMlung$`log2 fold change`>0 &
                    TRMlung$`P adj`<0.05, ]$`Gene ID`


TRMlung
TRMtumor


mice$time = gsub("(T[0|1|3|7])-(HTO)([2|3|6|9])", "\\1", mice$HTO_classification)

mice$annotation = mice$celltype
mice$annotation = gsub("AlexisT", "Alexis.T", mice$annotation)
mice$annotation[grep("undef|unknown|active|Cd8|naive", mice$annotation)] = "others"
mice$annotation = gsub("Alexis.", "", mice$annotation)


dittoPlot(mice, var = "Cxcr6", split.by = "celltype", group.by  = "time")
mus.genes = rownames(mice)
grep("Cxcl2", mus.genes)

head(mus.genes)


human2musGENE<-function(x){
  
  mart = read.csv("mart_export_human_mouse_pig.txt", sep = "\t", header = T)
  t(head(mart))
  idx = match(x, mart$Gene.name)
  
  if(sum(is.na(idx))){
    warning(paste(sum(is.na(idx)),"human gene names not found in mart: "), toString(x[is.na(idx)]))
  }
  
  idx = idx[complete.cases(idx)]
  mouse_genes = mart$Mouse.gene.name[idx]
  incomplete = which(mouse_genes==""|is.na(mouse_genes))
  
  warning(length(incomplete),"/" , length(x), " mouse gene homologs not found: ", 
          toString(mart$Gene.name[idx[incomplete]]))
  
  mouse_genes = mouse_genes[-incomplete]
  # mouse_genes_full = mouse_genes
  # mouse_genes_full[incomplete] = mart$Mouse.gene.stable.ID[idx[incomplete]]
  # mouse_genes = mouse_genes[which(mouse_genes!="")]
  return(mouse_genes)  
}

mus2humanGENE<-function(x){
  
  mart = read.csv("mart_export_mouse_human.txt", sep = "\t", header = T)
  t(head(mart))
  idx = match(x, mart$Gene.name)
  
  if(sum(is.na(idx))){
    warning(paste(sum(is.na(idx)),"mouse gene names not found in mart: "), toString(x[is.na(idx)]))
  }
  
  idx = idx[complete.cases(idx)]
  human_genes = mart$Human.gene.name[idx]
  incomplete = which(human_genes==""|is.na(human_genes))
  if(length(incomplete)){
    warning(length(incomplete),"/" , length(x), " human gene homologs not found: ", 
            toString(mart$Gene.name[idx[incomplete]]))
    human_genes = human_genes[-incomplete]
  }
  
  # mouse_genes_full = mouse_genes
  # mouse_genes_full[incomplete] = mart$Mouse.gene.stable.ID[idx[incomplete]]
  # mouse_genes = mouse_genes[which(mouse_genes!="")]
  return(human_genes)  
}

mice.trmSkin_homo = unique(human2musGENE(TRM.markers))
TRMlung.mice.homo = unique(human2musGENE(TRMlung))
TRMtumor.mice.homo = unique(human2musGENE(TRMtumor))

mice = AddModuleScore(mice, features = list("TRM" = mice.trmSkin_homo), name = "TRM_human_skin", seed = 10,assay = "RNA")
mice = AddModuleScore(mice, features = list("TRMlung" = TRMlung.mice.homo), name = "TRM_human_Lung", seed = 10,assay = "RNA")
mice = AddModuleScore(mice, features = list("TRMtumor" = TRMtumor.mice.homo), name = "TRM_human_Tumor", seed = 10,assay = "RNA")
mice = AddModuleScore(mice, features = list("TRM_union" = union(TRMlung.mice.homo,TRMtumor.mice.homo)), name = "TRM_human_total", seed = 10,assay = "RNA")


table(mice$reducedCellType, mice$celltype.SCINA)
DimPlot(mice, group.by = "celltype", label = T)+
  DimPlot(mice, group.by = "celltype.SCINA", label = T)

mice.alexis  = subset(mice, cells = grep("Alexis", mice$celltype))

dittoPlot(mice, var = "TRM_human_skin1", group.by = "reducedCellType")+geom_boxplot(fill="white", width = 0.5)+
  ggtitle("TRM Signature from HUMAN Dermis(JEM paper)")+NoLegend()+
  dittoPlot(mice.alexis, var = "TRM_human_skin1", group.by = "celltype")+geom_boxplot(fill="white", width = 0.5)+
  ggtitle("Reduced Annotation")+NoLegend()

dittoDotPlot(mice,vars = intersect(TRMlung.mice.homo, rownames(mice)), group.by = "reducedCellType" )
dittoDotPlot(mice,vars = intersect(TRMtumor.mice.homo, rownames(mice)),  group.by = "reducedCellType" )

table(mice$time)

table(mice.alexis$celltype)
Idents(mice.alexis) = "celltype"
library(ggrepel)
library(dplyr)
plotVolCano<-function(genes){
  genes$padj = genes$p_val_adj
  genes$log2FoldChange = genes$avg_log2FC
  genes$pvalue  = genes$p_val
  genes$Gene = rownames(genes)
  genes[which(genes[,"pvalue"] == 0), "pvalue"] <- .Machine$double.xmin
  
  genes$Significant <- ifelse(genes$padj < 0.05 & abs(genes$log2FoldChange) > 0, "FDR < 0.05", "Not Sig")
  head(genes)
  
  # idx = genes %>%  arrange(padj, -abs(log2FoldChange)) %>% 
  #   filter(Significant != "Not Sig") %>% 
  #   slice_head(n = 30) %>% 
  #   select(Gene)
  # genes$plot = "no"
  # genes[idx$Gene, ]$plot = "yes"
  
  ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("red", "grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel(
      data = subset(genes, genes$Significant != "Not Sig"),
      aes(label = Gene),
      size = 3.5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
}

table(mice.alexis$time)
# 
# T0   T1   T3   T7 
# 1016  419 1202 3492 
Idents(mice.alexis) = "celltype"

t0.mice.trm = FindMarkers(subset(mice.alexis, cells = grep("T0", mice.alexis$time)), ident.1 = "Alexis.Tmem")
t1.mice.trm = FindMarkers(subset(mice.alexis, cells = grep("T1", mice.alexis$time)), ident.1 = "Alexis.Tmem" )
t3.mice.trm = FindMarkers(subset(mice.alexis, cells = grep("T3", mice.alexis$time)), ident.1 = "Alexis.Tmem" )
t7.mice.trm = FindMarkers(subset(mice.alexis, cells = grep("T7", mice.alexis$time)), ident.1 = "Alexis.Tmem" )


p0 = plotVolCano(t0.mice.trm)+ggtitle("Alexis Tmem vs rest TCells in T0")
p1 = plotVolCano(t1.mice.trm)+ggtitle("Alexis Tmem vs rest TCells in T1")
p3 = plotVolCano(t3.mice.trm)+ggtitle("Alexis Tmem vs rest TCells in T3")
p7 = plotVolCano(t7.mice.trm)+ggtitle("Alexis Tmem vs rest TCells in T7")


gridExtra::grid.arrange(p0,p1,p3,p7, ncol = 4)

top_mice.trm = Reduce(intersect, 
                      list(rownames(t0.mice.trm)[t0.mice.trm$avg_log2FC>0&t0.mice.trm$p_val_adj<0.05],
                           rownames(t1.mice.trm)[t1.mice.trm$avg_log2FC>0&t1.mice.trm$p_val_adj<0.05],
                           rownames(t3.mice.trm)[t3.mice.trm$avg_log2FC>0&t3.mice.trm$p_val_adj<0.05],
                           rownames(t7.mice.trm)[t7.mice.trm$avg_log2FC>0&t7.mice.trm$p_val_adj<0.05]))
top_mice.trm
cat(sort(top_mice.trm), sep = ", ")


# 
# temp = Reduce(intersect, 
#        list(rownames(t0.mice.trm)[abs(t0.mice.trm$avg_log2FC)>0&t0.mice.trm$p_val_adj<0.05],
#             rownames(t1.mice.trm)[abs(t1.mice.trm$avg_log2FC)>0&t1.mice.trm$p_val_adj<0.05],
#             rownames(t3.mice.trm)[abs(t3.mice.trm$avg_log2FC)>0&t3.mice.trm$p_val_adj<0.05],
#             rownames(t7.mice.trm)[abs(t7.mice.trm$avg_log2FC)>0&t7.mice.trm$p_val_adj<0.05]))
# # setdiff(temp, top_mice.trm)
# # [1] "Rps15a" "Rps19"  "Rpl12"  "Rplp1"  "Rps18" 

intersect(top_mice.trm,TRMlung.mice.homo )
intersect(top_mice.trm,TRMtumor.mice.homo )

table(mice$time)

# trem_mice_signature 
saveRDS(list(mouse = sort(top_mice.trm), 
             human_homolog = mus2humanGENE(sort(top_mice.trm))), 
        file = "mice_TRM_sig.rds")
## Alternate DEG at each TRM-Tx vs rest of all cells

mice$pheno = paste(mice$celltype, mice$time, sep = "_")
Idents(mice) = "pheno"

Alexis.Tmem = subset(mice, cells = grep("Alexis.Tmem", mice$celltype))

table(Alexis.Tmem$pheno)

trm.mice.t0 = FindMarkers(Alexis.Tmem, ident.1 = "Alexis.Tmem_T0")
trm.mice.t1 = FindMarkers(Alexis.Tmem, ident.1 = "Alexis.Tmem_T1")
trm.mice.t3 = FindMarkers(Alexis.Tmem, ident.1 = "Alexis.Tmem_T3")
trm.mice.t7 = FindMarkers(Alexis.Tmem, ident.1 = "Alexis.Tmem_T7")



p0 = plotVolCano(trm.mice.t0)+ggtitle("T0 vs rest Tx within Tmem")
p1 = plotVolCano(trm.mice.t1)+ggtitle("T1 vs rest Tx within Tmem")
p3 = plotVolCano(trm.mice.t3)+ggtitle("T3 vs rest Tx within Tmem")
p7 = plotVolCano(trm.mice.t7)+ggtitle("T7 vs rest Tx within Tmem")


gridExtra::grid.arrange(p0,p1,p3,p7, ncol = 4)


## Alternate 3 Tmem vs rest

tmem.v.rest = FindMarkers(mice,ident.1 = "Alexis.Tmem", group.by = "celltype")
tmem.v.rest$Gene = rownames(tmem.v.rest)
mice.tmem.v.rest = rownames(tmem.v.rest)[tmem.v.rest$avg_log2FC>2 & tmem.v.rest$p_val_adj<0.05] 


mice = AddModuleScore(mice, features = list(mice.tmem.v.rest= mice.tmem.v.rest), name = "mice.tmem.v.rest", seed = 10,assay = "RNA")
FeaturePlot(mice, features = "mice.tmem.v.rest1", min.cutoff = 0.1)+ ggtitle("Tmem vs Rest Signature", subtitle = "from mice")
human2musGENE(mice.tmem.v.resthomo)

mice = AddModuleScore(mice, features = list("TMEM" = top_mice.trm), name = "TMEM", seed = 10,assay = "RNA")
FeaturePlot(mice, features = "TMEM1", min.cutoff = 0.1, order = T)+ggtitle("T.mem signature common Tx", subtitle = "from mice")

f1 = FeaturePlot(mice, features = "TMEM1", min.cutoff = 0.1, order = T)+ggtitle("T.mem signature common Tx", subtitle = "from mice")
f2 = FeaturePlot(mice, features = "TRM_human_Lung1", min.cutoff = 0.05, order = T)+
  ggtitle("Human Lung TRM",subtitle = "from literature")
f3 = FeaturePlot(mice, features = "TRM_human_Tumor1", min.cutoff = 0.05, order = T)+
  ggtitle("Human Tumor TRM",subtitle = "from literature")

d1 = DimPlot(mice, group.by = "celltype", label = T,label.box = T, label.size = 3)+NoLegend()

gridExtra::grid.arrange(d1, f1,f2,f3,ncol = 2)


# Alternate 4 , REdeine TRM interest to include cluster 17
p1 = DimPlot(mice, group.by = "seurat_clusters", label = T, label.box = T)+NoLegend()
p2 = DimPlot(mice, group.by = "isCT", label = T, label.box = T)+NoLegend()
p3 = DimPlot(mice, group.by = "celltype")


mice.alexis  = subset(mice, cells = grep("Alexis", mice$celltype))
DimPlot(mice.alexis)

mice.alexis= RunUMAP(mice.alexis, dim = 1:15)

DimPlot(mice.alexis)

clust17markers = FindMarkers(mice, ident.1 = c(5,17), group.by = "seurat_clusters")
clust17markers$Gene = rownames(clust17markers)
mice.clust17 = rownames(clust17markers)[clust17markers$avg_log2FC>0.8 & clust17markers$p_val_adj<0.05] 


mice = AddModuleScore(mice, features = list(clust17mar= mice.clust17), name = "clust17mar", seed = 10,assay = "RNA")
f17 = FeaturePlot(mice, features = "clust17mar1",min.cutoff = 0.1)


gridExtra::grid.arrange(p3,p2,p1,f17,  ncol = 2)



## End of mice





### crap

## Read Human Transcript Data 


# normalize by median
MedianNorm<-function(x){
  x/median(x, na.rm=T);
}


QuantileNormalize <- function(data){
  data = t(preprocessCore::normalize.quantiles(t(data), copy=FALSE))
  varCol <- apply(data, 2, var, na.rm=T);
  constCol <- (varCol == 0 | is.na(varCol));
  constNum <- sum(constCol, na.rm=T);
  if(constNum > 0){
    print(paste("After quantile normalization", constNum, "features with a constant value were found and deleted."));
    data <- data[,!constCol, drop=FALSE];
    colNames <- colnames(data);
    rowNames <- rownames(data);
  }
  return(data)
}

# Log Transform
# generalize log, tolerant to 0 and negative values
LogNorm<-function(x, min.val){
  # log10((x + sqrt(x^2 + min.val^2))/2)
  log2(x+1)
  
}

# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}



mart_export <- read.delim("C:/Projects/CXCR6/TRM/HUMAN TRM TUMOR/mart_export.txt")

mart_export = mart_export[, c("Gene.stable.ID", "Transcript.stable.ID","Transcript.stable.ID.version", "Gene.name")]
mart_export = mart_export[!duplicated(mart_export), ]
rownames(mart_export) = mart_export$Transcript.stable.ID.version
mart_export = mart_export[mart_export$Gene.name != "",]
mart_export = mart_export[complete.cases(mart_export),]

# mart_export = mart_export[!duplicated(mart_export$Gene.name), ]
a = data.frame(table(mart_export$Transcript.stable.ID.version))
a = a[order(a$Freq, decreasing = T),]
head(a)

t(mart_export[which(mart_export$Transcript.stable.ID.version == "ENST00000636378.2"), ])


matching_rows = intersect(rownames(human_t) , rownames(mart_export))
mart_export = mart_export[matching_rows, ]
human_match = human_t[matching_rows, ]
# rownames(human_match) = mart_export[matching_rows,]$Gene.name

which(mart_export$Gene.name=="A1BG")
human_match[mart_export$Transcript.stable.ID.version[which(mart_export$Gene.name=="A1BG")], ]

human_match$Gene = mart_export[rownames(human_match),]$Gene.name


human_raw = aggregate(.~Gene,data=human_match,FUN=sum)
rownames(human_raw) = human_raw$Gene
human_raw = human_raw[,-1]


total_expression = rowSums(human_raw)

dim(human_raw)
boxplot(human_raw)


data = t(human_raw)
dim(data)
data<-t(apply(data, 1, MedianNorm))
data<-QuantileNormalize(data)
dim(data)
boxplot(t(data))

min.val <- min(abs(data[data!=0]))/10
data<-log1p(data)
dim(data)
boxplot(t(data))

data<-apply(data, 1, AutoNorm)
dim(data)
boxplot(data)



human_aggregate<- t(data)
head(rownames(human_aggregate))
head(colnames(human_aggregate))

dim(human_aggregate)
human_aggregate = data.frame(human_aggregate)
intersect(TRM.markers, rownames(human_aggregate))

getSignature<-function(fpkm,gene.list){
  mat_interest = fpkm[intersect(gene.list, rownames(fpkm)) , ] 
  # mat_interest =log2(mat_interest[complete.cases(mat_interest),]+1)
  mat_interest = mat_interest[abs(rowSums(mat_interest))>0, ]
  signature = apply(scale(t(mat_interest), center = F),1, mean)
  return(signature)
}


plot.meta = data.frame(CXCR6 = getSignature(human_aggregate, "CXCR6"))
plot.meta = cbind(plot.meta, human_t.meta[rownames(plot.meta), ])
plot.meta$TRM_Sig = getSignature(human_aggregate, TRM.markers)

plot.meta$patient = as.factor(plot.meta$patient)
plot.meta$population = as.factor(plot.meta$population)


cor(plot.meta$TRM_Sig, plot.meta$CXCR6)
cor(getSignature(human_raw, TRM.markers), unlist(unlist(human_raw["CXCR6", ])))
cor(getSignature(human_aggregate, TRM.markers), unlist(unlist(human_aggregate["CXCR6", ])))


plot.data = plot.meta[,-which(colnames(plot.meta)=="CXCR6")]
head(plot.data)
plot.data = reshape::melt(plot.data)
head(plot.data)




sp <- ggscatter(plot.meta, x = "CXCR6", y = "TRM_Sig", color = "population",
                add = "reg.line",  
                # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)

ggplot(plot.data, aes(x = population, y = value, fill = population))+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.05)+ylab("TRM")+
  theme_classic()+theme(legend.position = "none")+ggtitle("Human Tumour Cells")+
  stat_compare_means(label = "p.format", method = "wilcox.test")+
  sp+ stat_cor(aes(color = population))


dittoPlot(mice, var = "Cxcr6", split.by = "celltype", group.by  = "time", boxplot.width = 0.3,
          plots = c("vlnplot", "boxplot", "jitter"),jitter.size = 0.2)+
  ggtitle("Mice Lung")+
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref = "T7")

# Define Gene signature from Mice Alexis Tmem CXCR6+

table(mice$reducedCellType, mice$celltype)

mice.tmem = subset(mice, cells = grep("Alexis.Tmem", mice$celltype))
mice.tmem@active.assay = "RNA"
mice.tmem = NormalizeData(mice.tmem)
mice.tmem = FindVariableFeatures(mice.tmem)
mice.tmem = ScaleData(mice.tmem)
mice.tmem = RunPCA(mice.tmem)
ElbowPlot(mice.tmem)
mice.tmem = RunUMAP(mice.tmem, dim = 3:10)
DimPlot(mice.tmem, group.by = "time")

table(mice.tmem$time)

t0.mice.tmem = subset(mice.tmem, cells = grep("T0", mice.tmem$time))
t1.mice.tmem = subset(mice.tmem, cells = grep("T1", mice.tmem$time))
t3.mice.tmem = subset(mice.tmem, cells = grep("T3", mice.tmem$time))
t7.mice.tmem = subset(mice.tmem, cells = grep("T7", mice.tmem$time))

library(scWGCNA)

t0.mice.tmem = NormalizeData(t0.mice.tmem)
t0.mice.tmem = FindVariableFeatures(t0.mice.tmem)
t0.mice.tmem = ScaleData(t0.mice.tmem)
t0.mice.tmem = RunPCA(t0.mice.tmem,npcs = 20)
T0.pells = calculate.pseudocells(s.cells = t0.mice.tmem, # Single cells in Seurat object
                                 seeds=0.5, # Fraction of cells to use as seeds to aggregate pseudocells
                                 nn = 5, # Number of neighbors to aggregate
                                 reduction = "pca", # Reduction to use
                                 dims = 1:10) # The dimensions to use
t0.mice.tmem.scWGCNA = run.scWGCNA(p.cells = T0.pells, # Pseudocells (recommended), or Seurat single cells
                                   s.cells = t0.mice.tmem, # single cells in Seurat format  
                                   is.pseudocell = T # We are using single cells twice this time
                                   # features = VariableFeatures(t0.mice.tmem)
)
