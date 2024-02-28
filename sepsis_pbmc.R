library(Seurat)
library(org.Mm.eg.db)
library(dittoSeq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(reshape2)


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


sepsis = readRDS("C:/Projects/sepsis/pbmc.rds")

# sepsis[['SCT']] <- NULL

head(sepsis@meta.data)

table(sepsis$phenotype)


am = list()
am$mice = readRDS("TRM/Mice Lung/top_mice_am.rds")
am$human_homolog = mus2humanGENE(sort(am$mice))


sepsis= AddModuleScore(sepsis, features = list("AM" = am$human_homolog), name = "AM")
colnames(sepsis@meta.data) = gsub("AM1", "AM", colnames(sepsis@meta.data))
colnames(sepsis@meta.data)

table(sepsis$predicted.celltype.l1)

plot_df = sepsis@meta.data[intersect(grep("other", sepsis$predicted.celltype.l1, invert = T),
                                     grep("Control|Bac", sepsis$phenotype)), ]

plot_df$condition  = factor(plot_df$phenotype , levels = c("Control", "Bac-SEP"))

gsep =ggplot(plot_df[plot_df$predicted.celltype.l1=="Mono",], aes(x = condition, y = AM, fill = condition))+
  # facet_wrap(~predicted.celltype.l1, ncol = 6)+
  geom_violin(show.legend = F)+
  geom_jitter(width = 0.2, size = 0.1, show.legend = F,)+
  geom_boxplot(outlier.shape = NA, show.legend = F, width =0.5, fill = "#FFFFFFAA")+
  scale_fill_hue(direction = -1)+
  # geom_jitter(width = 0.1, show.legend = F)+
  stat_compare_means(ref.group = "Control",   method = "t.test", label = "p.format")+
  ylab("AM signature from mice")+
  xlab("Monocytes")+
  theme_classic()

pdf("sepsis_AM.pdf", width = 6, height = 5)
gridExtra::grid.arrange(gsep, ncol = 1, top=textGrob("Sepsis PBMC Dataset, Reyes et al. "))
dev.off()
