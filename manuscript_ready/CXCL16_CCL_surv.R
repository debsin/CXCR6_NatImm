library(stringr)
library(ggplot2)
library(ggridges)
library(ggExtra)
library(psych)
library(survival)
library(gridExtra)
# List files

c_types = c("LUAD")
full_names = c("Lung Adenocarcinoma")
names(full_names)  = c_types

gene_symbols = NULL

c_type = c_types[1]


c_path_mat = file.path("C:/Projects/CXCR6/TCGA", c_type, paste0("TCGA-",c_type,".htseq_fpkm.tsv.gz"))
mrna = read.csv(gzfile(c_path_mat), sep='\t', header = T, row.names = 1)

# Read survival data
c_path_surv = file.path("C:/Projects/CXCR6/TCGA", c_type,paste0("TCGA-",c_type,".survival.tsv"))
surv = read.csv(c_path_surv, sep="\t", header = T, row.names = 1)

head(surv)

head(colnames(mrna))
head(rownames(mrna))
mrna[1:5,1:4]

colnames(mrna) = gsub(pattern = '[.]', '-', colnames(mrna))


# clean gene ids 
rownames(mrna) <- str_replace(rownames(mrna),
                              pattern = ".[0-9]+$",
                              replacement = "")

# Check only in tumor patients from TCGA barcode, Sample type decoded in  4th field
# Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.

tumor.ids = which(substring(colnames(mrna), 14, 14) == "0") 
surv.tuumor.ids =  which(substring(rownames(surv), 14, 14) == "0")

intersect.ids = intersect(colnames(mrna)[tumor.ids], 
                          rownames(surv)[surv.tuumor.ids])



# ENSG00000172215 = CXCR6
if(is.null(gene_symbols)){
  require("EnsDb.Hsapiens.v86")
  
  options(ensembldb.seqnameNotFound = "NA")
  edb <- EnsDb.Hsapiens.v86
  gene_symbols = mapIds(edb,  keys= rownames(mrna), keytype = "GENEID", column = "SYMBOL")
}

# TRM/ CXCR6 path
gene.path = "C:/Projects/external/jem_20190249_tables1/"

lungTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T)
lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[lungTRM_noTRM$`log2 fold change`>2 & lungTRM_noTRM$`P adj`<0.01]
# lungTRM_noTRM.genes = lungTRM_noTRM$`Gene ID`[order(lungTRM_noTRM$`log2 fold change`, decreasing = T)[1:10]]


tumorTRM_noTRM = readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS4.xlsx"),skip = 2,col_names = T)
tumorTRM_noTRM.genes = tumorTRM_noTRM$`Gene ID`[tumorTRM_noTRM$`log2 fold change`>2 & tumorTRM_noTRM$`P adj`<0.01]
# tumorTRM_noTRM.genes = tumorTRM_noTRM$`Gene ID`[order(tumorTRM_noTRM$`log2 fold change`, decreasing = T)[1:10]]

tumorTRM_noTRM.genes = setdiff(tumorTRM_noTRM.genes, "CXCR6")
lungTRM_noTRM.genes = setdiff(lungTRM_noTRM.genes, "CXCR6")



# TRM.markers = intersect(lungTRM_noTRM.genes, tumorTRM_noTRM.genes)
# TRM.markers= union(conserved_TRM.genes, c("RBPJ", "ITGAE","ZNF683"))

trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS3.xlsx"),skip = 2,col_names = T,sheet = 2)
trm_table = trm_table[!is.na(trm_table$`Log2 fold change in human lung`),]
table(trm_table$`Log2 fold change in human lung`>0)
TRM.skin = trm_table$`Cheuk, et al. 2017 human skin TRM signature`[trm_table$`Log2 fold change in human lung`>0]
# TRM.markers= c("RBPJ", "ITGAE","ZNF683")

tumorTRM_noTRM.genes.ex = setdiff(tumorTRM_noTRM.genes, lungTRM_noTRM.genes)
lungTRM_noTRM.genes.ex = setdiff(lungTRM_noTRM.genes, tumorTRM_noTRM.genes)
TRM.markers = setdiff(intersect(tumorTRM_noTRM.genes, lungTRM_noTRM.genes), "CXCR6")
TRM.skin.ex = setdiff(TRM.skin, union(tumorTRM_noTRM.genes.ex,lungTRM_noTRM.genes.ex ))
TRM.markers_union = setdiff(union(tumorTRM_noTRM.genes, lungTRM_noTRM.genes), "CXCR6")


getSignature<-function(fpkm,gene.list){
  require("EnsDb.Hsapiens.v86")
  options(ensembldb.seqnameNotFound = "NA")
  edb <- EnsDb.Hsapiens.v86
  
  ens.idx = grep("ENSG", gene.list)
  non.ens.idx = grep("ENSG", gene.list, invert = T)
  
  gene_symbols = c()
  if(length(non.ens.idx)>0)
    gene_symbols = mapIds(edb,  keys= gene.list[non.ens.idx], column = "GENEID", keytype = "SYMBOL")
  if(length(ens.idx)>0)
    gene_symbols = c(gene_symbols, gene.list[ens.idx])
  
  notfound = anyNA(gene_symbols)| anyNA(match(gene_symbols,rownames(fpkm)))
  if(notfound){
    # warning('Not found')
    warning('Not found: ', paste(gene.list[is.na(match(gene_symbols,rownames(fpkm)))],collapse=", "))
  }
    
  
  gene_symbols = gene_symbols[!is.na(gene_symbols)]
  mat_interest = fpkm[gene_symbols, , drop = F] 
  mat_interest = log2(mat_interest[complete.cases(mat_interest),]+1)
  mat_interest = mat_interest[rowSums(mat_interest)>0, ]
  signature = apply(scale(t(mat_interest), center = F),1, mean)
  return(signature)
}


tumorTRM_noTRM.signature = getSignature(mrna[,intersect.ids], tumorTRM_noTRM.genes.ex)
lungTRM_noTRM.signature = getSignature(mrna[,intersect.ids],lungTRM_noTRM.genes.ex)
dermisTRM.signature = getSignature(mrna[,intersect.ids],TRM.skin.ex)
commonTRM.signature = getSignature(mrna[,intersect.ids],TRM.markers)
unionTRM.signature = getSignature(mrna[, intersect.ids], TRM.markers_union)
CXCR6.signature_pseudo = getSignature(mrna[, intersect.ids], c("CXCR6"))
CCR2.signature_pseudo = getSignature(mrna[, intersect.ids], c("CCR2"))

# CXCL16 correlation with CD64, CD14, and CD68 expression 
CXCL16.signature_pseudo = getSignature(mrna[, intersect.ids], c("CXCL16"))
CD64.signature_pseudo = getSignature(mrna[, intersect.ids], c("FCGR1A")) # CD64
CD14.signature_pseudo = getSignature(mrna[, intersect.ids], c("CD14"))
CD68.signature_pseudo = getSignature(mrna[, intersect.ids], c("CD68"))
CCL2.signature_pseudo = getSignature(mrna[, intersect.ids], c("CCL2"))
CCL6.signature_pseudo = getSignature(mrna[, intersect.ids], c("CCL6"))

CCL7.signature_pseudo = getSignature(mrna[, intersect.ids], c("CCL7"))
CCL15.signature_pseudo = getSignature(mrna[, intersect.ids], c("ENSG00000275718"))
CCL23.signature_pseudo = getSignature(mrna[, intersect.ids], c("CCL23"))
CCL.signature = getSignature(mrna[, intersect.ids], c("CCL2", "CCL7", "CCL23", "ENSG00000275718"))


#  "CCL2", "CCL23", "CCL7", "CCL9", "CCL15","FCGR1A", "CD14", "CD68"


df.signatures = data.frame(tumor.lungTRM = tumorTRM_noTRM.signature,
                           # CXCR6 = CXCR6.signature_pseudo,
                           # CCR2 = CCR2.signature_pseudo,
                           CXCL16 = CXCL16.signature_pseudo,
                           CD64 = CD64.signature_pseudo,
                           CD14 = CD14.signature_pseudo,
                           CD68 = CD68.signature_pseudo,
                           CCL2 = CCL2.signature_pseudo,
                           CCL7 = CCL7.signature_pseudo,
                           # CCL15 = CCL15.signature_pseudo,
                           CCL23 = CCL23.signature_pseudo, 
                           CCL.signature = CCL.signature)

pairs.panels(df.signatures,
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pdf(file.path("C:/Projects/CXCR6/manuscript_ready/TRM_CCL_panel.pdf"), width = 8, height = 8)

pairs.panels(df.signatures[,c(1:2,9)],
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
pairs.panels(df.signatures[,c(1:2,6:9)],
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pairs.panels(df.signatures[,2:5],
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB", 
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pairs.panels(df.signatures[,1:8],
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB", 
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


pairs.panels(df.signatures[,c(1:5,9)],
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB", 
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

dev.off()

pdf(file.path("C:/Projects/CXCR6/manuscript_ready/TRM_CCL_panel_7.pdf"), width = 8, height = 8)
pairs.panels(df.signatures[,1:7],
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB", 
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

pairs(

## grading survivals


surv.data <- surv[intersect.ids,c("OS.time","OS")] # survival data
mySurv <- Surv(surv.data$OS.time, surv.data$OS, type='right')
print(paste(full_names[c_type], "Total events:", length(which(mySurv[,"status"]==1))))


#define low both CXCR6/TRM 
QT_LOW = 0.333
QT_HI = 0.666

# low  = intersecting low of different signatures
# high  = intersecting high of above signatures
plotSurvSig<-function(sigName = NULL){
  
  colID = grep(paste(sigName, collapse = "|"),
               colnames(df.signatures), ignore.case = T)
  
  
  
  low = list()
  high = list()
  for(i in colID){
    low[[paste0(i)]] = which(df.signatures[,i] < quantile(df.signatures[,i], QT_LOW))
    high[[paste0(i)]] = which(df.signatures[,i] > quantile(df.signatures[,i], QT_HI))
  }
  low = Reduce(intersect, low)
  high = Reduce(intersect, high)
  
  df <- data.frame(survival.time = surv.data$OS.time, censor =surv.data$OS)
  head(df)
  df$group= NA
  df$group[low] = paste0(colnames(df.signatures)[colID],"lo", collapse = "_")
  df$group[high] = paste0(colnames(df.signatures)[colID],"hi", collapse = "_")
  
  # df = df[complete.cases(df),]
  df = df[df$survival.time<3650,]
  
  df$group = factor(df$group, levels = sort(unique(df$group)))
  table(df$group)
  
  library(survminer)
  fit =  survfit(Surv(survival.time, censor) ~ group, data = df)
  # coxfit <- coxph(Surv(survival.time, censor) ~ TRM.grade, data = df, ties = 'exact')
  pg<-ggsurvplot(fit, data = df,
                 title=paste(colnames(df.signatures)[colID], collapse = "_"),xlab = c("Days"),
                 size = 0.5,                 # change line size
                 legend = c(0.75,0.85),          # legend position topright
                 palette =
                   c("red", "#2E9FDF"),# custom color palettes
                 # conf.int = TRUE,          # Add confidence interval
                 pval = T,              # Add p-value
                 # risk.table = TRUE,        # Add risk table
                 risk.table.col = "strata",# Risk table color by groups
                 surv.median.line = "hv",combine = T,
                 # test.for.trend = T,
                 legend.labs = levels(df$group),    # Change legend labels
                 risk.table.height = 0.25, # Useful to change when you have multiple groups
                 # ggtheme = theme_bw()      # Change ggplot2 theme
  )$plot
  # return(grid.arrange(ggplotGrob(pg$plot)))
  # p2<- grid.arrange(ggplotGrob(p[[2]]$plot))
  # 
  # return(coxfit)
  return(pg)
}



p = list()
# 
# p[[1]] = ggpubr::ggpar(plotSurvSig(c("CXCR6$")),font.legend = list(size = 9, color = "black"))
# p[[3]] = ggpubr::ggpar(plotSurvSig(c("^TRM")),font.legend = list(size = 9, color = "black"))
# p[[2]] = ggpubr::ggpar(plotSurvSig(c("CXCR6$","^lung")),font.legend = list(size = 9, color = "black"))
# p[[4]] = ggpubr::ggpar(plotSurvSig(c("CXCR6$","^tumor")),font.legend = list(size = 9, color = "black"))
# # p[[5]] = plotSurvSig(c("CXCR6_lung"))
# p[[6]] = plotSurvSig(c("CXCR6_tumor"))


p[[1]] = plotSurvSig(c("CXCL16$", "CCL2$"))
p[[2]] = plotSurvSig(c("CXCL16$", "CCL7"))
p[[3]] = plotSurvSig(c("CXCL16$", "CCL23"))
# p[[4]] = plotSurvSig(c("CXCL16$", "CCL15"))
p[[4]] = plotSurvSig(c("CXCL16$", "^CCL.sig"))

pdf(file.path("C:/Projects/CXCR6/manuscript_ready/CXCL16_CCL_surv.pdf"), width = 7, height = 6)
p[[1]]
p[[2]]
p[[3]]
p[[4]]
dev.off()

pair_hi<-function(sigName){
  colID = grep(paste(sigName, collapse = "|"),
               colnames(df.signatures), ignore.case = T)
  
  
  QT_HI = 0.66
  # low = list()
  high = list()
  for(i in colID){
    # low[[paste0(i)]] = which(df.signatures[,i] < quantile(df.signatures[,i], QT_LOW))
    high[[paste0(i)]] = which(df.signatures[,i] > quantile(df.signatures[,i], QT_HI))
  }
  low = Reduce(intersect, low)
  high = Reduce(intersect, high)
  
  hi = list()
  hi[[colnames(df.signatures)[colID[1]]]] = df.signatures[high, colnames(df.signatures)[colID[1]]]
  hi[[colnames(df.signatures)[colID[2]]]]  = df.signatures[high,  colnames(df.signatures)[colID[2]]]
  # df.signatures.hi = data.frame(CCL2.hi = CCL2.hi, CXCL16.hi = CXCL16.hi)
  df.signatures.hi= do.call(cbind, hi)
  colnames(df.signatures.hi) = paste0(colnames(df.signatures.hi) ,"_Hi")
  return(df.signatures.hi)
}

pdf(file.path("C:/Projects/CXCR6/manuscript_ready/CXCL16_CCL_corr_HI.pdf"), width = 8, height = 8)

pairs.panels(pair_hi(c("CXCL16$",  "CCL2")),
             stars = T, cex.cor = 0.8,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


pairs.panels(pair_hi(c("CXCL16$",  "CCL7")),
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pairs.panels(pair_hi(c("CXCL16$",  "CCL23")),
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pairs.panels(pair_hi(c("CXCL16$",  "^CCL.signature")),
             stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             main = paste(full_names[1], paste0('(',c_type,')')),
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

dev.off()
