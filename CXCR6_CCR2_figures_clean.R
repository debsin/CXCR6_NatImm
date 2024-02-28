library(stringr)
library(ggplot2)
library(ggridges)
library(ggExtra)
library(psych)
library(survival)
library(gridExtra)
# List files

c_types = c("LUAD","LIHC", "BRCA", "STAD", "LAML", "BLCA", "HNSC", "SKCM", "COAD", "READ",  "PRAD")
full_names = c("Lung Adenocarcinoma","Liver Cancer", "Breast Cancer", "Stomach Cancer","Acute myeloid leukemia",
               "Bladder Cancer", "Head and Neck Cancer", "Skin Cutaneous Melanoma",
               "Colon Cancer", "Rectal Cancer",  "Prostate Cancer")
names(full_names)  = c_types

gene_symbols = NULL

pdf(file.path("C:/Projects/TCGA/CCR2_CXCR6.pdf"), width = 8, height = 7)

# c_type = c_types[10]
for(c_type in c_types){
  dir = file.path("C:/Projects/TCGA/", c_type)
  setwd(dir)
  
  c_path_mat = paste0("TCGA-",c_type,".htseq_fpkm.tsv.gz")
  mrna = read.csv(gzfile(c_path_mat), sep='\t', header = T, row.names = 1)
  
  # Read survival data
  c_path_surv = paste0("TCGA-",c_type,".survival.tsv")
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
  
  conserved_TRM.genes = intersect(lungTRM_noTRM.genes, tumorTRM_noTRM.genes)
  
  
  tumorTRM_noTRM.genes = setdiff(tumorTRM_noTRM.genes, "CXCR6")
  lungTRM_noTRM.genes = setdiff(lungTRM_noTRM.genes, "CXCR6")
  
  
  # TRM.markers = intersect(lungTRM_noTRM.genes, tumorTRM_noTRM.genes)
  # TRM.markers= union(conserved_TRM.genes, c("RBPJ", "ITGAE","ZNF683"))
  
  trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS3.xlsx"),skip = 2,col_names = T,sheet = 2)
  trm_table = trm_table[!is.na(trm_table$`Log2 fold change in human lung`),]
  table(trm_table$`Log2 fold change in human lung`>0)
  TRM.markers = trm_table$`Cheuk, et al. 2017 human skin TRM signature`[trm_table$`Log2 fold change in human lung`>0]
  # TRM.markers= c("RBPJ", "ITGAE","ZNF683")
  
  tumorTRM_noTRM.genes.ex = setdiff(tumorTRM_noTRM.genes, lungTRM_noTRM.genes)
  lungTRM_noTRM.genes.ex = setdiff(lungTRM_noTRM.genes, tumorTRM_noTRM.genes)
  TRM.markers = setdiff(TRM.markers, union(tumorTRM_noTRM.genes.ex,lungTRM_noTRM.genes.ex ))
  
  getSignature<-function(fpkm,gene.list){
    require("EnsDb.Hsapiens.v86")
    options(ensembldb.seqnameNotFound = "NA")
    edb <- EnsDb.Hsapiens.v86
    gene_symbols = mapIds(edb,  keys= gene.list, column = "GENEID", keytype = "SYMBOL")
    
    notfound = is.na(gene_symbols)
    mat_interest = fpkm[gene_symbols, , drop = F] 
    mat_interest = log2(mat_interest[complete.cases(mat_interest),]+1)
    mat_interest = mat_interest[rowSums(mat_interest)>0, ]
    signature = apply(scale(t(mat_interest), center = F),1, mean)
    return(signature)
    
  }
  
  
  tumorTRM_noTRM.signature = getSignature(mrna[,intersect.ids], 
                                          setdiff(tumorTRM_noTRM.genes.ex, "CXCR6"))
  lungTRM_noTRM.signature = getSignature(mrna[,intersect.ids],
                                         setdiff(lungTRM_noTRM.genes.ex, "CXCR6"))
  TRM.signature = getSignature(mrna[,intersect.ids],
                               setdiff(TRM.markers, "CXCR6"))
  
  
  cor(tumorTRM_noTRM.signature, lungTRM_noTRM.signature)
  
  df.signatures = data.frame(samples = intersect.ids,
                             CXCR6 = getSignature(mrna[, intersect.ids], c("CXCR6")),
                             CCR2 = getSignature(mrna[, intersect.ids], c("CCR2")),
                             # tumorTRM = tumorTRM_noTRM.signature,
                             # lungTRM = lungTRM_noTRM.signature,
                             TRM = TRM.signature
                             # CXCR6_tumorTRM = getSignature(mrna[, intersect.ids], 
                             #                               union(tumorTRM_noTRM.genes.ex, "CXCR6")),
                             # CXCR6_lungTRM = getSignature(mrna[, intersect.ids], 
                             #                               union(lungTRM_noTRM.genes.ex, "CXCR6")),
                             # CXCR6_TRM = getSignature(mrna[, intersect.ids], 
                             #                          union(TRM.markers, "CXCR6"))
                             # gDT = gdt.signature,
                             # CD8TMem = cd8Mem.signature,
                             # CD4TMem = cd4Mem.signature
  )
  
  
  head(df.signatures)
  
  sig.plot.df = reshape::melt(df.signatures)
  
  head(sig.plot.df)
  
  
  # ggplot(sig.plot.df, aes(x = samples, y = value, group = variable))+geom_line(aes(col = variable))
  
  # pairs(df.signatures[,-1], pch = 16)
  
  
  
  cor(df.signatures[,-1])
  
  
  
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
  
  
  p[[1]] = plotSurvSig(c("CXCR6$"))
  # p[[2]] = plotSurvSig(c("^TRM"))
  p[[2]] = plotSurvSig(c("^CCR2"))
  p[[3]] = plotSurvSig(c("CCR2$", "CXCR6$"))
  # p[[2]] = ggpubr::ggpar(plotSurvSig(c("CXCR6$","^lung")),font.legend = list(size = 9, color = "black"))
  # p[[4]] = ggpubr::ggpar(plotSurvSig(c("CXCR6$","^tumor")),font.legend = list(size = 9, color = "black"))
  # 
  
  
  # pdf(file.path(dir ,paste0(c_type,"_CXCR6.pdf"), width = 11, height = 10)
  df.signature = df.signatures[intersect.ids,-1]
  sp <- ggscatter(df.signature, x = "CCR2", y = "CXCR6",size = 0.5,title  = full_names[c_type],
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
  )+font("title", size = 18, face = "bold.italic")
  # Add correlation coefficient
  cor.plot = sp + stat_cor(method = "spearman", cor.coef.name = "rho")+
    xlab("CCR2 expression")+
    ylab("CXCR6 expression")
  
  mar.plot = ggMarginal(cor.plot, type = "density",fill = "#00AFBB88")
  
  
  
  lay <- rbind(c(1,1,1,4,4,4),
               c(1,1,1,4,4,4),
               c(1,1,1,4,4,4),
               c(2,2,2,3,3,3),
               c(2,2,2,3,3,3),
               c(2,2,2,3,3,3))
  grid.arrange(mar.plot,ggplotGrob(p[[3]]),
               ggplotGrob(p[[1]]),
               ggplotGrob(p[[2]]) ,ncol = 2)
  
  
  
  
  # survp = arrange_ggsurvplots(p,ncol = 2, nrow= 2, byrow=FALSE,print = F)
  #   pairs.panels(df.signatures[,-1],stars = T,show.points = T,
  #                method = "pearson", # correlation method
  #                hist.col = "grey", main = full_names[i],cex.cor = 0.4, breaks = 20,
  #                density = TRUE,  # show density plots
  #                ellipses = F # show correlation ellipses 
  #   )
  #   
  
  
  # dev.off()
  
}

dev.off()
