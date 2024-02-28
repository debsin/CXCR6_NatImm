
# List files

dir = "C:/Projects/TCGA/"
setwd(dir)

tcga.path = "luad/denseDataOnlyDownload.tsv"

  
tcga.df = read.table(tcga.path, row.names = 1, header = T, sep = "\t")
tcga.df = tcga.df[,-1]
tcga.df = tcga.df[complete.cases(tcga.df),]

tumor.tcga = which(substring(rownames(tcga.df), 14, 14) %in% c("0")) 
tcga.df = tcga.df[tumor.tcga,]

plot(density(scale(tcga.df$CXCR6, scale = F)))
     

tcga.df$CXCR6 = scale(tcga.df$CXCR6)
shapiro.test(tcga.df$CXCR6)

cutoff = median(tcga.df$CXCR6)
hist(tcga.df$CXCR6, probability = T, breaks = 50)
lines(density(tcga.df$CXCR6,), col = "red", lwd = 2)
abline(v = cutoff,col="blue", lwd = 2, lty = "dashed")

tcga.df$group = as.factor(ifelse(tcga.df$CXCR6 > cutoff, "high", "low"))

library("survival")
library("survminer")
library(ggplot2)
library(stringr)


fit <- survfit(Surv(OS.time, OS) ~ group, data = tcga.df)
print(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))



# Read mRNA expression data
c_type="LUAD"


fpkm.df = read.csv("luad/gdc_fpkm_luad.tsv", sep='\t', header = T, row.names = 1)
fpkm.df = fpkm.df[complete.cases(fpkm.df),]

# colnames(mrna) = gsub(pattern = '[.]', '-', colnames(mrna))

tumor.ids = which(substring(rownames(fpkm.df), 14, 14) %in% c("0")) 


# Test for normality


scale.genes = as.data.frame(cbind(FPKM = fpkm.df$ENSG00000172215.5[tumor.ids],
                                  log2_FPKM = log2(fpkm.df$ENSG00000172215.5[tumor.ids]+1)))

par(mfrow = c(2,2))

qqnorm(scale.genes$FPKM,xlab = "FPKM", main = "CXCR6")
qqline(scale.genes$FPKM)

qqnorm(scale.genes$log2_FPKM, xlab = "log2(FPKM+1)", main = "CXCR6")
qqline(scale.genes$log2_FPKM)

sh = shapiro.test(scale.genes$log2_FPKM)

boxplot(scale.genes)

hist(scale.genes$log2_FPKM, 
     main=paste("Shapiro test: stat",round(sh$statistic, 3), "pval",round(sh$p.value,3)), 
     xlab="log2(FPKM+1)", 
     border="light blue", 
     col="blue", 
     las=1, 
     breaks=50)

par(mfrow = c(1,1))

# clean gene ids 

colnames(fpkm.df) <- str_replace(colnames(fpkm.df),
                              pattern = ".[0-9]+$",
                              replacement = "")

map_ids = readRDS("mapIDs.Rds")
genes_interest = c("CXCR6")
get_rows = na.omit(map_ids$ensembl_gene_id[match(genes_interest, map_ids$hgnc_symbol)])
if(length(get_rows) != length(genes_interest)) stop("Do something! unexpected number of genes retrived.")




intersect.ids = tumor.ids

surv.data <- fpkm.df[intersect.ids,c("OS.time","OS")] # survival data
mySurv <- Surv(fpkm.df$OS.time, fpkm.df$OS, type='right')
print(paste("Total events:", length(which(mySurv[,"status"]==1))))

check.genes = as.data.frame(fpkm.df[intersect.ids, get_rows])
head(check.genes)
colnames(check.genes) = genes_interest

# check.genes = scale(check.genes,center = T,scale = T)
check.genes = log2(check.genes+1)


df <- data.frame(survival.time = surv.data$OS.time, censor =surv.data$OS, check.genes)
head(df)

hist(df$CXCR6, probability = T, breaks = 50)
lines(density(df$CXCR6,), col = "red", lwd = 2)

# density to find reasonable cutoff
itp1 = median(df$CXCR6)- 0.07
abline(v=itp1, col="blue", lwd = 2, lty = "dashed")



## Some processing to get the staging right
df$CXCR6.grade <- df$CXCR6
df$CXCR6.grade <- "high"
df$CXCR6.grade[df$CXCR6 < itp1] <- "low"
table(df$CXCR6.grade)




plot_list = list()


plot_list[[1]] = ggplot(df, aes(x=CXCR6)) + geom_histogram(aes(y=..density..), bins = 50, colour="black", fill="white")+
  geom_density(alpha=.2, fill="brown", size=0.7) + geom_vline(aes(xintercept=itp1),
                                                              color="forestgreen", linetype="longdash", size=1)+ theme_bw()+
  ggtitle("LUAD",subtitle = "")

# Drawing survival curves
fit =  survfit(Surv(survival.time, censor) ~ CXCR6.grade, data = df)
plot_list[[2]] = ggsurvplot(fit, data = df,  title="CXCR6",
                            size = 0.5,                 # change line size
                            legend = c(0.85,0.85),          # legend position topright 
                            palette = 
                              c("#E7B800", "#2E9FDF"),# custom color palettes
                            conf.int = TRUE,          # Add confidence interval
                            pval = TRUE,              # Add p-value
                            risk.table = TRUE,        # Add risk table
                            risk.table.col = "strata",# Risk table color by groups
                            surv.median.line = "v",
                            legend.labs = 
                              c("high", "low"),    # Change legend labels
                            risk.table.height = 0.25, # Useful to change when you have multiple groups
                            ggtheme = theme_bw()      # Change ggplot2 theme
)


gridExtra::grid.arrange(grobs=list(plot_list[[1]], plot_list[[2]]$plot, plot_list[[2]]$table), heights = c(2,3,1))


