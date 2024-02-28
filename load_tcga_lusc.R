
# List files

dir = "C:/Projects/TCGA/lusc/"
setwd(dir)


rna = read.table(gzfile("TCGA.LUSC.sampleMap_HiSeqV2.gz"), sep = "\t", header = T, row.names = 1)

os = read.table("survival_LUSC_survival.txt", sep = "\t", header = T, row.names = 1)

table(substr(rownames(os),14,15))

library(survival)
plot(survfit(Surv(os$OS.time, os$OS) ~ 1, data = rna), 
     xlab = "Days", 
     ylab = "Overall survival probability")

my.data = os

my.data$group = ifelse((os$OS.time < median(os$OS.time)), "low", "high")
km_AG_fit <- survfit(Surv(os$OS.time, os$OS) ~ 1, data=rna["CXCR6",, drop = F])
autoplot(km_AG_fit)




survdiff(Surv(time, status) ~ sex, data = rna["CXCR6",])


# read RNA file 
rna <- read.table('LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', header=T,row.names=1,sep='\t')
# and take off first row cause we don't need it
rna <- rna[-1,]



# and read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
clinical <- as.data.frame(t(read.delim('LUSC.merged_only_clinical_clin_format.txt', header = T, row.names = 1, sep='\t')))



# first I remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]

# see the values
table(substr(colnames(rna),14,14))

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == '1')
t_index <- which(substr(colnames(rna),14,14) == '0')


# apply voom function from limma package to normalize the data
library(limma)
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}


rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub('\\.','-',substr(colnames(rna),1,12))

# and check how data look, they should look normally-ish distributed
hist(log(t(apply(rna,1,as.numeric))))
hist(rna_vm)

# we can remove the old "rna" cause we don't need it anymore
rm(rna)

# z = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)


# calculate z-scores
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])
# set the rownames keeping only gene name
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]])



# match the patient ID in clinical data with the colnames of z_rna
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
sum(clinical$IDs %in% colnames(z_rna)) # we have 529 patients that we could use

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(clinical))

# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}

# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')


# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}


# create vector for death censoring
table(clinical$patient.vital_status)
# alive dead
# 343   161

all_clin$death_event <- ifelse(clinical$patient.vital_status == 'alive', 0,1)

#finally add row.names to clinical
rownames(all_clin) <- clinical$IDs



# create event vector for RNASeq data
event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))

# since we need the same number of patients in both clinical and RNASeq data take the indices for the matching samples
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(all_clin))
ind_clin <- which(rownames(all_clin) %in% colnames(z_rna))

# pick your gene of interest
ind_gene <- which(rownames(z_rna) == 'CXCR6')

# check how many altered samples we have
table(event_rna[ind_gene,])

# run survival analysis
library(survival)
s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

# extraect the p.value
pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

# plot the data
plot(survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]),
     col=c(1:3), frame=F, lwd=2,main=paste('LUSC',rownames(z_rna)[ind_gene],sep='\n'))

# add lines for the median survival
x1 <- ifelse ( is.na(as.numeric(summary(s)$table[,'median'][1])),'NA',as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if(x1 != 'NA' & x2 != 'NA'){
   lines(c(0,x1),c(0.5,0.5),col='blue')
   lines(c(x1,x1),c(0,0.5),col='black')
   lines(c(x2,x2),c(0,0.5),col='red')
}

# add legend
legend(1800,0.995,legend=paste('p.value = ',pv[[1]],sep=''),bty='n',cex=1.4)
legend(max(as.numeric(as.character(all_clin$death_days)[ind_clin]),na.rm = T)*0.7,0.94,
       legend=c(paste('NotAltered=',x1),paste('Altered=',x2)),bty='n',cex=1.3,lwd=3,col=c('black','red'))





















sample.path= "gdc_sample_sheet.2021-09-09.tsv"

exp.path = "gexp"

clinical.path  = "clinical.cart.2021-09-09.tar.gz"



samplesheet = read.csv(file.path(dir,sample.path), sep = "\t", header = T)

untar(tarfile  = file.path(dir, clinical.path), files = "clinical.tsv",exdir = dir)

# read each file and concatenate
fpkm.list = list()
for(i in 1:nrow(samplesheet)){
  file.names = file.path(dir, exp.path, samplesheet$File.ID[i], samplesheet$File.Name[i])
  fpkm.list[[samplesheet$Sample.ID[i]]] = read.table(gzfile(file.names), row.names = 1)
  colnames(fpkm.list[[samplesheet$Sample.ID[i]]]) = samplesheet$Sample.ID[i]
}

fpkm  = do.call(cbind, fpkm.list)
rm(fpkm.list)

gene_query = gsub(rownames(fpkm),pattern = ".[0-9]+$", replacement = "")
rownames(fpkm) = gene_query

# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# gene_symbols = getBM(filters = "ensembl_gene_id", uniqueRows = F,
#                 attributes = c("hgnc_symbol","ensembl_gene_id"),
#                 values = gene_query, mart = mart)

library("EnsDb.Hsapiens.v86")
options(ensembldb.seqnameNotFound = "NA")
edb <- EnsDb.Hsapiens.v86
gene_symbols = mapIds(edb,  keys= gene_query, column = "SYMBOL", keytype = "GENEID")

notfound = is.na(gene_symbols)
rownames(fpkm)[!notfound] = make.unique(gene_symbols[!notfound])


mapIds(edb,  keys= "CXCR6", column = "GENEID", keytype = "SYMBOL")

grep("ENSG00000172215", gene_query)
rownames(fpkm)[28048]


# read clinical data

meta.data = read.csv(file.path(dir, "clinical.tsv"), sep = "\t", header = T)
meta.data$SampleID = samplesheet$Sample.ID[match(meta.data$case_id, samplesheet$File.ID)]
# days_to_death days_to_last_follow_up vital_status

