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
  gene_symbols = mapIds(edb,  keys= gene.list, column = "GENEID", keytype = "SYMBOL")
  
  notfound = is.na(gene_symbols)
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

## SUPPLEMENT Figure on LUAD and TRM
df.signatures = data.frame(CXCR6 = CXCR6.signature_pseudo,
                           tumor.lungTRM = tumorTRM_noTRM.signature,
                           lungTRM = lungTRM_noTRM.signature,
                           intersect.TRM = commonTRM.signature,
                           union.TRM = unionTRM.signature,
                           dermisTRM = dermisTRM.signature)


### Main figure on  CXCR6 and CCL2, etc
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


df.signatures2 = data.frame(tumor.lungTRM = tumorTRM_noTRM.signature,
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


pairs.panels2<-function (x, smooth = TRUE, scale = FALSE, density = TRUE, ellipses = TRUE,
                         digits = 2, method = "pearson", pch = 20, lm = FALSE, cor = TRUE, 
                         jiggle = FALSE, factor = 2, hist.col = "cyan", show.points = TRUE, 
                         rug = TRUE, breaks = "Sturges", cex.cor = 1, wt = NULL, smoother = FALSE, 
                         stars = FALSE, ci = FALSE, alpha = 0.05, ...) 
{
  "panel.hist.density" <- function(x, ...) {
    usr <- par("usr")
    on.exit(par("usr"))
    par(usr = c(usr[1], usr[2], 0, 1.5))
    tax <- table(x)
    if (length(tax) < 11) {
      breaks <- as.numeric(names(tax))
      y <- tax/max(tax)
      interbreak <- min(diff(breaks)) * (length(tax) - 
                                           1)/41
      rect(breaks - interbreak, 0, breaks + interbreak, 
           y, col = hist.col)
    }
    else {
      h <- hist(x, breaks = breaks, plot = FALSE)
      breaks <- h$breaks
      nB <- length(breaks)
      y <- h$counts
      y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col = hist.col)
    }
    if (density) {
      tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", 
                               adjust = 1.2), silent = TRUE)
      if (!inherits(tryd, "try-error")) {
        d$y <- d$y/max(d$y)
        lines(d)
      }
    }
    if (rug) 
      rug(x)
  }
  "panel.cor" <- function(x, y, prefix = "", ...) {
    usr <- par("usr")
    on.exit(par("usr"))
    par(usr = c(0, 1, 0, 1))
    if (is.null(wt)) {
      r <- cor(x, y, use = "pairwise", method = method)
    }
    else {
      r <- cor.wt(data.frame(x, y), w = wt[, c(1:2)])$r[1, 
                                                        2]
    }
    txt <- format(c(round(r, digits), 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (stars) {
      # pval <- r.test(sum(!is.na(x * y)), r)$p
      # symp <- symnum(pval, corr = FALSE, cutpoints = c(0, 
      #                                                  0.001, 0.01, 0.05, 1), symbols = c("***", "**", 
      #                                                                                     "*", " "), legend = FALSE)
      pval = cor.test(x,y)$p.value
      # pval <- r.test(sum(!is.na(x * y)), r)$p
      pval <- signif(pval, 2)
      
    }
    cex <- cex.cor * 0.4/(max(strwidth("0.12***"), strwidth(txt)))
    pval= paste0("p=", pval)
    if (scale) {
      cex1 <- cex * abs(r)
      if (cex1 < 0.2) 
        cex1 <- 0.2
      text(0.5, 0.5,  bquote(.(txt)^.(pval)), cex = cex1)
    }
    else {
      text(0.5, 0.5,  bquote(.(txt)^.(pval)), cex = cex)
    }
  }
  "panel.smoother" <- function(x, y, pch = par("pch"), col.smooth = "red", 
                               span = 2/3, iter = 3, ...) {
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    r = cor(x, y, use = "pairwise", method = method)
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (smoother) {
      smoothScatter(x, y, add = TRUE, nrpoints = 0)
    }
    else {
      if (show.points) 
        points(x, y, cex = 0.5, pch = pch, ...)
    }
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      if (smooth & ci) {
        lml <- loess(y ~ x, degree = 1, family = "symmetric")
        tempx <- data.frame(x = seq(min(x, na.rm = TRUE), 
                                    max(x, na.rm = TRUE), length.out = 47))
        pred <- predict(lml, newdata = tempx, se = TRUE)
        if (ci) {
          upperci <- pred$fit + confid * pred$se.fit
          lowerci <- pred$fit - confid * pred$se.fit
          polygon(c(tempx$x, rev(tempx$x)), c(lowerci, 
                                              rev(upperci)), col = adjustcolor("light grey", 
                                                                               alpha.f = 0.8), border = NA)
        }
        lines(tempx$x, pred$fit, col = col.smooth, ...)
      }
      else {
        if (smooth) 
          lines(stats::lowess(x[ok], y[ok], f = span, 
                              iter = iter), col = col.smooth)
      }
    }
    if (ellipses) 
      draw.ellipse(xm, ym, xs, ys, r, col.smooth = col.smooth, 
                   ...)
  }
  "panel.lm" <- function(x, y, pch = par("pch"), col.lm = "red", 
                         ...) {
    ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin, xmin), max(ymax, xmax))
    xlim <- ylim
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (smoother) {
      smoothScatter(x, y, add = TRUE, nrpoints = 0)
    }
    else {
      if (show.points) {
        points(x, y, pch = pch, ylim = ylim, xlim = xlim, cex = 0.5, 
               ...)
      }
    }
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      lml <- lm(y ~ x)
      if (ci) {
        tempx <- data.frame(x = seq(min(x, na.rm = TRUE), 
                                    max(x, na.rm = TRUE), length.out = 47))
        pred <- predict.lm(lml, newdata = tempx, se.fit = TRUE)
        upperci <- pred$fit + confid * pred$se.fit
        lowerci <- pred$fit - confid * pred$se.fit
        polygon(c(tempx$x, rev(tempx$x)), c(lowerci, 
                                            rev(upperci)), col = adjustcolor("light grey", 
                                                                             alpha.f = 0.8), border = NA)
      }
      if (ellipses) {
        xm <- mean(x, na.rm = TRUE)
        ym <- mean(y, na.rm = TRUE)
        xs <- sd(x, na.rm = TRUE)
        ys <- sd(y, na.rm = TRUE)
        r = cor(x, y, use = "pairwise", method = method)
        draw.ellipse(xm, ym, xs, ys, r, col.smooth = col.lm, 
                     ...)
      }
      abline(lml, col = col.lm, ...)
    }
  }
  "draw.ellipse" <- function(x = 0, y = 0, xs = 1, ys = 1, 
                             r = 0, col.smooth, add = TRUE, segments = 51, ...) {
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    if (!is.na(r)) {
      if (abs(r) > 0) 
        theta <- sign(r)/sqrt(2)
      else theta = 1/sqrt(2)
      shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
                                                              theta, -theta, theta), ncol = 2, byrow = TRUE)
      ellipse <- unit.circle %*% shape
      ellipse[, 1] <- ellipse[, 1] * xs + x
      ellipse[, 2] <- ellipse[, 2] * ys + y
      if (show.points) 
        points(x, y, pch = 19, col = col.smooth, cex = 1.5)
      lines(ellipse, ...)
    }
  }
  "panel.ellipse" <- function(x, y, pch = par("pch"), col.smooth = "red", 
                              ...) {
    segments = 51
    usr <- par("usr")
    on.exit(par("usr"))
    par(usr = c(usr[1] - abs(0.05 * usr[1]), usr[2] + abs(0.05 * 
                                                            usr[2]), 0, 1.5))
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    r = cor(x, y, use = "pairwise", method = method)
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (smoother) {
      smoothScatter(x, y, add = TRUE, nrpoints = 0)
    }
    else {
      if (show.points) {
        points(x, y, pch = pch, ...)
      }
    }
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    if (!is.na(r)) {
      if (abs(r) > 0) 
        theta <- sign(r)/sqrt(2)
      else theta = 1/sqrt(2)
      shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
                                                              theta, -theta, theta), ncol = 2, byrow = TRUE)
      ellipse <- unit.circle %*% shape
      ellipse[, 1] <- ellipse[, 1] * xs + xm
      ellipse[, 2] <- ellipse[, 2] * ys + ym
      points(xm, ym, pch = 19, col = col.smooth, cex = 1.5)
      if (ellipses) 
        lines(ellipse, ...)
    }
  }
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (missing(cex.cor)) 
    cex.cor <- 1
  for (i in 1:ncol(x)) {
    if (is.character(x[[i]])) {
      x[[i]] <- as.numeric(as.factor(x[[i]]))
      colnames(x)[i] <- paste(colnames(x)[i], "*", sep = "")
    }
  }
  n.obs <- nrow(x)
  confid <- qt(1 - alpha/2, n.obs - 2)
  if (!lm) {
    if (cor) {
      pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, 
            lower.panel = panel.smoother, pch = pch, ...)
    }
    else {
      pairs(x, diag.panel = panel.hist.density, upper.panel = panel.smoother, 
            lower.panel = panel.smoother, pch = pch, ...)
    }
  }
  else {
    if (!cor) {
      pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, 
            lower.panel = panel.lm, pch = pch,  ...)
    }
    else {
      pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, 
            lower.panel = panel.lm, pch = pch,  ...)
    }
  }
}


pdf(file.path("C:/Projects/CXCR6/manuscript_ready/FigSupp_TRM_LUAD_panel.pdf"), width = 5, height = 5)
pairs.panels2(df.signatures[,c(1:3,6)],
              stars = T,
              method = "pearson", # correlation method
              hist.col = "#00AFBB",
              # title = paste(full_names[1], paste0('(',c_type,')')),
              # sub = paste(dim(df.signatures)[1], "patients"),
              main = paste(full_names[1], paste0('(',c_type,')')), 
              density = TRUE,  # show density plots
              ellipses = TRUE # show correlation ellipses
)
title("", sub = paste(nrow(df.signatures), "patients"))
dev.off()

pdf(file.path("C:/Projects/CXCR6/manuscript_ready/Fig9_TRM_CCL_panel_7.pdf"), width = 8, height = 8)
pairs.panels2(df.signatures2[,1:7],
              stars = T,
              method = "pearson", # correlation method
              hist.col = "#00AFBB", 
              main = paste(full_names[1], paste0('(',c_type,')')),
              density = TRUE,  # show density plots
              ellipses = TRUE # show correlation ellipses
)
title(main = "", sub = paste(nrow(df.signatures2), "patients"))
dev.off()