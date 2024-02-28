library(dittoSeq)
library(Seurat)
library(ggplot2)
library(dplyr)
covid_MR1 = readRDS("C:/Projects/PBMC_covid_sc_/covid_meta_MR1.rds")


head(covid_MR1)

g1 = ggplot(covid_MR1, aes(y =  MR1, x = Condition, fill = Condition))+
  geom_violin(scale = "count", trim = F,  show.legend = F)+
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
  theme_classic2()+
  scale_fill_brewer(type = "qual", palette = "Set1")+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

table(covid_MR1$MR1_raw_count, covid_MR1$patient_id)
table(covid_MR1$initial_clustering)


ggplot(covid_MR1[covid_MR1$Condition=="Healthy", ], 
       aes(y =  MR1, x = initial_clustering, fill = initial_clustering))+
  geom_violin(scale = "count", trim = F)+
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5)+
  theme_classic2()+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
head(covid_MR1)

covid_MR1_clean = covid_MR1[covid_MR1$initial_clustering %in% c("B_cell",   "CD4" ,  "CD8",  "CD14" , "CD16"  , 
                                                                "DCs", "NK_16hi", "NK_56hi" ,"Treg"  , "gdT", "pDC"), ]
covid_MR1_clean$cell_type = covid_MR1_clean$initial_clustering
covid_MR1_clean$cell_type = gsub("NK_16hi|NK_56hi", "NK", covid_MR1_clean$cell_type)
covid_MR1_clean$cell_type = gsub("CD14|CD16", "Monocytes", covid_MR1_clean$cell_type)
covid_MR1_clean$cell_type = gsub("CD4|Treg", "CD4", covid_MR1_clean$cell_type)




# ggplot(covid_MR1_clean, aes(fill=cell_type, y=MR1, x=Condition)) + 
#   geom_bar(position="fill", stat="identity")
# 
# ggplot(covid_MR1_clean[covid_MR1_clean$Condition=="Healthy", ], aes(fill=cell_type, y=MR1_raw_count, x=cell_type)) + 
#   geom_bar(position="dodge", stat="identity")
# 
# table(covid_MR1[covid_MR1$Condition=="Healthy",]$initial_clustering)
# covid_heal = covid_MR1[covid_MR1$Condition=="Healthy",]
# table(covid_heal$MR1_raw_count, covid_heal$initial_clustering)
# 
# table(covid_MR1_clean[covid_MR1_clean$Condition=="Healthy",]$cell_type)

covid_means <- covid_MR1_clean[covid_MR1_clean$Condition=="Healthy",] %>% 
  group_by( cell_type) %>% 
  summarise( count = n(), mean = mean(MR1), prop = sum(MR1_raw_count>0)/count)
head(covid_means)
# making the plot
g2 = 
  ggplot(covid_MR1_clean[covid_MR1_clean$Condition=="Healthy",] , 
         aes(y =  MR1, x = cell_type, fill = cell_type))+
  geom_violin(scale = "count", trim = F,  show.legend = F)+
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
  theme_classic2()+
  xlab("Cells from Healthy individuals")+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggplot(covid_means, aes(x = cell_type, y = mean, fill = cell_type, label = count)) +
# geom_bar(stat = "identity", position = "dodge", show.legend = F)+
# # geom_label(show.legend = F)+
# ylab("Mean MR1 expression")+
# theme_classic2()+
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
g3 = ggplot(covid_means, aes(x = cell_type, y = prop, fill = cell_type, label = count)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = F)+
  geom_label(show.legend = F)+
  ylab("Proportion of MR1>0 expressing cells")+
  xlab("Cells from Healthy individuals")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dim(covid_MR1_clean[covid_MR1_clean$Condition %in% c("Healthy", "Mild", "Moderate", "Severe") & 
                      covid_MR1_clean$initial_clustering=="DCs",])
table(covid_MR1$Condition)

g4 = ggplot(covid_MR1[covid_MR1$Condition %in% c("Healthy", "Mild", "Moderate", "Severe") & 
                        # substr(covid_MR1$patient_id, start = 1, stop = 2)=="MH" &
                        covid_MR1$initial_clustering=="DCs" & 
                        covid_MR1$MR1_raw_count>0,] , 
            aes(y =  MR1, x = Condition, fill = Condition))+
  geom_violin(scale = "count", trim = F,  show.legend = F)+
  geom_jitter(size = 0.5,  width = 0.1,show.legend = F)+
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5,  show.legend = F)+
  stat_compare_means(ref.group = "Healthy", label = "p.format")+
  theme_classic2()+
  xlab("DCs, Status")+
  scale_y_log10()+ylab("MR1 (log scale)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("C:/Projects/PBMC_covid_sc_/MR1.pdf", width = 10, height = 10)
gridExtra::grid.arrange(g1,g2,g3,g4, ncol = 2)
dev.off()
ggsave("C:/Projects/PBMC_covid_sc_/MR1.pdf", plot = gg, width = 12, height = 9)
