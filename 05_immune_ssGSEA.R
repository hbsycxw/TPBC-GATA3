# immune ----
rm(list = ls())
library(GSVA)
library(tidyverse)
load("./ssGSEA/Cell_maker_ssGSEA.Rda")
load("data/BRCA_TPBC_data.rda")
load("./data/TCGA_TPBC_data_tidy_group.rda")
load("./GDCdata/TCGA_BRCA_tpm.rda")
load("./GDCdata/gene_info.rda")
load("./GDCdata/TCGA_BRCA_tpm_anno.rda")

ssGSEA_BRCA <- TCGA_BRCA_tpm_anno
colnames(ssGSEA_BRCA) <- str_sub(colnames(ssGSEA_BRCA),1,15)

ssGSEA_TPBC <- ssGSEA_BRCA[,TCGA_TPBC_data$SAMPLE_ID]

ssGSEA_data <- gsva(as.matrix(ssGSEA_TPBC),
                    Cell_maker, 
                    method = "ssgsea")
save(ssGSEA_data,file = "./ssGSEA/TCGA_TPBC_ssGSEA_data.Rda")


data_1 <- TCGA_TPBC_data_tidy_group %>% select(PATIENT_ID,GATA3_group)
ssGSEA_data1 <- ssGSEA_data %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID") %>% mutate(PATIENT_ID = str_sub(PATIENT_ID,1,12)) %>% inner_join(data_1, "PATIENT_ID") 
ssGSEA_data1$GATA3_group <- factor(ssGSEA_data1$GATA3_group,levels = c("low", "high"))


ssGSEA_data2 <- ssGSEA_data1 %>%
  dplyr::select(1:29, 30) %>%
  pivot_longer(cols = 2:29,
               names_to= "celltype",
               values_to = "Immune_value")

library(ggplot2)
library(ggpubr)


p1 <- ggplot(ssGSEA_data2, 
       aes(x = celltype,
           y = Immune_value)) +
  geom_boxplot(aes(color = GATA3_group),
               outlier.shape = NA) +
  theme_classic(base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1, 
                                   colour = "black"))+
  stat_compare_means(aes(group = GATA3_group), 
                     label = "p.signif",
                     method = "wilcox.test",
                     label.x = 0.5,
                     label.y = 1.3,
                     bracket.size = 0.5) +
  scale_color_manual(values = c("#3B4992FF","#EE0000FF"))

# Estimate
rforge <-"http://r-forge.r-project.org"
install.packages("estimate",repos=rforge,dependencies=TRUE)
library(estimate)
load(file = "./GDCdata/TCGA_BRCA_tpm_anno.rda")
load(file = "./data/TCGA_TPBC_data_tidy_group.rda")
write.table(TCGA_BRCA_tpm_anno, file ="ssGSEA/TCGA_BRCA_tpm_anno.txt", sep ="\t", row.names =T, col.names =T, quote =F)

estimate_data <- "./ssGSEA/TCGA_BRCA_tpm_anno.txt"
commonGenes <- "./ssGSEA/commonGenes.gct"


filterCommonGenes(input.f = estimate_data, 
                  output.f = commonGenes, 
                  id = "GeneSymbol")
estimateScore(input.ds = "./ssGSEA/commonGenes.gct",
              output.ds="./ssGSEA/estimateScore.gct", 
              platform="agilent")


scores <- read.table("./ssGSEA/estimateScore.gct", skip = 2, header = T)
groups <- TCGA_TPBC_data_tidy_group %>% select(PATIENT_ID,GATA3_group)
scores1 <- scores %>%
  column_to_rownames("NAME") %>%
  dplyr::select(-Description) %>%
  t() %>%
  as.data.frame() %>% rownames_to_column(var = "PATIENT_ID") %>% 
  mutate(PATIENT_ID = str_sub(PATIENT_ID,1,12)) %>% 
  mutate(PATIENT_ID = str_replace_all(PATIENT_ID,"\\.","-")) %>% 
  inner_join(groups, by = "PATIENT_ID")

save(scores1, file = "ssGSEA/ESTIMATE.rda")


# StromalScore
library(ggpubr)
library(ggsci)
my_comparisons_location <- list(c("low", "high"))
p3 <- ggboxplot(scores1, x="GATA3_group", y="StromalScore", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_location,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))
# ImmuneScore
p4 <- ggboxplot(scores1, x="GATA3_group", y="ImmuneScore", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_location,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))

# ESTIMATEScore
p5 <- ggboxplot(scores1, x="GATA3_group", y="ESTIMATEScore", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_location,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))

## immune related genes 
## form https://www.immport.org/home
library(GSVA)
# 
immune_genes <- read_delim("ssGSEA/GeneList_immport.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

type <- split(immune_genes,immune_genes$Category)
immune_genes_list <- lapply(type, function(x){
  dd = x$Symbol
  unique(dd)
})

groups <- TCGA_TPBC_data_tidy_group %>% select(PATIENT_ID,GATA3_group)
colnames(TCGA_BRCA_tpm_anno) <- str_sub(colnames(TCGA_BRCA_tpm_anno),1,12)

exp_immunegenes <- log(TCGA_BRCA_tpm_anno[,TCGA_TPBC_data_tidy_group$PATIENT_ID]+0.001) %>% as.matrix()

immune_gsva <- gsva(exp_immunegenes, immune_genes_list)
immune_gsva1 <- immune_gsva %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID") %>% inner_join(groups, "PATIENT_ID") 
immune_gsva1$GATA3_group <- factor(immune_gsva1$GATA3_group,levels = c("low", "high"))

immune_gsva2 <- immune_gsva1 %>%
  dplyr::select(1:18, 19) %>%
  pivot_longer(cols = 2:18,
               names_to= "immune_genes",
               values_to = "Immune_values")


p6 <- ggplot(immune_gsva2, 
       aes(x = immune_genes,
           y = Immune_values)) +
  geom_boxplot(aes(color = GATA3_group),
               outlier.shape = NA) +
  geom_jitter(data=immune_gsva2, 
              mapping=aes(x= immune_genes,
                          y= Immune_values,
                          fill = GATA3_group,
                          color= GATA3_group),
              size = 0.8,
              alpha= 0.5,
              shape = 16,
              stroke = 0.15, show.legend = FALSE, 
              position = position_jitterdodge(jitter.height= 0.05, 
                                              jitter.width = 0.5, 
                                              dodge.width = 0.8))+
  theme_classic(base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1, 
                                   colour = "black"))+
  stat_compare_means(aes(group = GATA3_group), 
                     label = "p.signif",
                     method = "wilcox.test",
                     label.x = 0.5,
                     label.y = 1.3,
                     bracket.size = 0.5) +
  scale_color_manual(values = c("#3B4992FF","#EE0000FF"))


# GSVA DEG limma ---
## DEG
library(ggplot2)
library(ggprism)

group_list <- groups %>% arrange(desc(GATA3_group))
group_list$GATA3_group <-factor(group_list$GATA3_group,levels = c("low","high"))
table(group_list$GATA3_group)
es.max <- immune_gsva[,group_list$PATIENT_ID]
library(limma)
design <- model.matrix(~ 0 + factor(group_list$GATA3_group))
colnames(design)<-levels(factor(group_list$GATA3_group))
rownames(design)<-colnames(es.max)
contrast.matrix<-makeContrasts("high-low",
                               levels = design)
contrast.matrix

fit <-lmFit(es.max,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

res <- decideTests(fit2,p.value=0.05)
summary (res)

tmp <- topTable(fit2,coef = 1,number = 200)
DEG <- na.omit(tmp)
head(DEG[1:4,])

DEG_sig <- DEG %>% rownames_to_column("geneset") %>% 
  mutate(gene_label=ifelse(geneset %in% c("TCRsignalingPathway"), geneset , "")) %>% mutate(sig=ifelse((logFC >= 0.2 & adj.P.Val < 0.05), "up", ifelse((logFC <= -0.2 & adj.P.Val < 0.05), "down", "notsig")))
colours <- setNames(c( "#009392", "darkgrey", "#CF597E"), c("down", "notsig", "up"))

G_M <- tmp %>%
  rownames_to_column("ID") %>%
  dplyr::select(1,4)


G_M$group <- if_else(G_M$t < -4.28, "1",
                     if_else(G_M$t > -2, "3",
                             "2"))
table(G_M$group)


sortG_M <- G_M[order(G_M$t),]
sortG_M$ID <- factor(sortG_M$ID, 
                     levels = sortG_M$ID)
head(sortG_M)

p7 <- ggplot(sortG_M, 
       aes(ID, t, fill = group)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  scale_fill_manual(values = c("steelblue", "gray70", "darkred"),
                    guide = FALSE) + 
  geom_hline(yintercept = c(-4,-2), 
             color="white",
             linetype = 2, 
             size = 0.5)+
  xlab('') + 
  ylab('t value of GSVA score, high versus low GATA3 group') + 
  theme_prism(border = T) +
  theme(
    axis.ticks.y = element_blank()
  )

# PD-L1 PD-1 CTLA-4 LAG3
## PD-L1 ~ CD274 ; PD-1 ~ PDCD1; CTLA-4 ~ CTLA4
immunotherapy_genes <- c("CD274","PDCD1","CTLA4","LAG3")
groups <- TCGA_TPBC_data_tidy_group %>% select(SAMPLE_ID,GATA3_group)
load(file = "./GDCdata/TCGA_BRCA_tpm_anno.rda")
exp_immunotherapy <- log(TCGA_BRCA_tpm_anno[immunotherapy_genes,] + 0.001) %>% t() %>% as.data.frame() %>%
  rownames_to_column("SAMPLE_ID") %>% mutate(SAMPLE_ID = str_sub(SAMPLE_ID,1,15)) %>% inner_join(groups, "SAMPLE_ID")

exp_immunotherapy$GATA3_group <- factor(exp_immunotherapy$GATA3_group,levels = c("low", "high"))

df <- exp_immunotherapy %>% 
  pivot_longer(cols = 2:5,
               names_to= "gene_symbol",
               values_to = "expression_value")

p8 <- ggplot(df, 
       aes(x = gene_symbol,
           y = expression_value)) +
  geom_boxplot(aes(color = GATA3_group),
               outlier.shape = NA) +
  theme_classic(base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1, 
                                   colour = "black"),
        legend.position = "top")+
  geom_jitter(data=df, 
              mapping=aes(x= gene_symbol,
                          y= expression_value,
                          fill = GATA3_group,
                          color= GATA3_group),
              size = 1.2,
              alpha= 0.5,
              shape = 16,
              stroke = 0.15, show.legend = FALSE, 
              position = position_jitterdodge(jitter.height= 0.05, 
                                              jitter.width = 0.5, 
                                              dodge.width = 0.8))+
  stat_compare_means(aes(group = GATA3_group), 
                     label = "p.signif",
                     method = "wilcox.test",
                     label.x = 0.5,
                     label.y = 5,
                     bracket.size = 0.5) +
  scale_color_manual(values = c("#3B4992FF","#EE0000FF")) +
  labs(x = NULL) +
  labs(y = "Gene Expression value")


# IPS TCIA
TCIA_ClinicalData_BRCA <- read_delim("TCIA/TCIA-ClinicalData_BRCA.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

groups <- TCGA_TPBC_data_tidy_group %>% select(PATIENT_ID,GATA3_group)

TCIA_TPBC_GATA3 <- TCIA_ClinicalData_BRCA %>% filter(barcode %in% groups$PATIENT_ID) %>% left_join(groups, by = c("barcode" = "PATIENT_ID")) %>% mutate(ips_ctla4_pos = ips_ctla4_pos_pd1_neg + ips_ctla4_pos_pd1_pos, ips_pd1_pos = ips_ctla4_neg_pd1_pos + ips_ctla4_pos_pd1_pos)

library(ggpubr)
library(ggsci)
my_comparisons_TCIA <- list(c("low", "high"))
p9 <- ggboxplot(TCIA_TPBC_GATA3, x = "GATA3_group" , y = "ips_ctla4_pos",
          color = "GATA3_group",
          add = "jitter", size = 1) +
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_TCIA,
                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS"))) +
  xlab("GATA3 Group") +
  ylab('IPS CTLA4 POS') 

p10 <- ggboxplot(TCIA_TPBC_GATA3, x = "GATA3_group" , y = "ips_pd1_pos",
                color = "GATA3_group",
                add = "jitter", size = 1) +
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_TCIA,
                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS"))) +
  xlab("GATA3 Group") +
  ylab('IPS PD-1 POS')

save(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,file = 'Figure/Fig5.RData')
