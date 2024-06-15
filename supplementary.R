# 01  ----
rm(list = ls())
library(tidyverse)
library(survminer)
library(survival)

load("GDCdata/TCGA_BRCA_sample_info.rda")
load("data/TCGA_BRCA_exp_anno.rda")
TCGA_BRCA_clinical_patient <- read_delim("data/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", delim = "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)

df_1 <- TCGA_BRCA_exp_anno %>% t() %>% as.data.frame() %>% dplyr::select(GATA3) %>% rownames_to_column("PATIENT_ID")
df_1$PATIENT_ID <- str_sub(df_1$PATIENT_ID,1,12)
df_2 <- TCGA_BRCA_clinical_patient %>% dplyr::select(PATIENT_ID,PAM50 = SUBTYPE,OS_STATUS,OS_MONTHS) %>% left_join(df_1,"PATIENT_ID") %>% drop_na(PAM50)

df_2$GATA3 <- log(df_2$GATA3 + 0.001)

df_2$OS_STATUS <- ifelse(df_2$OS_STATUS == "0:LIVING",0,1)

df_Basal <- df_2 %>% filter(PAM50 == "BRCA_Basal")
df_Her2 <- df_2 %>% filter(PAM50 == "BRCA_Her2")
df_LumA <- df_2 %>% filter(PAM50 == "BRCA_LumA")
df_LumB <- df_2 %>% filter(PAM50 == "BRCA_LumB")

library(survival)
library(survminer)

res.cut <- surv_cutpoint(df_Basal, 
                         time = "OS_MONTHS", 
                         event = "OS_STATUS", 
                         variables = "GATA3")
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$OS_MONTHS, res.cat$OS_STATUS)
group <- res.cat[,3]
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group )
p1 <- ggsurvplot(fit,
                 data = survival_dat,
                 pval = TRUE,
                 linetype = "solid",
                 surv.median.line = NULL,
                 title = NULL,
                 ylab = "Overall Survival (%)",
                 xlab = " Time (Months)",
                 legend.title = "Group",
                 legend.labs = c("GATA3-High","GATA3-Low"),
                 risk.table = F,
                 tables.height = 0.3,
                 palette = c("#ED0000FF","#0099B4FF"),
                 ggtheme = theme_classic(),
                 tables.theme = theme_survminer()
)

res.cut <- surv_cutpoint(df_Her2, 
                         time = "OS_MONTHS", 
                         event = "OS_STATUS", 
                         variables = "GATA3")
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$OS_MONTHS, res.cat$OS_STATUS)
group <- res.cat[,3]
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group )
p2 <- ggsurvplot(fit,
                 data = survival_dat,
                 pval = TRUE,
                 linetype = "solid",
                 surv.median.line = NULL,
                 title = NULL,
                 ylab = "Overall Survival (%)",
                 xlab = " Time (Months)",
                 legend.title = "Group",
                 legend.labs = c("GATA3-High","GATA3-Low"),
                 risk.table = F,
                 tables.height = 0.3,
                 palette = c("#ED0000FF","#0099B4FF"),
                 ggtheme = theme_classic(),
                 tables.theme = theme_survminer()
) 

res.cut <- surv_cutpoint(df_LumA, 
                         time = "OS_MONTHS", 
                         event = "OS_STATUS", 
                         variables = "GATA3")
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$OS_MONTHS, res.cat$OS_STATUS)
group <- res.cat[,3]
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group )
p3 <- ggsurvplot(fit,
                 data = survival_dat,
                 pval = TRUE,
                 linetype = "solid",
                 surv.median.line = NULL,
                 title = NULL,
                 ylab = "Overall Survival (%)",
                 xlab = " Time (Months)",
                 legend.title = "Group",
                 legend.labs = c("GATA3-High","GATA3-Low"),
                 risk.table = F,
                 tables.height = 0.3,
                 palette = c("#ED0000FF","#0099B4FF"),
                 ggtheme = theme_classic(),
                 tables.theme = theme_survminer()
) 

res.cut <- surv_cutpoint(df_LumB, 
                         time = "OS_MONTHS", 
                         event = "OS_STATUS", 
                         variables = "GATA3")
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$OS_MONTHS, res.cat$OS_STATUS)
group <- res.cat[,3]
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group )
p4 <- ggsurvplot(fit,
                 data = survival_dat,
                 pval = TRUE,
                 linetype = "solid",
                 surv.median.line = NULL,
                 title = NULL,
                 ylab = "Overall Survival (%)",
                 xlab = " Time (Months)",
                 legend.title = "Group",
                 legend.labs = c("GATA3-High","GATA3-Low"),
                 risk.table = F,
                 tables.height = 0.3,
                 palette = c("#ED0000FF","#0099B4FF"),
                 ggtheme = theme_classic(),
                 tables.theme = theme_survminer()
) 

splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2
splots[[3]] <- p3
splots[[4]] <- p4


res <- arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 2) #定义行数和列数

ggsave("FigS1.pdf", res)

# 02 
# oncopredict

library(oncoPredict)
library(tidyverse)


## 
trainingExprData <- readRDS(file='./oncopredict/CTRP2_Expr (TPM, not log transformed).rds')
dim(trainingExprData) #
trainingExprData[1:4,1:4]

#IMPORTANT note: here I do e^IC50 since the IC50s are actual ln values/log transformed already, and the calcPhenotype function Paul #has will do a power transformation (I assumed it would be better to not have both transformations)
trainingPtype <- readRDS(file="./oncopredict/CTRP2_Res.rds")
trainingPtype[1:4,1:4]

## 
load(file = "data/TCGA_BRCA_exp_anno.rda")
load(file = "./GDCdata/gene_info.rda")
load(file = "./data/TCGA_TPBC_data_tidy_group.rda")
gene_info <- as.data.frame(gene_info) %>% filter(gene_type == "protein_coding") 

colnames(TCGA_BRCA_exp_anno) <- str_sub(colnames(TCGA_BRCA_exp_anno),1,12)

# 
tpm_in <- TCGA_BRCA_exp_anno[gene_info$gene_name,TCGA_TPBC_data_tidy_group$PATIENT_ID] %>%
  as.matrix() 

testExprData = tpm_in
testExprData[1:4,1:4]

#  https://github.com/cran/oncoPredict/blob/master/R/CALCPHENOTYPE.R
batchCorrect<-"eb"
powerTransformPhenotype<-TRUE
removeLowVaryingGenes<-0.2
removeLowVaringGenesFrom<-"homogenizeData"
minNumSamples=10
selection<- 1
printOutput=TRUE
pcr=FALSE
report_pc=FALSEcc=FALSE
rsq=FALSE
percent=80


# calcPhenotype


calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=testExprData,
              batchCorrect=batchCorrect,
              powerTransformPhenotype=powerTransformPhenotype,
              removeLowVaryingGenes=removeLowVaryingGenes,
              removeLowVaringGenesFrom=removeLowVaringGenesFrom,
              minNumSamples=minNumSamples,
              selection=selection,
              printOutput=printOutput,
              pcr=pcr,
              report_pc=report_pc,
              percent=percent,
              rsq=rsq)


DrugPredictions <- read_csv("calcPhenotype_Output/DrugPredictions.csv")
DrugPredictions <- DrugPredictions %>% rename(PATIENT_ID = "...1")
DrugPredictions_group <- DrugPredictions %>% left_join(TCGA_TPBC_data_tidy_group,"PATIENT_ID") 

library(ggpubr)
library(ggsci)
my_comparisons_group <- list(c("low", "high"))
p5 <- ggboxplot(DrugPredictions_group, x="GATA3_group", y="tamoxifen", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))

p6 <- ggboxplot(DrugPredictions_group, x="GATA3_group", y="fulvestrant", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))

p7 <- ggboxplot(DrugPredictions_group, x="GATA3_group", y="gefitinib", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))
df  <-  DrugPredictions_group %>% select(GATA3_group, 2:546) %>% 
  pivot_longer(-1, names_to = "grp", values_to = "val")
df_wilcox <- df %>% group_by(grp) %>% 
  wilcox_test(val ~ GATA3_group, detailed = TRUE)

df_wilcox_p <- df_wilcox %>% filter(p < 0.05)

library(patchwork)
p5 + p6  + plot_layout(nrow=1) + plot_annotation(tag_levels = 'A')


