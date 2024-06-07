## 00 加载包 ----
library(tidyverse)
library(readr)
library(maftools)

library(ggpubr)
library(ggsci)

library(compareGroups) # 制作表1

library(limma)
library(pheatmap)
library(DESeq2)

library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment) #assay
library(mice) # 缺失值可视化

## 01 TCGA数据整理 ----
### 01-01 表达数据整理 ----
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts")

GDCdownload(query,
            method = "api",
            files.per.chunk = 10)

TCGA_BRCA_exp_RAW <- GDCprepare(
  query = query,
  save = TRUE, 
  save.filename = "./GDCdata/TCGA_BRCA_exp_RAW.rda"
)

TCGA_BRCA_exp_RAW <- data
TCGA_BRCA_sample_info <- as.data.frame(TCGA_BRCA_exp_RAW@colData)

TCGA_BRCA_tpm <- assay(TCGA_BRCA_exp_RAW, i = "tpm_unstrand") %>% as.data.frame()
TCGA_BRCA_counts <- assay(TCGA_BRCA_exp_RAW, i = "unstranded") %>% as.data.frame()

save(TCGA_BRCA_sample_info, file = "./GDCdata/TCGA_BRCA_sample_info.rda")
save(TCGA_BRCA_tpm, file = "./GDCdata/TCGA_BRCA_tpm.rda")
save(TCGA_BRCA_counts, file = "./GDCdata/TCGA_BRCA_counts.rda")

gene_info <- TCGA_BRCA_exp_RAW@rowRanges 

save(gene_info, file = "./GDCdata/gene_info.rda")

### 01-02 临床数据整理 ----
rm(list = ls())
library(TCGAbiolinks)
load("GDCdata/TCGA_BRCA_sample_info.rda")

query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
GDCdownload(query, method = "api")

TCGA_BRCA_clinical <- GDCprepare(query = query,
                                 save = TRUE,
                                 save.filename = "TCGA_BRCA_clinical.rda")
save(TCGA_BRCA_clinical, file = "GDCdata/TCGA_BRCA_clinical.rda")
TCGA_BRCA_clinical_patient_original <- TCGA_BRCA_clinical$clinical_patient_brca

TCGA_BRCA_clinical_patient <- read_delim("data/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", delim = "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)

TCGA_BRCA_clinical_sample <- read_delim("data/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt", delim = "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)

write.csv(TCGA_BRCA_clinical_patient_original, file ="./data/TCGA_BRCA_clinical_patient_original.csv")

# 按照指南规范编辑HER2状态
TCGA_BRCA_clinical_patient_IHC <- read_csv("./data/TCGA_BRCA_clinical_patient_tidy3.csv")

TCGA_BRCA_clinical <- TCGA_BRCA_clinical_patient_IHC %>%
  left_join(TCGA_BRCA_clinical_patient, by = "PATIENT_ID") %>% 
  left_join(TCGA_BRCA_clinical_sample, by = "PATIENT_ID") %>% 
  left_join(TCGA_BRCA_clinical_patient_IHC, by = "PATIENT_ID") %>% 
  select(
    PATIENT_ID,SAMPLE_ID,
    SUBTYPE,AGE,AJCC_PATHOLOGIC_TUMOR_STAGE,
    PATH_M_STAGE,PATH_N_STAGE,PATH_T_STAGE,
    OS_STATUS,OS_MONTHS,
    DSS_STATUS,DSS_MONTHS,
    PFS_STATUS,PFS_MONTHS,
    DFS_STATUS,DFS_MONTHS,
    PR_STATUS_BY_IHC = pr_status_by_ihc.x, 
    ER_STATUS_BY_IHC = er_status_by_ihc.x, 
    HER2_STATUS = her2_status.x,
    CANCER_TYPE_DETAILED,
    ANEUPLOIDY_SCORE,
    MSI_SCORE_MANTIS,
    MSI_SENSOR_SCORE,TMB_NONSYNONYMOUS) 


TCGA_BRCA_clinical$OS_STATUS <- ifelse(TCGA_BRCA_clinical$OS_STATUS == "0:LIVING",0,1)
TCGA_BRCA_clinical$PFS_STATUS <- ifelse(TCGA_BRCA_clinical$PFS_STATUS == "0:CENSORED",0,1) 
TCGA_BRCA_clinical$DSS_STATUS <- ifelse(TCGA_BRCA_clinical$DSS_STATUS == "0:ALIVE OR DEAD TUMOR FREE",0,1) 
TCGA_BRCA_clinical$DFS_STATUS <- ifelse(TCGA_BRCA_clinical$DFS_STATUS == "0:DiseaseFree",0,1) 

# library(VIM) #探索缺失值
# aggr(TCGA_BRCA_clinical)

# TCGA_BRCA_clinical$Tumor_Sample_Barcode <- TCGA_BRCA_clinical$SAMPLE_ID
# 
# # 删除重复样本
# patient_duplicated <- TCGA_BRCA_clinical$PATIENT_ID %>% table() %>% as.data.frame() %>% filter(Freq >=2 ) %>% dplyr::rename(PATIENT_ID = ".")
# 
# TCGA_BRCA_clinical <- TCGA_BRCA_clinical %>% filter(!PATIENT_ID %in% patient_duplicated$PATIENT_ID)
# 
TCGA_BRCA_clinical_tidy <- TCGA_BRCA_clinical #余609

save(TCGA_BRCA_clinical_tidy,file = "./data/TCGA_BRCA_clinical_tidy.rda")
### 01-03 突变数据整理 ----
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation")
GDCdownload(query, method = "api")

TCGA_BRCA_snp <- GDCprepare(query = query,
                            save = TRUE,
                            save.filename = "./GDCdata/TCGA_BRCA_snp.rda")
TCGA_BRCA_snp <- data
save(TCGA_BRCA_snp, file = "./GDCdata/TCGA_BRCA_snp.rda")
# 过滤
Variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
             "In_Frame_Ins","Missense_Mutation","Nonsense_Mutation",
             "Nonstop_Mutation", "Splice_Region","Splice_Site",
             "Translation_Start_Site")

TCGA_BRCA_mutation_filter <- TCGA_BRCA_snp %>% mutate(VAF = t_alt_count/t_depth) %>% filter(VAF >= 0.05) %>%  filter(Variant_Classification %in% Variants) 

TCGA_BRCA_mutation_filter$Tumor_Sample_Barcode <- str_sub(TCGA_BRCA_mutation_filter$Tumor_Sample_Barcode, 1, 15) 

save(TCGA_BRCA_mutation_filter, file = "./GDCdata/TCGA_BRCA_mutation_filter.rda")

# 01-04 甲基化数据 ----
illuminaMethyl450_hg38_GDC <- read_delim("data/xena/illuminaMethyl450_hg38_GDC", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE) %>% dplyr::rename(id = `#id` )

save(illuminaMethyl450_hg38_GDC, file = "data/illuminaMethyl450_hg38_GDC.rda")

TCGA_BRCA_methylation450_tsv <- read_delim("data/xena/TCGA-BRCA.methylation450.tsv.gz", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
TCGA_BRCA_methylation450_tsv <- TCGA_BRCA_methylation450_tsv %>% column_to_rownames("Composite Element REF") 
colnames(TCGA_BRCA_methylation450_tsv) <- str_sub(colnames(TCGA_BRCA_methylation450_tsv),1,15)
keep <- colnames(TCGA_BRCA_methylation450_tsv) %in% TCGA_TPBC_data_tidy_group$SAMPLE_ID
TCGA_TPBC_meth <- TCGA_BRCA_methylation450_tsv[,keep]
save(TCGA_TPBC_meth, file = "data/TCGA_TPBC_meth.rda")

# 02 TCGA数据探索 ----
rm(list = ls())

load(file = "./GDCdata/TCGA_BRCA_tpm.rda")
load(file = "./GDCdata/gene_info.rda")
load(file = "./data/TCGA_BRCA_clinical_tidy.rda")

### GATA3 ----
gene_info <- as.data.frame(gene_info) %>% select(gene_name,gene_id)

TCGA_BRCA_exp_anno <- TCGA_BRCA_tpm %>% rownames_to_column("gene_id") %>% 
  left_join(gene_info,"gene_id") %>% 
  select(-gene_id) %>% aggregate(. ~ gene_name, . , mean) %>% 
  column_to_rownames("gene_name")

save(TCGA_BRCA_exp_anno, file = "data/TCGA_BRCA_exp_anno.rda")

TCGA_BRCA_GATA3_exp_tpm <- TCGA_BRCA_tpm %>% rownames_to_column("gene_id") %>% left_join(gene_info, by ="gene_id") %>% select(-gene_id) %>% filter(gene_name == "GATA3") %>% column_to_rownames("gene_name") %>% t() %>% as.data.frame() 

TCGA_BRCA_GATA3_exp_tpm$GATA3 <- as.numeric(TCGA_BRCA_GATA3_exp_tpm$GATA3)
TCGA_BRCA_GATA3_exp_tpm <- log(TCGA_BRCA_GATA3_exp_tpm + 0.001) %>% rownames_to_column("Tumor_Sample_Barcode") 

TCGA_BRCA_GATA3_exp_tpm$Tumor_Sample_Barcode <- str_sub(TCGA_BRCA_GATA3_exp_tpm$Tumor_Sample_Barcode, 1, 15)

TCGA_BRCA_clinical_GATA3 <- TCGA_BRCA_clinical_tidy %>% left_join(TCGA_BRCA_GATA3_exp_tpm, by = c("SAMPLE_ID"="Tumor_Sample_Barcode")) 

TCGA_TPBC_data <- TCGA_BRCA_clinical_GATA3 %>% filter(ER_STATUS_BY_IHC == "Positive") %>% filter(PR_STATUS_BY_IHC == "Positive") %>% filter(HER2_STATUS == "Positive") 
#drop_na一定要放在最后进行

load("GDCdata/TCGA_BRCA_sample_info.rda")
TCGA_normal_group <- TCGA_BRCA_sample_info %>% filter(shortLetterCode == "NT") %>% 
  mutate(Tumor_Sample_Barcode = str_sub(sample,1,15))
TCGA_BRCA_normal_data <- TCGA_BRCA_GATA3_exp_tpm %>% filter(Tumor_Sample_Barcode %in% TCGA_normal_group$Tumor_Sample_Barcode) %>% mutate(tumor_group = "NT")

save(TCGA_BRCA_normal_data, file = "./data/BRCA_normal.rda")
save(TCGA_TPBC_data, file = "./data/BRCA_TPBC_data.Rdata")

# 03 METABRIC数据整理 ----
METABRIC_exp_raw  <- read_delim("data/brca_metabric/data_mrna_agilent_microarray.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

METABRIC_exp <- METABRIC_exp_raw %>% select(-Entrez_Gene_Id) %>% aggregate(. ~ Hugo_Symbol, . , mean) %>%  column_to_rownames(var = "Hugo_Symbol")

save(METABRIC_exp, file = "./data/METABRIC_exp.Rdata")

## 临床数据
METABRIC_patient <- read_delim("data/brca_metabric/data_clinical_patient.txt", 
                               delim = "\t", escape_double = FALSE, 
                               comment = "#", trim_ws = TRUE) 


## 样本数据
METABRIC_sample <- read_delim("data/brca_metabric/data_clinical_sample.txt", 
                              delim = "\t", escape_double = FALSE, 
                              comment = "#", trim_ws = TRUE) 

## 合并临床与样本数据
METABRIC_clinical <- METABRIC_patient %>% left_join(METABRIC_sample, by = "PATIENT_ID")

METABRIC_clinical$OS_STATUS <- ifelse(METABRIC_clinical$OS_STATUS == "0:LIVING",0,1)
METABRIC_clinical$RFS_STATUS <- ifelse(METABRIC_clinical$RFS_STATUS == "0:Not Recurred",0,1)

METABRIC_clinical <- METABRIC_clinical %>% dplyr::rename(Tumor_Sample_Barcode = PATIENT_ID) 

METABRIC_exp_GATA3_TPBC <- METABRIC_exp["GATA3",] %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Tumor_Sample_Barcode") %>% left_join(METABRIC_clinical,by = "Tumor_Sample_Barcode") %>% filter(ER_STATUS == "Positive") %>% filter(PR_STATUS == "Positive") %>% filter(HER2_STATUS == "Positive") 

save(METABRIC_exp_GATA3_TPBC, file = "./data/METABRIC_exp_GATA3_TPBC.rda")

