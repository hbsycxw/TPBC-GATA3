# 02 GATA3 in different Goups
rm(list = ls())
load("./data/BRCA_normal.rda")
load("./data/BRCA_TPBC_data.rda")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
df_tumor <- TCGA_TPBC_data %>% select(SAMPLE_ID,GATA3) %>% mutate(tumor_group = "TPBC")

df_normal <- TCGA_BRCA_normal_data %>% rename(SAMPLE_ID = Tumor_Sample_Barcode)

df_Basal <- df_Basal %>% select(SAMPLE_ID = PATIENT_ID,GATA3,tumor_group = PAM50)

df <- rbind(df_tumor,df_normal)
df$tumor_group <- factor(df$tumor_group, levels = c("NT", "TPBC"))

compaired <- list(c("NT", "TPBC"))   


p1 <- ggboxplot(df,
          x = "tumor_group", y = "GATA3",
          color = "tumor_group",
          add = "jitter", 
          add.params = list(color = "tumor_group",alpha=0.5,shape = 16)
) +   
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",  
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns"))) +
  scale_color_aaas() + 
  theme(legend.position = "none",
        axis.line = element_line(size = 1.0)) +
  labs(x = "Tumour Group") +
  labs(y = "GATA3 mRNA expression") 

## 

data_rppa_zscores <- read_delim("data/brca_tcga_pan_can_atlas_2018/data_rppa_zscores.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)

TCGA_RPPA_GATA3 <- data_rppa_zscores %>% filter( Composite.Element.REF %in% c("GATA3|GATA3") ) %>% column_to_rownames("Composite.Element.REF") %>% t() %>% as.data.frame() %>% rownames_to_column("SAMPLE_ID") %>% rename("GATA3_RPPA" = "GATA3|GATA3")


TCGA_RPPA_GATA3_group <- TCGA_RPPA_GATA3 %>% inner_join(df,"SAMPLE_ID")

p2 <- ggscatter(TCGA_RPPA_GATA3_group, x = "GATA3", y = "GATA3_RPPA", 
          color = "steelblue",fill = "lightgray", size = 2.5,
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = 0.9),
          cor.method = "spearman") + 
  theme(axis.line = element_line(size = 1.0)) +
  labs(x = "GATA3 mRNA expression") +
  labs(y = "GATA3 RPPA levels")

# 
load(file = "data/TCGA_BRCA_exp_anno.rda")

TPBC_ER_HER2_exp <- TCGA_BRCA_exp_anno %>% rownames_to_column("gene_name") %>% filter(gene_name %in% c("ESR1","ERBB2")) %>% column_to_rownames("gene_name") %>% t() %>% as.data.frame() %>% rownames_to_column("SAMPLE_ID") %>% mutate(SAMPLE_ID = str_sub(SAMPLE_ID,1,15)) %>% filter(SAMPLE_ID %in% TCGA_TPBC_data$SAMPLE_ID) %>% column_to_rownames("SAMPLE_ID")

TPBC_ER_HER2_exp_log <- log(TPBC_ER_HER2_exp + 0.001)

p3 <- ggscatter(TPBC_ER_HER2_exp_log, x = "ESR1", y = "ERBB2", 
                color = "coral2",fill = "lightgray", size = 2.5,
                add = "reg.line", conf.int = TRUE, 
                add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
                cor.coef = T,
                cor.coeff.args = list(label.x.npc = "left", label.y.npc = 0.9),
                cor.method = "spearman") + 
  theme(axis.line = element_line(size = 1.0)) +
  labs(x = "ESR1 mRNA expression") +
  labs(y = "ERBB2 mRNA expression")

# 

p1 + p2 + p3 + plot_layout(nrow=1) + plot_annotation(tag_levels = 'A')

save(p1,p2,p3, file = "Fig2.rda" )