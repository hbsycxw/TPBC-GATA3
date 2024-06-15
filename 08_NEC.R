# 

rm(list = ls())

library(tidyverse)
library(GSVA)
library(readxl)
library(ggpubr)
library(ggsci)

library(HGNChelper)

GSE123845_group <- read_excel("./GEO/GSE123845/GSE123845_group.xlsx")
X41467_2020_19933_MOESM3_ESM <- read_excel("GEO/GSE123845/41467_2020_19933_MOESM3_ESM.xlsx")

GSE123845_group_tidy <- GSE123845_group %>% left_join(X41467_2020_19933_MOESM3_ESM, by = c("SAMPLE_ID" = "sample_id")) %>% filter(er_status_diagnosis == 1, pr_status_diagnosis == 1, her2_status_diagnosis == 1) %>% select(SAMPLE_ID, NAC_response)
save(GSE123845_group, file = "./GEO/GSE123845_group_tidy.rda")


GSE123845_exp_tpm <- read_csv("GEO/GSE12385/GSE123845_exp_tpm_matrix.csv") %>% aggregate(. ~ gene_id, . , mean)

GSE123845_exp_tpm_GATA3 <- GSE123845_exp_tpm %>% filter(gene_id == "GATA3") %>% column_to_rownames("gene_id") %>% t() %>% as.data.frame() %>% rownames_to_column("SAMPLE_ID")

GSE123845_group_GATA3 <- GSE123845_group_tidy %>% left_join(GSE123845_exp_tpm_GATA3,"SAMPLE_ID")


# 
df <- GSE123845_group_GATA3 %>% rename(tumor_group = NAC_response)
compaired <- list(c("pCR", "RD"))   
ggboxplot(df,
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

# ssGSEA
load("./ssGSEA/Cell_maker_ssGSEA.Rda")
ssGSEA_BRCA <- GSE123845_exp_tpm %>% column_to_rownames("gene_id")

ssGSEA_TPBC <- ssGSEA_BRCA[,GSE123845_group_tidy$SAMPLE_ID]

ssGSEA_data <- gsva(as.matrix(ssGSEA_TPBC),
                    Cell_maker, 
                    method = "ssgsea")


ssGSEA_data1 <- ssGSEA_data %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SAMPLE_ID") %>% inner_join(GSE123845_group_GATA3, "SAMPLE_ID") 
ssGSEA_data1$NAC_response <- factor(ssGSEA_data1$NAC_response,levels = c("pCR", "RD"))

ssGSEA_data2 <- ssGSEA_data1 %>%
  pivot_longer(cols = 2:29,
               names_to= "celltype",
               values_to = "Immune_value")

ggplot(ssGSEA_data2, 
       aes(x = celltype,
           y = Immune_value)) +
  geom_boxplot(aes(color = NAC_response),
               outlier.shape = NA) +
  theme_classic(base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1, 
                                   colour = "black"))+
  stat_compare_means(aes(group = NAC_response), 
                     label = "p.signif",
                     method = "wilcox.test",
                     label.x = 0.5,
                     label.y = 1.3,
                     bracket.size = 0.5) +
  scale_color_manual(values = c("#3B4992FF","#EE0000FF"))


df_GATA3_immunecell <- ssGSEA_data1 %>% select(SAMPLE_ID, `Activated CD4 T cell`, GATA3) %>% rename(Activated_CD4_T_cell = `Activated CD4 T cell`)
ggscatter(df_GATA3_immunecell, x = "Activated_CD4_T_cell", y = "GATA3", 
          color = "coral2",fill = "lightgray", size = 2.5,
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(label.x.npc = 0.6,label.y.npc = 0.9),
          cor.method = "spearman") + 
  theme(axis.line = element_line(size = 1.0)) +
  labs(x = "Activated CD4 T cell levels") +
  labs(y = "GATA3 mRNA expression")
