rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)

load(file = "./data/TCGA_TPBC_data_tidy_group.rda")
load(file = "./GDCdata/TCGA_BRCA_counts.rda")
load(file = "./GDCdata/gene_info.rda")

colnames(TCGA_BRCA_counts) <- str_sub(colnames(TCGA_BRCA_counts), 1,15)
gene_info <- gene_info %>% as.data.frame() %>% dplyr::select(gene_name) %>% rownames_to_column("gene_symbol")

exp <- TCGA_BRCA_counts[,!duplicated(colnames(TCGA_BRCA_counts))] %>% t() %>% as.data.frame() %>% rownames_to_column("SAMPLE_ID") %>% filter(SAMPLE_ID %in% TCGA_TPBC_data_tidy_group$SAMPLE_ID) %>% column_to_rownames("SAMPLE_ID") %>% t() %>% as.data.frame()

TCGA_TPBC_group_list <- TCGA_TPBC_data_tidy_group %>% arrange(GATA3)

exp <- exp[,TCGA_TPBC_group_list$SAMPLE_ID]

dim(exp)
exp = exp[apply(exp, 1, function(x) sum(x > 1) > 9), ]
dim(exp)
exp[1:4,1:4]

keep_exp <- rowSums(exp>0) >= floor(0.75*ncol(exp))
table(keep_exp)
exp <- exp[keep_exp,]
dim(exp)
exp[1:4,1:4]

group_list = factor(TCGA_TPBC_group_list$GATA3_group,levels = c("low","high"))

library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=group_list)

dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = colData,
  design = ~ condition)
dds <- DESeq(dds)
save(dds,file = "./DEG/TCGA_TPBC_DESeq2_dds.rda")


res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] 
DEG <- as.data.frame(resOrdered)
head(DEG)


logFC_cutoff <- 1
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
DEG <- DEG %>% rownames_to_column("gene_symbol") %>% left_join(gene_info,"gene_symbol")



genelist <- bitr(DEG$gene_name, 
                 fromType="SYMBOL",
                 toType="ENTREZID", 
                 OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_name"="SYMBOL"))


geneList = DEG[,3]

names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)

geneList = sort(geneList, decreasing = TRUE)

gsego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               ont          = "CC",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05, #
               verbose      = FALSE)


go_result <- setReadable(gsego, OrgDb = org.Hs.eg.db)

go_result_df <- as.data.frame(go_result)



gsekk <- gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 1000, 
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
gsekk_result_df <- as.data.frame(gsekk)


kegmt_h <- read.gmt("./GSEA/h.all.v2023.1.Hs.entrez.gmt")
KEGG_h <-GSEA(geneList,TERM2GENE = kegmt_h)
KEGG_h_result_df <- as.data.frame(KEGG_h)

kegmt_c2 <- read.gmt("./GSEA/c2.all.v2023.1.Hs.entrez.gmt")
KEGG_c2 <-GSEA(geneList,TERM2GENE = kegmt_c2)
KEGG_c2_result_df <- as.data.frame(KEGG_c2)

kegmt_c5 <- read.gmt("./GSEA/c5.all.v2023.1.Hs.entrez.gmt")
KEGG_c5 <-GSEA(geneList,TERM2GENE = kegmt_c5)
KEGG_c5_result_df <- as.data.frame(KEGG_c5)

kegmt_h_c5 <- read.gmt("./GSEA/h_c5.all.v2023.1.Hs.entrez.gmt")
KEGG_h_c5 <-GSEA(geneList,TERM2GENE = kegmt_h_c5)
KEGG_h_c5_result_df <- as.data.frame(KEGG_h_c5)


gseaplot2(gsekk,geneSetID = c("hsa04062","hsa04630","hsa04014"))

gseaplot2(gsekk,geneSetID = 1:5)
gseaplot2(gsekk,geneSetID = 6:11)
gseaplot2(KEGG_h,geneSetID = c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_COMPLEMENT","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_IL2_STAT5_SIGNALING"))
gseaplot2(KEGG_c2,geneSetID = c("REACTOME_DISEASES_OF_DNA_REPAIR","KEGG_MISMATCH_REPAIR"))
gseaplot2(KEGG_c5,geneSetID = c("GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION","GOBP_ACTIVATION_OF_IMMUNE_RESPONSE"))

gseaplot2(KEGG_h_c5,geneSetID = c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_COMPLEMENT","GOBP_ACTIVATION_OF_IMMUNE_RESPONSE"))


rm(list = ls())
library(ggpubr)
library(ggsci)
library(tidyverse)
# MSI ----
load(file = "./data/TCGA_TPBC_data_tidy_group.rda")
load(file = "data/METABRIC_TPBC_GATA3_group.rda")
TCGA_TPBC_data_tidy_group$GATA3_group <- factor(TCGA_TPBC_data_tidy_group$GATA3_group,levels = c("low","high"))
my_comparisons_location <- list(c("low", "high"))

p1 <- ggboxplot(TCGA_TPBC_data_tidy_group, x="GATA3_group", y="MSI_SENSOR_SCORE", 
                color="GATA3_group",
                add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_location,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1), method = "wilcox.test",
                                        symbols = c("***", "**", "*", "NS"))) +
  labs(x = "GATA3 Group") +
  labs(y = "MSI SENSOR SCORE") 

# TMB ----

ggboxplot(TCGA_TPBC_data_tidy_group, x="GATA3_group", y="TMB_NONSYNONYMOUS", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_location,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1), method = "wilcox.test",
                                        symbols = c("***", "**", "*", "NS")))
METABRIC_TPBC_GATA3_group$GATA3_group <- factor(METABRIC_TPBC_GATA3_group$GATA3_group,levels = c("low","high"))

p2 <- ggboxplot(METABRIC_TPBC_GATA3_group, x="GATA3_group", y="TMB_NONSYNONYMOUS", 
                color="GATA3_group",
                add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_location,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1), method = "wilcox.test",
                                        symbols = c("***", "**", "*", "NS"))) +
  labs(x = "GATA3 Group") +
  labs(y = "TMB") 

p1 + p2
# MMR ----
MMR_genes <- c("MLH1","MSH2","MSH6","PMS2")
groups <- TCGA_TPBC_data_tidy_group %>% select(SAMPLE_ID,GATA3_group)
load(file = "./GDCdata/TCGA_BRCA_tpm_anno.rda")
exp_mmr <- log(TCGA_BRCA_tpm_anno[MMR_genes,] + 0.001) %>% t() %>% as.data.frame() %>%
  rownames_to_column("SAMPLE_ID") %>% mutate(SAMPLE_ID = str_sub(SAMPLE_ID,1,15)) %>% inner_join(groups, "SAMPLE_ID")

exp_mmr$GATA3_group <- factor(exp_mmr$GATA3_group,levels = c("low", "high"))

df <- exp_mmr %>% 
  pivot_longer(cols = 2:5,
               names_to= "gene_symbol",
               values_to = "expression_value")

ggplot(df, 
       aes(x = gene_symbol,
           y = expression_value)) +
  geom_boxplot(aes(color = GATA3_group),
               outlier.shape = NA) +
  theme_classic(base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1, 
                                   colour = "black"))+
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
