# mutation ----
library(maftools)
library(RColorBrewer)
library(paletteer)
rm(list = ls())
load(file = "./GDCdata/TCGA_BRCA_mutation_filter.rda")
load(file = "./data/TCGA_TPBC_data_tidy_group.rda")
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
Variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
             "In_Frame_Ins","Missense_Mutation","Nonsense_Mutation",
             "Nonstop_Mutation", "Splice_Region","Splice_Site",
             "Translation_Start_Site")
TCGA_TPBC_mutation <- TCGA_BRCA_mutation_filter %>% filter(Tumor_Sample_Barcode %in% TCGA_TPBC_data_tidy_group$SAMPLE_ID)
high_group_list <- TCGA_TPBC_data_tidy_group %>% filter(GATA3_group == "high")
low_group_list <- TCGA_TPBC_data_tidy_group %>% filter(GATA3_group == "low")
TCGA_TPBC_mutation_high <- TCGA_BRCA_mutation_filter %>% filter(Tumor_Sample_Barcode %in% high_group_list$SAMPLE_ID)
TCGA_TPBC_mutation_low <- TCGA_BRCA_mutation_filter %>% filter(Tumor_Sample_Barcode %in% low_group_list$SAMPLE_ID)

high_group_maf <- read.maf(maf = TCGA_TPBC_mutation_high,
                           vc_nonSyn = Variants)
low_group_maf <- read.maf(maf = TCGA_TPBC_mutation_low,
                          vc_nonSyn = Variants)

p1 <- coBarplot(m1 = high_group_maf, m2 = low_group_maf, m1Name = "High", m2Name = "Low",colors = vc_cols,borderCol = "black")


# lollipopPlot2(m1 = high_group_maf, m2 = low_group_maf, gene = "GATA3", AACol1 = "HGVSp", AACol2 = "HGVSp", m1_name = "High", m2_name = "Low",colors = vc_cols)

# CNV ----

#  ----
load(file = "data/illuminaMethyl450_hg38_GDC.rda")
load(file = "data/TCGA_TPBC_meth.rda")
meth450_GATA3_anno <- illuminaMethyl450_hg38_GDC %>% filter(gene == "GATA3")
TCGA_TPBC_GATA3_meth <- TCGA_TPBC_meth[meth450_GATA3_anno$id,] %>% t() %>% as.data.frame() %>% mutate(GATA3_meth_level = rowMeans(.)) %>% rownames_to_column("SAMPLE_ID")

TCGA_TPBC_GATA3_meth_exp <- TCGA_TPBC_GATA3_meth %>% inner_join(TCGA_TPBC_data_tidy_group,"SAMPLE_ID") %>% select(GATA3_meth_level,GATA3,GATA3_group) 

# 
library(ggpubr)
p2 <- ggscatter(TCGA_TPBC_GATA3_meth_exp, x = "GATA3", y = "GATA3_meth_level", 
          color = "GATA3_group",fill = "lightgray", size = 2.5,
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = 0.1),
          cor.method = "spearman") 

my_comparisons_group <- list(c("low", "high"))
p3 <- ggboxplot(TCGA_TPBC_GATA3_meth_exp, x="GATA3_group", y="GATA3_meth_level", 
          color="GATA3_group",
          add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))


#  ----

data_rppa_zscores <- read_delim("data/brca_tcga_pan_can_atlas_2018/data_rppa_zscores.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)

TCGA_RPPA_GATA3 <- data_rppa_zscores %>% filter( Composite.Element.REF %in% c("GATA3|GATA3") ) %>% column_to_rownames("Composite.Element.REF") %>% t() %>% as.data.frame() %>% rownames_to_column("SAMPLE_ID") %>% dplyr::rename("GATA3_RPPA" = "GATA3|GATA3")

groups <- TCGA_TPBC_data_tidy_group %>% select(SAMPLE_ID,GATA3_group)

TCGA_RPPA_GATA3_group <- TCGA_RPPA_GATA3 %>% inner_join(groups,"SAMPLE_ID")


# ggpubr ----
my_comparisons_RPPA <- list(c("low", "high"))
p4 <- ggboxplot(TCGA_RPPA_GATA3_group, x = "GATA3_group" , y = "GATA3_RPPA",
          color = "GATA3_group",
          add = "jitter", size = 1) +
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_RPPA,
                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))



TCGA_TPBC_maf <- read.maf(maf = TCGA_TPBC_mutation,vc_nonSyn = Variants)
p5 <- oncoplot(maf = TCGA_TPBC_maf, top = 20, colors = vc_cols,showTitle = F)


save(p1,p2,p3,p4, file =  "./Figure/Fig6.rda") 
save(p5, file = "./Figure/Fig6p5.rda")