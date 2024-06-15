rm(list = ls())
library(readxl)
library(ggstatsplot)
library(gmodels)
library(ggsci)
TPBC_H_score_R <- read_excel("TPBC_H_score_R.xlsx") 
TPBC_H_score_R$groups <- factor(TPBC_H_score_R$groups, levels = c("L","H"))

my_comparisons_group <- list(c("L", "H"))

p1 <- ggboxplot(TPBC_H_score_R, x="groups", y="TILs", 
                color="groups",
                add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))

CrossTable(TPBC_H_score_R$groups, TPBC_H_score_R$Stage, digits = 4, expected = T, chisq = T, fisher = T, mcnemar = F, format = "SPSS")

p2 <- ggbarstats(TPBC_H_score_R, Stage,groups,title = NULL,results.subtitle = FALSE) +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_text(aes(label = "*"), 
            x = 1.5,
            y = 1.05,
            vjust = -0.5, 
            size = 4) +
  geom_segment(aes(x = 1.0, xend = 2.0, y = 1.05, yend = 1.05), linetype = "solid") 

p1 / p2 


p3 <- ggboxplot(TPBC_H_score_R, x="groups", y="Age", 
                color="groups",
                add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))

p4 <- ggboxplot(TPBC_H_score_R, x="groups", y="`Ki-67`", 
                color="groups",
                add = "jitter",size = 1) + 
  scale_color_aaas() +
  stat_compare_means(comparisons = my_comparisons_group,
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.001,0.01,0.05,1),
                                        symbols = c("***", "**", "*", "NS")))
CrossTable(TPBC_H_score_R$groups, TPBC_H_score_R$histological_grade, digits = 4, expected = T, chisq = T, fisher = T, mcnemar = F, format = "SPSS")
p5 <- ggbarstats(TPBC_H_score_R, histological_grade ,groups,title = NULL,results.subtitle = FALSE ) +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_text(aes(label = "NS"), 
            x = 1.5,
            y = 1.05,
            vjust = -0.5, 
            size = 4) +
  geom_segment(aes(x = 1.0, xend = 2.0, y = 1.05, yend = 1.05), linetype = "solid") 

p3 + p4 + p5 

library(cowplot)
plot_grid( p3, p4, p5, nrow = 1, align = "vh")

# groups ~ T-Stage 
CrossTable(TPBC_H_score_R$groups, TPBC_H_score_R$T_Stage, digits = 4, expected = T, chisq = T, fisher = T, mcnemar = F, format = "SPSS")
p6 <- ggbarstats(TPBC_H_score_R, T_Stage,groups,title = NULL,results.subtitle = FALSE) +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_text(aes(label = "NS"), 
            x = 1.5,
            y = 1.05,
            vjust = -0.5, 
            size = 4) +
  geom_segment(aes(x = 1.0, xend = 2.0, y = 1.05, yend = 1.05), linetype = "solid") 

# groups ~ N-Stage 
CrossTable(TPBC_H_score_R$groups, TPBC_H_score_R$N_Stage, digits = 4, expected = T, chisq = T, fisher = T, mcnemar = F, format = "SPSS")
p7 <- ggbarstats(TPBC_H_score_R, N_Stage,groups,title = NULL,results.subtitle = FALSE) +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_text(aes(label = "NS"), 
            x = 1.5,
            y = 1.05,
            vjust = -0.5, 
            size = 4) +
  geom_segment(aes(x = 1.0, xend = 2.0, y = 1.05, yend = 1.05), linetype = "solid") 

(p3 | p4) / (p5 + p6 + p7)