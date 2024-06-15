rm(list = ls())
load("./data/BRCA_TPBC_data.rda")
TCGA_BRCA_aclbi <- read_csv("data/TCGA_BRCA_aclbi.csv") %>% rename(PATIENT_ID = `...1`)

#  ----
TCGA_TPBC_data_tidy <- TCGA_BRCA_aclbi %>% filter(PATIENT_ID %in% TCGA_TPBC_data$PATIENT_ID) %>% 
  select(PATIENT_ID,Age,T,N,M,Stage,newTumor,Radiation,Neoadjuvant,Therapy) %>% left_join(TCGA_TPBC_data,"PATIENT_ID")

TCGA_TPBC_data_tidy <- TCGA_TPBC_data_tidy %>% 
  mutate(
    Stage = case_when(
      Stage == "I" ~ "I+II",
      Stage == "II" ~ "I+II",
      Stage == "III" ~ "III+IV",
      Stage == "IV" ~ "III+IV"
  )) %>% 
  mutate(
    T = case_when(
      T == "T1" ~ "T1+T2",
      T == "T2" ~ "T1+T2",
      T == "T3" ~ "T3+T4",
      T == "T4" ~ "T3+T4"
    )) %>% 
  mutate(
    N = case_when(
      N == "N0" ~ "N0+N1",
      N == "N1" ~ "N0+N1",
      N == "N2" ~ "N2+N3",
      N == "N3" ~ "N2+N3"))

TCGA_TPBC_data_tidy$Stage <- factor(TCGA_TPBC_data_tidy$Stage, levels = c("I+II","III+IV"))
TCGA_TPBC_data_tidy$T <- factor(TCGA_TPBC_data_tidy$T, levels = c("T1+T2","T3+T4"))
TCGA_TPBC_data_tidy$N <- factor(TCGA_TPBC_data_tidy$N, levels = c("N0+N1","N2+N3"))
TCGA_TPBC_data_tidy$M <- factor(TCGA_TPBC_data_tidy$M, levels = c("M0","M1"))

# KM ----
library(survival)
library(survminer)

res.cut <- surv_cutpoint(TCGA_TPBC_data_tidy, 
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
           ylab = "OS Survival (%)",
           xlab = " Time (Months)",
           legend.title = "Group",
           legend.labs = c("GATA3-High","GATA3-Low"),
           risk.table = F,
           tables.height = 0.3,
           palette = c("#ED0000FF","#0099B4FF"),
           ggtheme = theme_classic(),
           tables.theme = theme_survminer()
) 

# cox ----

TCGA_TPBC_data_tidy_group <- cbind(TCGA_TPBC_data_tidy,as.data.frame(group)) %>% rename(GATA3_group = group) 
TCGA_TPBC_data_tidy_group$GATA3_group <- factor(TCGA_TPBC_data_tidy_group$GATA3_group, levels = c("low","high"))
save(TCGA_TPBC_data_tidy_group, file = "./data/TCGA_TPBC_data_tidy_group.rda")

df_1 <- TCGA_TPBC_data_tidy_group %>% select(OS_MONTHS, OS_STATUS,AGE,GATA3_group,Stage,T,N,M,MSI_SENSOR_SCORE,TMB_NONSYNONYMOUS)

Uni_cox <- data.frame()
for(i in colnames(df_1[,3:ncol(df_1)])){
  cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ df_1[,i], data = df_1)
  coxSummary <- summary(cox)
  Uni_cox  <- rbind(Uni_cox , 
                    cbind(gene_name = i,
                          HR = coxSummary[["coefficients"]][,"exp(coef)"],
                          z = coxSummary[["coefficients"]][,"z"],
                          Pvalue = coxSummary[["coefficients"]][,"Pr(>|z|)"],
                          Lower = coxSummary[["conf.int"]][,"lower .95"],
                          Upper = coxSummary[["conf.int"]][,"upper .95"]))
}


multi_cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~GATA3_group+AGE+Stage, data = df_1)
# multi_cox <- step(multi_cox, direction = "both")
multiSummary <- summary(multi_cox)

multioutput <- cbind.data.frame(
  coef = multiSummary[["coefficients"]][,"coef"],
  HR = multiSummary[["conf.int"]][,"exp(coef)"],
  Pvalue = multiSummary[["coefficients"]][,"Pr(>|z|)"],
  Lower = multiSummary[["conf.int"]][,"lower .95"],
  Upper = multiSummary[["conf.int"]][,"upper .95"])

write.csv(Uni_cox, file = "./cox/unicox_GATA3.csv")
write.csv(multioutput, file = "./cox/multicox_GATA3.csv")
save(Uni_cox,multioutput, file = "./cox/uni_multicox.rda")

#  ----
library(forestploter)
# 
unicox_GATA3_forest <- read_csv("cox/unicox_GATA3_forest.csv")
unicox_GATA3_forest$Characteristics <- ifelse(is.na(unicox_GATA3_forest$`Total (N)`),
                                              unicox_GATA3_forest$Characteristics,
                                              paste0("    ", unicox_GATA3_forest$Characteristics ))
unicox_GATA3_forest$HR <- as.numeric(unicox_GATA3_forest$HR)

unicox_GATA3_forest$`Total (N)` <- ifelse(is.na(unicox_GATA3_forest$`Total (N)`),"",unicox_GATA3_forest$`Total (N)`)

unicox_GATA3_forest$`HR (95%CI)` <- ifelse(is.na(unicox_GATA3_forest$`HR (95%CI)`),"",unicox_GATA3_forest$`HR (95%CI)`)

unicox_GATA3_forest$`P value` <- ifelse(is.na(unicox_GATA3_forest$`P value`),"",unicox_GATA3_forest$`P value`)

# 
unicox_GATA3_forest$` ` <- paste(rep(" ", 12), collapse = " ")

# 
multicox_GATA3_forest <- read_csv("cox/multicox_GATA3_forest.csv")
multicox_GATA3_forest$Characteristics <- ifelse(is.na(multicox_GATA3_forest$`Total (N)`),
                                              multicox_GATA3_forest$Characteristics,
                                              paste0("    ", multicox_GATA3_forest$Characteristics ))
multicox_GATA3_forest$HR <- as.numeric(multicox_GATA3_forest$HR)

multicox_GATA3_forest$`Total (N)` <- ifelse(is.na(multicox_GATA3_forest$`Total (N)`),"",multicox_GATA3_forest$`Total (N)`)

multicox_GATA3_forest$`HR (95%CI)` <- ifelse(is.na(multicox_GATA3_forest$`HR (95%CI)`),"",multicox_GATA3_forest$`HR (95%CI)`)

multicox_GATA3_forest$`P value` <- ifelse(is.na(multicox_GATA3_forest$`P value`),"",multicox_GATA3_forest$`P value`)

multicox_GATA3_forest$` ` <- paste(rep(" ", 12), collapse = " ")

#  ----

tm <- forest_theme(
  base_size = 10,        

  ci_pch = 15,           
  ci_col = "blue4",    
  ci_fill = "blue4",      
  ci_alpha = 0.8,        
  ci_lty = 1,            
  ci_lwd = 1.5,         
  ci_Theight = 0.2,      
  

  refline_lwd = 2,         
  refline_lty = "dashed",  
  refline_col = "grey20",  
  

  vertline_lwd = 1,        
  vertline_lty = "dashed", 
  vertline_col = "grey20",  
  

  footnote_cex = 0.6,            
  footnote_fontface = "italic", 
  footnote_col = "red4"         
)


p2_1 <- forest(unicox_GATA3_forest[,c(1, 2, 8, 6, 7)],  
       est = unicox_GATA3_forest$HR,           
       lower = unicox_GATA3_forest$Lower,    
       upper = unicox_GATA3_forest$Upper,   
       ci_column = 3,             
       ref_line = 1,               
       arrow_lab = c("Low risk", "High Risk"),   
       xlim = c(-1, 10),           
       ticks_at = c(-0.5, 1, 3, 5, 7),  
       theme = tm                  
       )  
p2 <- edit_plot(p2_1, row = c(1,4,7), which = "background",
                gp = gpar(fill = "yellow"))

p3 <- forest(multicox_GATA3_forest[,c(1, 2, 8, 6, 7)],  
               est = multicox_GATA3_forest$HR,           
               lower = multicox_GATA3_forest$Lower,   
               upper = multicox_GATA3_forest$Upper,   
               ci_column = 3,             
               ref_line = 1,               
               arrow_lab = c("Low risk", "High Risk"),  
               xlim = c(-1, 10),          
               ticks_at = c(-0.5, 1, 3, 5, 7),   
               theme = tm                
)  





# tableone ----
library(compareGroups)

restab <- compareGroups(GATA3_group ~ AGE + Stage + T + N + M + MSI_SCORE_MANTIS + TMB_NONSYNONYMOUS, 
                        data = df_1, 
                        method = c(AGE=NA,MSI_SCORE_MANTIS=NA,TMB_NONSYNONYMOUS=NA),
                        ref = 1
) 

creatab <- createTable(restab)

export2word(creatab, file="./TCGA_TPBC_table1.docx") # 只能在根目录输出三线表格式

# METABRIC ----
load("./data/METABRIC_exp_GATA3_TPBC.rda")
res.cut <- surv_cutpoint(METABRIC_exp_GATA3_TPBC, time = "OS_MONTHS", 
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
           ylab = "OS Survival (%)",
           xlab = " Time (Months)",
           legend.title = "Group",
           legend.labs = c("GATA3-High","GATA3-Low"),
           risk.table = F,
           tables.height = 0.3,
           palette = c("#ED0000FF","#0099B4FF"),
           ggtheme = theme_classic(),
           tables.theme = theme_survminer()
) 

# 分组分析 ----
# 病理学分期
compaired <- list(c("I+II", "III+IV")) 
df_stage <- TCGA_TPBC_data_tidy_group %>% select(Stage, GATA3) %>% drop_na(Stage)
p5 <- ggboxplot(df_stage,
          x = "Stage", y = "GATA3",
          color = "Stage",
          add = "jitter", 
          add.params = list(color = "Stage",alpha=0.5,shape = 16)
) +   #添加抖动点
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns"))) +
  scale_color_aaas() + 
  theme(legend.position = "none",
        axis.line = element_line(size = 1.0)) +
  labs(x = "Stage Group") +
  labs(y = "GATA3 mRNA expression") 

# T分期
compaired <- list(c("T1+T2", "T3+T4")) 
df_T <- TCGA_TPBC_data_tidy_group %>% select("T", GATA3) %>% drop_na("T")
p6 <- ggboxplot(df_T,
                x = "T", y = "GATA3",
                color = "T",
                add = "jitter", 
                add.params = list(color = "T",alpha=0.5,shape = 16)
) +   #添加抖动点
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns"))) +
  scale_color_aaas() + 
  theme(legend.position = "none",
        axis.line = element_line(size = 1.0)) +
  labs(x = "T stage") +
  labs(y = "GATA3 mRNA expression") 

# N分期
compaired <- list(c("N0+N1", "N2+N3")) 
df_N <- TCGA_TPBC_data_tidy_group %>% select("N", GATA3) %>% drop_na("N")
p7 <- ggboxplot(df_N,
                x = "N", y = "GATA3",
                color = "N",
                add = "jitter", 
                add.params = list(color = "N",alpha=0.5,shape = 16)
) +   #添加抖动点
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns"))) +
  scale_color_aaas() + 
  theme(legend.position = "none",
        axis.line = element_line(size = 1.0)) +
  labs(x = "N stage") +
  labs(y = "GATA3 mRNA expression") 

