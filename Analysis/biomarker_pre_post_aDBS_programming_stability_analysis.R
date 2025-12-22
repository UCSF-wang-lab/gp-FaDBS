library(tidyverse)
library(lmerTest)
library(emmeans)
library(rstatix)

#### Helper Function ####
outlier_threshold <- function(df,variableName,quantileVal)
{
  iqr_val <- IQR(df[[variableName]],na.rm = TRUE)
  quantile_val <- quantile(df[[variableName]],probs = quantileVal,na.rm = TRUE)
  ifelse(quantileVal>.50, return(quantile_val+1.5*iqr_val), return(quantile_val-1.5*iqr_val))
}

##### Load data ####
ddata <- read.csv("/Users/USER/Documents/all_patients_biomarker_power.csv")
data <- data %>% mutate(VisitName = if_else(str_detect(VisitName,"dbsOpt|ptClinic"),"Post optimized DBS",if_else(str_detect(VisitName,"aDBS"),"aDBS","Pre optimized DBS"))) %>% 
  mutate(VisitName = as_factor(VisitName)) %>% 
  mutate(VisitName = fct_relevel(VisitName,c("Pre optimized DBS","Post optimized DBS","aDBS"))) %>% 
  mutate(Side = ifelse(str_detect(Side,"left"),"Left","Right")) %>% 
  mutate(Side = as_factor(Side)) %>% 
  mutate(Side = fct_relevel(Side,c("Left","Right"))) %>% 
  mutate(SubjectID = as_factor(SubjectID)) %>% 
  mutate(SubjectID = fct_relevel(SubjectID,c("P1","P2","P3","P4","P5")))

P1 <- data %>% filter(SubjectID == "P1",MedState == "ON") %>% 
  group_by(VisitName) %>% filter(PowerBand>outlier_threshold(.,"PowerBand",0.25),PowerBand<outlier_threshold(.,"PowerBand",0.75))

P2 <- data %>% filter(SubjectID == "P2",MedState == "ON") %>% 
  group_by(VisitName) %>% filter(PowerBand>outlier_threshold(.,"PowerBand",0.25),PowerBand<outlier_threshold(.,"PowerBand",0.75))

p3 <- data %>% filter(SubjectID == "P3",MedState == "ON") %>% 
  group_by(VisitName) %>% filter(PowerBand>outlier_threshold(.,"PowerBand",0.25),PowerBand<outlier_threshold(.,"PowerBand",0.75))

P4 <- data %>% filter(SubjectID == "P4",MedState == "ON") %>% 
  group_by(VisitName) %>% filter(PowerBand>outlier_threshold(.,"PowerBand",0.25),PowerBand<outlier_threshold(.,"PowerBand",0.75))

P4_right_blank <- setNames(data.frame(matrix(ncol = ncol(data), nrow = 0)), colnames(data))
P4_right_blank[1,] <- NaN
P4_right_blank[1,"SubjectID"] <- "P4"
P4_right_blank[1,"VisitName"] <- "Pre optimized DBS"
P4_right_blank[1,"MedState"] <- "ON"
P4_right_blank[1,"StimState"] <- "OFF"
P4_right_blank[1,"Side"] <- "Right"
P4_right_blank[1,"FS"] <- 500
P4_right_blank[1,"NFFT"] <- 256
P4_right_blank[1,"OverlapPercent"] <- 50
P4_right_blank[1,"FFTInterval"] <- 256
P4_right_blank[1,"Bitshift"] <- 3

P5 <- data %>% filter(SubjectID == "P5",MedState == "ON") %>% 
  group_by(VisitName) %>% filter(PowerBand>outlier_threshold(.,"PowerBand",0.25),PowerBand<outlier_threshold(.,"PowerBand",0.75))

combined_data <- bind_rows(P1,P2,P3,P4,P5) %>% filter(VisitName %in% c("Post optimized DBS","aDBS"))
combined_data <- bind_rows(P1,P2,P3,P4,P5)

##### Stats ####
PowerBand_normality <- combined_data %>% group_by(SubjectID,VisitName,Side) %>% sample_n(size = min(n(),5000),replace = FALSE) %>% shapiro_test(PowerBand)
PowerBandStats <- combined_data %>% group_by(SubjectID,Side) %>% kruskal_test(PowerBand ~ VisitName) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
PowerBandStatsPostHoc <- combined_data %>% group_by(SubjectID,Side) %>% dunn_test(PowerBand ~ VisitName) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
PowerBandEffectSize <- combined_data %>% group_by(SubjectID,Side) %>% cohens_d(PowerBand ~ VisitName)
