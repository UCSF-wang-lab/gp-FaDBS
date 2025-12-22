library(tidyverse)
library(lmerTest)
library(emmeans)
library(rstatix)

#### Helper functions ####
outlier_threshold <- function(df,variableName,quantileVal)
{
  iqr_val <- IQR(df[[variableName]],na.rm = TRUE)
  quantile_val <- quantile(df[[variableName]],probs = quantileVal,na.rm = TRUE)
  ifelse(quantileVal>.50, return(quantile_val+1.5*iqr_val), return(quantile_val-1.5*iqr_val))
}

#### Load in data, filter, and format ####
# P1
P1_gait_data <- read.csv('/Users/USER/Documents/P1_aggregate_gait_metrics.csv')

P1_clinicOpt_to_aDBS <- P1_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P2
P2_gait_data <- read.csv('Users/USER/Documents/P2_aggregate_gait_metrics.csv');

P2_clinicOpt_to_aDBS <- P2_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P3
P3_gait_data <- read.csv('Users/USER/Documents/P3_aggregate_gait_metrics.csv')

P3_clinicOpt_to_aDBS <- P3_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P4
P4_gait_data <- read.csv('Users/USER/Documents/P4_aggregate_gait_metrics.csv')

P4_clinicOpt_to_aDBS <- P4_gait_data %>% filter(VisitName %in% c("dbsOptUnilateral","aDB")) %>% mutate(Side = as_factor(Side))

# P5
P5_gait_data <- read.csv('Users/USER/Documents/P5_aggregate_gait_metrics.csv')

P5_clinicOpt_to_aDBS <- P5_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# Combine data sets
clinicOpt_to_aDBS <- bind_rows(P1_clinicOpt_to_aDBS,P2_clinicOpt_to_aDBS,P3_clinicOpt_to_aDBS,P4_clinicOpt_to_aDBS,P5_clinicOpt_to_aDBS)

clinicOpt_to_aDBS_w_symmetry <- clinicOpt_to_aDBS %>% group_by(SubjectID,VisitName,SettingNum,MedState,StimState,TrialNum,GaitCycle) %>% 
  mutate(StepLengthSymmetry = (first(StepLength)-last(StepLength))/sum(StepLength),
         StepTimeSymmetry = (first(StepTime)-last(StepTime))/sum(StepTime)) %>% 
  mutate(VisitName = ifelse(str_detect(VisitName,"dbsOpt"),"cDBS","aDBS")) %>% 
  mutate(VisitName = as_factor(VisitName)) %>% 
  mutate(VisitName = fct_relevel(VisitName,c("cDBS","aDBS")))

#### Variance and Coefficient of variation calculation ####
summary_df <- clinicOpt_to_aDBS %>% group_by(SubjectID,VisitName,Side) %>% 
  summarise(StepLengthMean = mean(StepLength,na.rm = TRUE),
            StepTimeMean = mean(StepTime,na.rm = TRUE),
            StepLengthVar = var(StepLength,na.rm = TRUE),
            StepTimeVar = var(StepTime,na.rm = TRUE),
            StepLengthCV = sd(StepLength,na.rm = TRUE)/mean(StepLength,na.rm = TRUE),
            StepTimeCV = sd(StepTime,na.rm = TRUE)/mean(StepTime,na.rm = TRUE)) %>% 
  mutate(VisitName = ifelse(str_detect(VisitName,"dbsOpt"),"cDBS","aDBS")) %>% 
  mutate(VisitName = as_factor(VisitName)) %>% 
  mutate(VisitName = fct_relevel(VisitName,c("cDBS","aDBS")))

group_step_length = clinicOpt_to_aDBS_w_symmetry %>% 
  filter(Side == "L") %>% 
  group_by(SubjectID) %>% 
  filter(StepLength>outlier_threshold(.,"StepLength",0.25),
         StepLength<outlier_threshold(.,"StepLength",0.75),
         StepLengthSymmetry < 0.25,StepLengthSymmetry > -0.25) %>% 
  mutate(StepLengthSymmetry = abs(StepLengthSymmetry))

group_step_time = clinicOpt_to_aDBS_w_symmetry %>% 
  filter(Side == "L") %>% 
  group_by(SubjectID) %>% 
  filter(StepTime>outlier_threshold(.,"StepTime",0.25),
         StepTime<outlier_threshold(.,"StepTime",0.75),
         StepTimeSymmetry < 0.25,StepTimeSymmetry > -0.25) %>% 
  mutate(StepTimeSymmetry = abs(StepTimeSymmetry))

group_step_CV = summary_df %>% 
  select(SubjectID,VisitName,Side,StepLengthCV,StepTimeCV)

#### Stats ####
StepLength_normality <- clinicOpt_to_aDBS_w_symmetry %>% group_by(SubjectID,Side) %>% filter(StepLength>outlier_threshold(.,"StepLength",0.25),StepLength<outlier_threshold(.,"StepLength",0.75)) %>% shapiro_test(StepLength)
StepLengthStats <- clinicOpt_to_aDBS_w_symmetry %>% group_by(SubjectID,Side) %>% filter(StepLength>outlier_threshold(.,"StepLength",0.25),StepLength<outlier_threshold(.,"StepLength",0.75)) %>% wilcox_test(StepLength ~ VisitName) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
StepLengthDiff <- clinicOpt_to_aDBS_w_symmetry %>% group_by(SubjectID,Side,VisitName) %>% filter(StepLength>outlier_threshold(.,"StepLength",0.25),StepLength<outlier_threshold(.,"StepLength",0.75)) %>% summarise(StepLengthMean = mean(StepLength)) 
StepLengthDiff2 <- StepLengthDiff %>% group_by(SubjectID,Side) %>% summarise(StepLengthPercentChange = ((StepLengthMean[2]-StepLengthMean[1])/StepLengthMean[1])*100)

StepTime_normality <- clinicOpt_to_aDBS_w_symmetry %>% group_by(SubjectID,Side) %>% filter(StepTime>outlier_threshold(.,"StepTime",0.25),StepTime<outlier_threshold(.,"StepTime",0.75)) %>% shapiro_test(StepTime) 
StepTimeStats <- clinicOpt_to_aDBS_w_symmetry %>% group_by(SubjectID,Side) %>% filter(StepTime>outlier_threshold(.,"StepTime",0.25),StepTime<outlier_threshold(.,"StepTime",0.75)) %>% wilcox_test(StepTime ~ VisitName) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
StepTimeDiff <- clinicOpt_to_aDBS_w_symmetry %>% group_by(SubjectID,Side,VisitName) %>% filter(StepTime>outlier_threshold(.,"StepTime",0.25),StepTime<outlier_threshold(.,"StepTime",0.75)) %>% summarise(StepTimeMean = mean(StepTime)) 
StepTimeDiff2 <- StepTimeDiff %>% group_by(SubjectID,Side) %>% summarise(StepTimePercentChange = ((StepTimeMean[2]-StepTimeMean[1])/StepTimeMean[1])*100)

StepLengthSymm_normality <- clinicOpt_to_aDBS_w_symmetry %>% filter(Side == "L") %>% group_by(SubjectID) %>% filter(StepLength>outlier_threshold(.,"StepLength",0.25),StepLength<outlier_threshold(.,"StepLength",0.75),StepLengthSymmetry < 0.25,StepLengthSymmetry > -0.25) %>% shapiro_test(StepLengthSymmetry) 
StepLengthSymmStats <- clinicOpt_to_aDBS_w_symmetry %>% filter(Side == "L") %>% group_by(SubjectID) %>% filter(StepLength>outlier_threshold(.,"StepLength",0.25),StepLength<outlier_threshold(.,"StepLength",0.75),StepLengthSymmetry < 0.25,StepLengthSymmetry > -0.25) %>% wilcox_test(StepLengthSymmetry ~ VisitName) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

StepTimeSymm_normality <- clinicOpt_to_aDBS_w_symmetry %>% filter(Side == "L") %>% group_by(SubjectID) %>% filter(StepTime>outlier_threshold(.,"StepTime",0.25),StepTime<outlier_threshold(.,"StepTime",0.75),StepTimeSymmetry < 0.25,StepTimeSymmetry > -0.25) %>% shapiro_test(StepTimeSymmetry) 
StepTimeSymmStats <- clinicOpt_to_aDBS_w_symmetry %>% filter(Side == "L") %>% group_by(SubjectID) %>% filter(StepTime>outlier_threshold(.,"StepTime",0.25),StepTime<outlier_threshold(.,"StepTime",0.75),StepTimeSymmetry < 0.25,StepTimeSymmetry > -0.25) %>% wilcox_test(StepTimeSymmetry ~ VisitName) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

# Aggregate subject symmetry
group_stats_step_length = group_step_length %>% ungroup() %>% 
  wilcox_test(StepLengthSymmetry ~ VisitName) %>% 
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
group_step_length_diff <- group_step_length %>% group_by(VisitName) %>% summarize(IQRMean = mean(IQR(StepLengthSymmetry)), StepLengthSymmMean = mean(StepLengthSymmetry))

group_stats_step_time = group_step_time %>% ungroup() %>% 
  wilcox_test(StepTimeSymmetry ~ VisitName) %>% 
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
group_step_time_diff <- group_step_time %>% group_by(VisitName) %>% summarize(IQRMean = mean(IQR(StepTimeSymmetry)), StepTimeSymmMean = mean(StepTimeSymmetry))

group_stats_step_length_var = group_step_CV %>% ungroup() %>% group_by(Side) %>%
  wilcox_test(StepLengthCV ~ VisitName) %>%
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

group_stats_step_time_var = group_step_CV %>% ungroup() %>% group_by(Side) %>%
  wilcox_test(StepTimeCV ~ VisitName) %>%
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
