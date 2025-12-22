library(tidyverse)
library(rstatix)
library(lmerTest)
library(ordinal)
library(FSA)

##### Helper functions #####
outlier_threshold <- function(df,variableName,quantileVal)
{
  iqr_val <- IQR(df[[variableName]], na.rm = TRUE)
  quantile_val <- quantile(df[[variableName]], probs = quantileVal, na.rm = TRUE)
  
  if (quantileVal > 0.5) {
    return(quantile_val + 1.5 * iqr_val)
  } else {
    return(quantile_val - 1.5 * iqr_val)
  }
}

##### In-clinic gait metrics #####
# Load data
clinic_data <-read.csv('/Users/USER/Documents/long-term_aDBS_gait_metrics.csv')
clinic_data <- clinic_data %>% 
  mutate(SubjectID = as.factor(SubjectID),
         SubjectID = fct_relevel(SubjectID,c("P2","P3","P4")),
         DBSCondition = factor(case_when(str_detect(DBSCondition,"Clinical") ~ "cDBS",str_detect(DBSCondition,"Ramp-Up") ~ "RU-aDBS",str_detect(DBSCondition,"Ramp-Down") ~ "RD-aDBS"),levels = c("cDBS","RU-aDBS","RD-aDBS")),
         GaitCycle = as.factor(GaitCycle),
         WalkType = as.factor(WalkType))

clinic_data_filt <- clinic_data %>% 
  group_by(SubjectID,DBSCondition) %>% 
  filter({
    lower <- outlier_threshold(cur_data(), "StepLength_L", 0.25)
    upper <- outlier_threshold(cur_data(), "StepLength_L", 0.75)
    StepLength_L >= lower & StepLength_L <= upper
  },
  {
    lower <- outlier_threshold(cur_data(), "StepLength_R", 0.25)
    upper <- outlier_threshold(cur_data(), "StepLength_R", 0.75)
    StepLength_R >= lower & StepLength_R <= upper
  },
  {
    lower <- outlier_threshold(cur_data(), "StepTime_L", 0.25)
    upper <- outlier_threshold(cur_data(), "StepTime_L", 0.75)
    StepTime_L >= lower & StepTime_L <= upper
  },
  {
    lower <- outlier_threshold(cur_data(), "StepTime_R", 0.25)
    upper <- outlier_threshold(cur_data(), "StepTime_R", 0.75)
    StepTime_R >= lower & StepTime_R <= upper
  })

gait_metrics <- clinic_data_filt %>% 
  mutate(StepLengthSymm_L = StepLengthSymm, StepTimeSymm_L = StepTimeSymm) %>% 
  select(-StepLengthSymm,-StepTimeSymm) %>% 
  pivot_longer(cols = ends_with("_L") | ends_with("_R"),names_to = "Metric",values_to = "Value") %>% 
  mutate(Side = case_when(grepl("_L$", Metric) ~ "Left",grepl("_R$", Metric) ~ "Right")) %>% 
  mutate(Metric = sub("_[LR]$","",Metric)) %>% 
  select(-WalkType)

summary_clinic_data <- gait_metrics %>% group_by(SubjectID,DBSCondition,Metric,Side) %>% 
  summarise(mean = mean(abs(Value/100),na.rm = TRUE),
            lower = mean(abs(Value/100),na.rm = TRUE)-(sd(abs(Value/100),na.rm = TRUE)/sqrt(n())),
            upper = mean(abs(Value/100),na.rm = TRUE)+(sd(abs(Value/100),na.rm = TRUE)/sqrt(n())),
            var = var(abs(Value/100),na.rm = TRUE),sd = sd(abs(Value/100),na.rm = TRUE),
            cv = sd(Value,na.rm = TRUE)/mean(Value,na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(SubjectID,Metric,Side) %>% 
  mutate(percentChange_C_to_RU = ((mean[2]-mean[1])/mean[1])*100,
         percentChange_C_to_RD = ((mean[3]-mean[1])/mean[1])*100,
         percentChange_RU_to_RD = ((mean[3]-mean[2])/mean[2])*100)

# Stats
SL_stat <- gait_metrics %>% filter(Metric == "StepLength") %>% group_by(SubjectID,Side) %>% kruskal_test(Value ~ DBSCondition) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
SL_stat2 <- gait_metrics %>% filter(Metric == "StepLength") %>% group_by(SubjectID,Side) %>% dunn_test(Value ~ DBSCondition) %>% mutate(signif = if_else(p.adj<0.05,if_else(p.adj<0.001,if_else(p.adj<0.0001,"***","**"),"*"),NA))

ST_stat <- gait_metrics %>% filter(Metric == "StepTime") %>% group_by(SubjectID,Side) %>% kruskal_test(Value ~ DBSCondition) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))
ST_stat2 <- gait_metrics %>% filter(Metric == "StepTime") %>% group_by(SubjectID,Side) %>% dunn_test(Value ~ DBSCondition) %>% mutate(signif = if_else(p.adj<0.05,if_else(p.adj<0.001,if_else(p.adj<0.0001,"***","**"),"*"),NA))

SL_Symm_stat <- kruskal.test(Value~DBSCondition,data = gait_metrics %>% filter(Metric == "StepLengthSymm"))
SL_Symm_stat2 <- dunnTest(Value~DBSCondition,data = gait_metrics %>% filter(Metric == "StepLengthSymm"),method = "holm")

ST_Symm_stat <- kruskal.test(Value~DBSCondition,data = gait_metrics %>% filter(Metric == "StepTimeSymm"))
ST_Symm_stat2 <- dunnTest(Value~DBSCondition,data = gait_metrics %>% filter(Metric == "StepTimeSymm"),method = "holm")

SL_CV_stat <- summary_clinic_data %>% filter(Metric == "StepLength") %>% group_by(Side) %>% kruskal_test(cv ~ DBSCondition) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

ST_CV_stat <- summary_clinic_data %>% filter(Metric == "StepTime") %>% group_by(Side) %>% kruskal_test(cv ~ DBSCondition) %>% mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

##### Long-term Rover #####
# Load data
rover_data <- read.csv('/Users/USER/Documents/P2_P3_P4_Long-Term_Rover_Metrics.csv')

data2 <- rover_data %>% 
  mutate(SubjectID = as.factor(SubjectID),
         SubjectID = fct_relevel(SubjectID,c("P2","P3","P4")),
         Condition = as.factor(Condition),
         Condition = fct_relevel(Condition,c("clinical","ramp_up","ramp_down"))) %>% 
  select(-Notes) %>% 
  pivot_longer(!c(SubjectID,Condition,Day),names_to = "metric",values_to = "value")

data3 <- data2 %>% group_by(SubjectID,Condition,metric) %>% 
  summarise(mean = mean(abs(value), na.rm = TRUE), sd = sd(abs(value), na.rm = TRUE), se = sd/sqrt(n())) %>% 
  ungroup() %>% 
  mutate(Condition = fct_recode(Condition, "C" = "clinical", "RU" = "ramp_up", "RD" = "ramp_down"))

data4 <- data2 %>% group_by(Condition,metric) %>% 
  summarise(mean = mean(abs(value), na.rm = TRUE), sd = sd(abs(value), na.rm = TRUE), se = sd/sqrt(n()), cv = sd/mean) %>% 
  ungroup() %>% 
  mutate(Condition = fct_recode(Condition, "C" = "clinical", "RU" = "ramp_up", "RD" = "ramp_down"))

# Stats
rover_KW_test <- data2 %>% group_by(SubjectID, metric) %>% kruskal_test(value ~ Condition) %>% 
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

rover_Dunn_test <- data2 %>% group_by(SubjectID, metric) %>% dunn_test(value ~ Condition) %>% 
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA))

##### Post-motor Diary #####
# Load data
motor_diary_data <- read.csv('/Users/USER/Documents/Post-testing Motor Diary.csv') %>% 
  mutate(Rigidity = as.factor(case_when(str_detect(Rigidity,"Better") ~ "Better",str_detect(Rigidity,"same") ~ "Same",str_detect(Rigidity,"Worse") ~ "Worse")),
         Tremor = as.factor(case_when(str_detect(Tremor,"Better") ~ "Better",str_detect(Tremor,"same") ~ "Same",str_detect(Tremor,"Worse") ~ "Worse")),
         Dyskinesia = as.factor(case_when(str_detect(Dyskinesia,"Better") ~ "Better",str_detect(Dyskinesia,"same") ~ "Same",str_detect(Dyskinesia,"Worse") ~ "Worse")),
         NumFall = factor(NumFall,levels = c("0","1","2-4","5+")),
         NumFreeze = factor(NumFreeze, levels = c("0","1","2-4","5+")))

aggregate_md_data <- motor_diary_data %>% filter(Setting != "") %>% 
  group_by(SubjectID,Setting) %>%
  pivot_longer(cols = c(5:7,17,20), names_to = "columns", values_to = "value") %>%
  count(columns, value) %>% 
  mutate(columns = as.factor(columns),
         Setting = factor(case_when(str_detect(Setting,"Clinical") ~ "cDBS",str_detect(Setting,"ramp-up") ~ "RU-aDBS",str_detect(Setting,"ramp-down") ~ "RD-aDBS"),levels = c("cDBS","RU-aDBS","RD-aDBS")))

filt_motor_diary_data <-motor_diary_data %>% filter(Setting != "") %>% 
  mutate(Setting = factor(case_when(str_detect(Setting,"Clinical") ~ "cDBS",str_detect(Setting,"ramp-up") ~ "RU-aDBS",str_detect(Setting,"ramp-down") ~ "RD-aDBS"),levels = c("cDBS","RU-aDBS","RD-aDBS")))

# Stats
rigidity_md_model <- clm(Rigidity ~ Setting, data = filt_motor_diary_data %>% filter(SubjectID %in% c("P2","P3","P4")) %>% mutate(Rigidity = factor(Rigidity,levels = c("Worse","Same","Better"))))
dyskinesia_md_model <- clm(Dyskinesia ~ Setting, data = filt_motor_diary_data %>% filter(SubjectID %in% c("P2","P3","P4")) %>% mutate(Dyskinesia = factor(Dyskinesia,levels = c("Worse","Same","Better"))))
tremor_md_model <- clm(Tremor ~ Setting, data = filt_motor_diary_data %>% filter(SubjectID %in% c("P2","P3","P4")) %>% mutate(Tremor = factor(Tremor,levels = c("Worse","Same","Better"))))
fall_md_model <- clm(NumFall ~ Setting, data = filt_motor_diary_data %>% filter(SubjectID %in% c("P2","P3","P4")) %>% mutate(NumFall = factor(NumFall,levels = c("5+","2-4","1","0"))))
freeze_md_model <- clm(NumFreeze ~ Setting, data = filt_motor_diary_data %>% filter(SubjectID %in% c("P2","P3","P4")) %>% mutate(NumFreeze = factor(NumFreeze,levels = c("5+","2-4","1","0"))))
