library(tidyverse)
library(gghalves)
library(Cairo)
library(export)
library(ggsignif)
library(cowplot)

#### Helper functions ####
outlier_threshold <- function(df,variableName,quantileVal)
{
  iqr_val <- IQR(df[[variableName]],na.rm = TRUE)
  quantile_val <- quantile(df[[variableName]],probs = quantileVal,na.rm = TRUE)
  ifelse(quantileVal>.50, return(quantile_val+1.5*iqr_val), return(quantile_val-1.5*iqr_val))
}

gen_plot <- function(df,variableName,plotType)
{
  if (plotType == "gait"){
    out_plot <- ggplot()+
      geom_half_violin(data = df %>% group_by(SubjectID,VisitName,Side) %>% 
                         filter(!!sym(variableName)>outlier_threshold(.,variableName,0.25),
                                !!sym(variableName)<outlier_threshold(.,variableName,0.75)),
                       aes(x = VisitName,y = !!sym(variableName), fill = Side, split = Side),
                       linewidth = 0.1,position = "identity")+
      scale_x_discrete()+
      scale_fill_manual(name = "Leg:",
                        labels = c("Left","Right"),
                        values = c("#fb8072","#80b1d3"))+
      facet_grid(~SubjectID,labeller = labeller(SubjectID = c("P1" = "Patient 1", "P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4","P5" = "Patient 5"))) + 
      theme_bw(base_size = 5)+
      theme(plot.title = element_text(hjust = 0.5,size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 6),
            axis.text = element_text(size = 6),
            legend.key.size = unit(0.5,"line"),
            strip.background = element_blank(),
            strip.text = element_text(size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.spacing = unit(0.25, "lines"),
            axis.line.x = element_line(linetype = "solid", colour = "black"),
            axis.line.y = element_line(linetype = "solid", colour = "black"))
  } else if (plotType == "var"){
    out_plot <- ggplot(data = summary_df,aes(x = VisitName,y = !!sym(variableName),color = Side,group = Side))+
      geom_point(position = position_dodge2(width = 0.5)) + 
      geom_line(position = position_dodge2(width = 0.5))+
      facet_grid(~SubjectID)+
      scale_x_discrete()+
      scale_color_manual(name = "Leg:",
                         labels = c("Left","Right"),
                         values = c("#fb8072","#80b1d3"))+
      facet_grid(~SubjectID,labeller = labeller(SubjectID = c("P1" = "Patient 1", "P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4","P5" = "Patient 5"))) + 
      theme_bw(base_size = 5)+
      theme(plot.title = element_text(hjust = 0.5,size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 6),
            axis.text = element_text(size = 6),
            legend.key.size = unit(0.5,"line"),
            legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
            strip.background = element_blank(),
            strip.text = element_text(size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.spacing = unit(0.25, "lines"),
            axis.line.x = element_line(linetype = "solid", colour = "black"),
            axis.line.y = element_line(linetype = "solid", colour = "black"))
  } else{
    out_plot <- ggplot(clinicOpt_to_aDBS_w_symmetry %>% 
                         filter(Side == "L") %>% 
                         group_by(SubjectID,VisitName) %>% 
                         filter(!!sym(str_remove(variableName,"Symmetry"))>outlier_threshold(.,str_remove(variableName,"Symmetry"),0.25),
                                !!sym(str_remove(variableName,"Symmetry"))<outlier_threshold(.,str_remove(variableName,"Symmetry"),0.75),
                                !!sym(variableName) < 0.25,
                                !!sym(variableName) > -0.25),
                       aes(x = VisitName,y = !!sym(variableName))) +
      geom_boxplot(outlier.shape = NA,show.legend = FALSE)+
      geom_hline(yintercept = 0,linetype = "dotted")+
      ylab("+ left / - right leg longer")+
      scale_x_discrete()+
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      facet_grid(~SubjectID,labeller = labeller(SubjectID = c("P1" = "Patient 1", "P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4","P5" = "Patient 5"))) + 
      theme_bw(base_size = 5)+
      theme(plot.title = element_text(hjust = 0.5,size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 6),
            axis.text = element_text(size = 6),
            strip.background = element_blank(),
            strip.text = element_text(size = 5.5),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.spacing = unit(0.25, "lines"),
            axis.line.x = element_line(linetype = "solid", colour = "black"),
            axis.line.y = element_line(linetype = "solid", colour = "black"))
  }
  return (out_plot)
}

#### Load in data, filter, and format ####
# P1
P1_gait_data <- read.csv('P1_gait_metrics.csv')

P1_clinicOpt_to_aDBS <- P1_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P2
P2_gait_data <- read.csv('P2_gait_metrics.csv');

P2_clinicOpt_to_aDBS <- P2_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P3
P3_gait_data <- read.csv('P3_gait_metrics.csv')

P3_clinicOpt_to_aDBS <- P3_gait_data %>% filter(VisitName %in% c("dbsOptBilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P4
P4_gait_data <- read.csv('P4_gait_metrics.csv')

P4_clinicOpt_to_aDBS <- P4_gait_data %>% filter(VisitName %in% c("dbsOptUnilateral","aDBS")) %>% mutate(Side = as_factor(Side))

# P5
P5_gait_data <- read.csv('P5_gait_metrics.csv')

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

#### Plots ####
# Step Length
step_length_plot <- gen_plot(clinicOpt_to_aDBS_w_symmetry,"StepLength","gait") +
  geom_signif(data = data.frame(SubjectID = c("P1","P2","P4","P5"),
                                start = c("cDBS","cDBS","cDBS","cDBS"),
                                end = c("aDBS","aDBS","aDBS","aDBS"),
                                y = c(0.90,0.90,0.90,0.90),
                                label = c("***","*","*","***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#fb8072", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_signif(data = data.frame(SubjectID = c("P1","P2","P3","P4","P5"),
                                start = c("cDBS","cDBS","cDBS","cDBS","cDBS"),
                                end = c("aDBS","aDBS","aDBS","aDBS","aDBS"),
                                y = c(0.85,0.85,0.85,0.85,0.85),
                                label = c("***","***","***","***","***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#80b1d3", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  ylab("meters")+
  ggtitle("Step Length")

# Step Time
step_time_plot <- gen_plot(clinicOpt_to_aDBS_w_symmetry,"StepTime","gait") + 
  geom_signif(data = data.frame(SubjectID = c("P1","P2","P3","P4","P5"),
                                start = c("cDBS","cDBS","cDBS","cDBS","cDBS"),
                                end = c("aDBS","aDBS","aDBS","aDBS","aDBS"),
                                y = c(0.85,0.85,0.85,0.85,0.85),
                                label = c("***","***","***","***","***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#fb8072", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_signif(data = data.frame(SubjectID = c("P1","P3","P4","P5"),
                                start = c("cDBS","cDBS","cDBS","cDBS"),
                                end = c("aDBS","aDBS","aDBS","aDBS"),
                                y = c(0.80,0.80,0.80,0.80),
                                label = c("*","***","***","***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#80b1d3", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  ylab("seconds")+
  ggtitle("Step Time")

# Step Length Symmetry
step_length_symmetry_plot <- gen_plot(clinicOpt_to_aDBS_w_symmetry,"StepLengthSymmetry","sym") + 
  geom_signif(data = data.frame(SubjectID = c("P2","P3","P4","P5"),
                                start = c("cDBS","cDBS","cDBS","cDBS"),
                                end = c("aDBS","aDBS","aDBS","aDBS"),
                                y = c(0.29,0.29,0.29,0.29),
                                label = c("***","***","***","***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  coord_cartesian(ylim = c(-0.25,0.31))+
  ggtitle("Step Length Symmetry")

# Step time symmetry
step_time_symmetry_plot <- gen_plot(clinicOpt_to_aDBS_w_symmetry,"StepTimeSymmetry","sym") + 
  geom_signif(data = data.frame(SubjectID = c("P1","P2","P4","P5"),
                                start = c("cDBS","cDBS","cDBS","cDBS"),
                                end = c("aDBS","aDBS","aDBS","aDBS"),
                                y = c(0.29,0.29,0.29,0.29),
                                label = c("***","***","***","**")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  ggtitle("Step Time Symmetry")

# Grouped symmetry plots
group_step_length_symmetry_plot = ggplot(data = group_step_length, aes(x = VisitName,y = StepLengthSymmetry)) + 
  geom_boxplot(outlier.shape = NaN) + 
  geom_point(data = group_step_length,aes(color = SubjectID,shape = SubjectID,fill = SubjectID), position = position_jitter(width = 0.25, height = 0), size = 0.5, alpha = 0.5) + 
  geom_signif(data = data.frame(start = c("cDBS"),
                                end = c("aDBS"),
                                y = c(0.26),
                                label = c("***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_fill_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c(21,22,23,24,25))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Absolute Asymmetry")+
  ggtitle("Grouped Step\nLength Symmetry ") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_step_time_symmetry_plot = ggplot(data = group_step_time, aes(x = VisitName,y = StepTimeSymmetry)) + 
  geom_boxplot(outlier.shape = NaN) + 
  geom_point(data = group_step_time,aes(color = SubjectID,shape = SubjectID,fill = SubjectID), position = position_jitter(width = 0.25, height = 0), size = 0.5, alpha = 0.5) + 
  geom_signif(data = data.frame(start = c("cDBS"),
                                end = c("aDBS"),
                                y = c(0.26),
                                label = c("*")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                    values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c(21,22,23,24,25))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Absolute Asymmetry")+
  ggtitle("Grouped Step\nTime Symmetry") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_step_length_CV_L_plot = ggplot(data = group_step_CV %>% filter(Side == "L"), aes(x = VisitName,y = StepLengthCV)) + 
  geom_boxplot(outlier.shape = NaN) + 
  geom_point(data = group_step_CV %>% filter(Side == "L"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_jitter(width = 0.25, height = 0), size = 0.5, alpha = 0.5) + 
  scale_color_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                    values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c(21,22,23,24,25))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+
  ggtitle("Left") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))
  
group_step_length_CV_R_plot = ggplot(data = group_step_CV %>% filter(Side == "R"), aes(x = VisitName,y = StepLengthCV)) + 
  geom_boxplot(outlier.shape = NaN) + 
  geom_point(data = group_step_CV %>% filter(Side == "R"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_jitter(width = 0.25, height = 0), size = 0.5, alpha = 0.5) + 
  scale_color_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                    values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c(21,22,23,24,25))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+
  ggtitle("Right") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_step_time_CV_L_plot = ggplot(data = group_step_CV %>% filter(Side == "L"), aes(x = VisitName,y = StepTimeCV)) + 
  geom_boxplot(outlier.shape = NaN) + 
  geom_point(data = group_step_CV %>% filter(Side == "L"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_jitter(width = 0.25, height = 0), size = 0.5, alpha = 0.5) + 
  geom_signif(data = data.frame(start = c("cDBS"),
                                end = c("aDBS"),
                                y = c(0.13),
                                label = c(""),
                                vjust = 1),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  scale_color_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                    values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c(21,22,23,24,25))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+ 
  ggtitle("Left") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_step_time_CV_R_plot = ggplot(data = group_step_CV %>% filter(Side == "R"), aes(x = VisitName,y = StepTimeCV)) + 
  geom_boxplot(outlier.shape = NaN) + 
  geom_point(data = group_step_CV %>% filter(Side=="R"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_jitter(width = 0.25, height = 0), size = 0.5, alpha = 0.5) + 
  scale_color_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                    values = c("#3CC1C8","#FAA41D","#ED2790","#6BBD46","#3B54A4"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 1","Patient 2", "Patient 3", "Patient 4", "Patient 5"),
                     values = c(21,22,23,24,25))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+ 
  ggtitle("Right") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

##### Figure #####
common_legend <- get_plot_component(step_length_plot + 
                                      guides(color = guide_legend(nrow = 1,ncol = 2)) +
                                      theme(legend.position = "bottom"),
                                    'guide-box-bottom',return_all = TRUE)

common_legend2 <- get_plot_component(group_step_length_symmetry_plot + 
                                       guides(color = guide_legend(keywidth = 0.5,nrow = 1,ncol = 5)) +
                                       theme(legend.position = "bottom", 
                                             legend.text = element_text(size = 5,margin = margin(0,0,0,0))),
                                     'guide-box-bottom',return_all = TRUE)

combined_legend <- plot_grid(common_legend2,NULL,common_legend,NULL,
                             ncol = 4,
                             rel_widths = c(1,0.4,1,0.35))

step_length_cv_plot <- plot_grid(group_step_length_CV_L_plot + theme(legend.position = "none"),
                                 NULL,
                                 group_step_length_CV_R_plot + theme(legend.position = "none"),
                                 ncol = 3,
                                 nrow = 1,
                                 rel_widths = c(1,0.01,1))

step_length_cv_title <- ggplot() + ggtitle("Grouped Step Length Variability") + theme(plot.title = element_text(hjust = 0.5,size = 8))

step_length_cv_w_title_plot <- plot_grid(step_length_cv_title,step_length_cv_plot,
                                         nrow = 2,
                                         rel_heights = c(0.1,1))

combined_step_length_metrics_CV_plot <- plot_grid(step_length_cv_w_title_plot,NULL,
                                                  step_length_plot + theme(legend.position = "none"),
                                                  ncol = 3,
                                                  rel_widths = c(0.5,0.01,1))

step_time_cv_plot <- plot_grid(group_step_time_CV_L_plot + theme(legend.position = "none"),
                               NULL,
                               group_step_time_CV_R_plot + theme(legend.position = "none"),
                               ncol = 3,
                               nrow = 1,
                               rel_widths = c(1,0.01,1))

step_time_cv_title <- ggplot() + ggtitle("Grouped Step Time Variability") + theme(plot.title = element_text(hjust = 0.5,size = 8))

step_time_cv_w_title_plot <- plot_grid(step_time_cv_title,step_time_cv_plot,
                                       nrow = 2,
                                       rel_heights = c(0.1,1))

combined_step_time_metrics_CV_plot <- plot_grid(step_time_cv_w_title_plot,NULL,
                                                step_time_plot + theme(legend.position = "none"),
                                                ncol = 3,
                                                rel_widths = c(0.5,0.01,1))

grouped_step_symmetry_grid_plot <- plot_grid(group_step_length_symmetry_plot + theme(legend.position = "none"),
                                             NULL,
                                             group_step_time_symmetry_plot + theme(legend.position = "none"),
                                             ncol = 3,
                                             rel_widths = c(1,0.01,1))



gait_metric_symmetry_grid_plot <- plot_grid(grouped_step_symmetry_grid_plot,
                                               NULL,
                                               step_length_symmetry_plot + theme(axis.text.x = element_text(angle = 45,hjust = 1)),
                                               NULL,
                                               step_time_symmetry_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
                                               ncol = 5,
                                               rel_widths = c(1,0.01,1,0.01,1))

paper_figure <- plot_grid(gait_metric_symmetry_grid_plot,NULL,
                          combined_step_length_metrics_CV_plot,NULL,
                          combined_step_time_metrics_CV_plot,NULL,
                          combined_legend,
                          NULL,
                          nrow = 8,
                          rel_heights = c(1,0.1,1,0.1,1,0.1,0.05,0.1))