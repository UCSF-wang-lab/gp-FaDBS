library(tidyverse)
library(cowplot)
library(export)
library(gghalves)
library(ggsignif)

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
clinic_data <-read.csv('Figure5D_ED2.csv')
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

# Plots
clinic_sl_plot <- ggplot()+
  geom_half_violin(data = gait_metrics %>% filter(Metric == "StepLength"),
                   aes(x = DBSCondition,y = Value, fill = Side, split = Side),
                   linewidth = 0.1,position = "identity")+
  geom_signif(data = data.frame(SubjectID = c("P2","P2","P3","P4"),
                                start = c("cDBS","RU-aDBS","RU-aDBS","cDBS"),
                                end = c("RU-aDBS","RD-aDBS","RD-aDBS","RD-aDBS"),
                                y = c(1.02,0.92,0.92,0.97),
                                label = c("***","*","***","***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#fb8072", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_signif(data = data.frame(SubjectID = c("P3"),
                                start = c("cDBS"),
                                end = c("RU-aDBS"),
                                y = c(1.02),
                                label = c("***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#80b1d3", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_signif(data = data.frame(SubjectID = c("P3","P4","P4"),
                                start = c("cDBS","cDBS","RU-aDBS"),
                                end = c("RD-aDBS","RU-aDBS","RD-aDBS"),
                                y = c(0.97,1.02,0.92),
                                label = c("* / ***","*** / ***","*** / ***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#8dd3c7", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  scale_x_discrete()+
  scale_fill_manual(name = "Leg:",
                    labels = c("Left","Right"),
                    values = c("#fb8072","#80b1d3"))+
  facet_grid(~SubjectID,labeller = labeller(SubjectID = c("P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4"))) + 
  ylab("meters") + 
  ggtitle("Step Length") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,2)),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

clinic_st_plot <- ggplot()+
  geom_half_violin(data = gait_metrics %>% filter(Metric == "StepTime"),
                   aes(x = DBSCondition,y = Value, fill = Side, split = Side),
                   linewidth = 0.1, position = "identity")+
  geom_signif(data = data.frame(SubjectID = c("P4"),
                                start = c("RU-aDBS"),
                                end = c("RD-aDBS"),
                                y = c(0.75),
                                label = c("***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#fb8072", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_signif(data = data.frame(SubjectID = c("P4"),
                                start = c("cDBS"),
                                end = c("RU-aDBS"),
                                y = c(0.81),
                                label = c("***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#80b1d3", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  geom_signif(data = data.frame(SubjectID = c("P2","P2","P2","P3","P3","P3","P4"),
                                start = c("cDBS","cDBS","RU-aDBS","cDBS","cDBS","RU-aDBS","cDBS"),
                                end = c("RU-aDBS","RD-aDBS","RD-aDBS","RU-aDBS","RD-aDBS","RD-aDBS","RD-aDBS"),
                                y = c(0.81,0.78,0.75,0.81,0.78,0.75,0.78),
                                label = c("** / *","* / *","*** / ***","* / ***","*** / ***","*** / ***","** / ***")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#8dd3c7", size = 0.25, textsize = 1.5, vjust = 0.5,tip_length = 0, manual = TRUE)+
  scale_x_discrete()+
  scale_fill_manual(name = "Leg:",
                    labels = c("Left","Right"),
                    values = c("#fb8072","#80b1d3"))+
  facet_grid(~SubjectID,labeller = labeller(SubjectID = c("P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4"))) + 
  ylab("seconds") + 
  ggtitle("Step Time") + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,2)),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_sl_CV_L_plot = ggplot(data = summary_clinic_data %>% filter(Metric == "StepLength", Side == "Left"), aes(x = DBSCondition, y = cv, group = SubjectID))+
  geom_point(data = summary_clinic_data %>% filter(Metric == "StepLength", Side == "Left"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_dodge2(width = 0.25), size = 1) + 
  geom_line(aes(color = SubjectID),position = position_dodge2(width = 0.25), linewidth = 0.25, show.legend = FALSE)+
  scale_color_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+
  ggtitle("Left") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_sl_CV_R_plot = ggplot(data = summary_clinic_data %>% filter(Metric == "StepLength", Side == "Right"), aes(x = DBSCondition, y = cv, group = SubjectID))+
  geom_point(data = summary_clinic_data %>% filter(Metric == "StepLength", Side == "Right"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_dodge2(width = 0.25), size = 1) + 
  geom_line(aes(color = SubjectID),position = position_dodge2(width = 0.25), linewidth = 0.25)+
  scale_color_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+
  ggtitle("Right") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_st_CV_L_plot = ggplot(data = summary_clinic_data %>% filter(Metric == "StepTime", Side == "Left"), aes(x = DBSCondition, y = cv, group = SubjectID))+
  geom_point(data = summary_clinic_data %>% filter(Metric == "StepTime", Side == "Left"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_dodge2(width = 0.25), size = 1) + 
  geom_line(aes(color = SubjectID),position = position_dodge2(width = 0.25), linewidth = 0.25)+
  scale_color_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+
  ggtitle("Left") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_st_CV_R_plot = ggplot(data = summary_clinic_data %>% filter(Metric == "StepTime", Side == "Right"), aes(x = DBSCondition, y = cv, group = SubjectID))+
  geom_point(data = summary_clinic_data %>% filter(Metric == "StepTime", Side == "Right"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_dodge2(width = 0.25), size = 1) + 
  geom_line(aes(color = SubjectID),position = position_dodge2(width = 0.25), linewidth = 0.25)+
  scale_color_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Coefficient of Variation")+
  ggtitle("Right") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_step_length_symmetry_plot = ggplot(data = summary_clinic_data %>% filter(Metric == "StepLengthSymm"), aes(x = DBSCondition, y = mean, group = SubjectID))+
  geom_point(data = summary_clinic_data %>% filter(Metric == "StepLengthSymm"),aes(color = SubjectID,shape = SubjectID, fill = SubjectID), position = position_dodge2(width = 0.25),size = 1) +
  geom_line(aes(color = SubjectID),position = position_dodge2(width = 0.25), linewidth = 0.25, show.legend = FALSE)+
  geom_errorbar(aes(ymin = lower,ymax = upper,color = SubjectID),width = 0.25,position = position_dodge2(width=0.25), linewidth = 0.25, show.legend = FALSE) + 
  geom_signif(data = data.frame(SubjectID = c("P2","P2","P2"),
                                start = c("cDBS","cDBS","RU-aDBS"),
                                end = c("RD-aDBS","RU-aDBS","RD-aDBS"),
                                y = c(0.18,0.17,0.16),
                                label = c("")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "black", size = 0.25, textsize = 2, vjust = 0.5,tip_length = 0, manual = TRUE)+
  scale_color_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Absolute Asymmetry")+
  ggtitle("Step Length Symmetry") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

group_step_time_symmetry_plot = ggplot(data = summary_clinic_data %>% filter(Metric == "StepTimeSymm"), aes(x = DBSCondition, y = mean, group = SubjectID))+
  geom_point(data = summary_clinic_data %>% filter(Metric == "StepTimeSymm"),aes(color = SubjectID, shape = SubjectID, fill = SubjectID), position = position_dodge2(width = 0.25),size = 1) +
  geom_line(aes(color = SubjectID),position = position_dodge2(width = 0.25), linewidth = 0.25, show.legend = FALSE)+
  geom_errorbar(aes(ymin = lower,ymax = upper,color = SubjectID),width = 0.25,position = position_dodge2(width=0.25), linewidth = 0.25, show.legend = FALSE) + 
  geom_signif(data = data.frame(SubjectID = c("P2","P2","P2"),
                                start = c("cDBS","cDBS","RU-aDBS"),
                                end = c("RD-aDBS","RU-aDBS","RD-aDBS"),
                                y = c(0.09,0.085,0.08),
                                label = c("")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "black", size = 0.25, textsize = 2, vjust = 0.5,tip_length = 0, manual = TRUE)+
  scale_color_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Absolute Asymmetry")+
  ggtitle("Step Time Symmetry") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

signif_color_plot <- ggplot(data = data.frame(x = factor(c("Lunch","Dinner","Lunch","Dinner","Lunch","Dinner"), levels=c("Lunch","Dinner")),
                                              y = c(1,1,2,2,3,3),
                                              signif = factor(c("Left Leg Only","Right Leg Only","Both Legs"), levels=c("Left Leg Only","Right Leg Only","Both Legs"))),
                            aes(x = x,y = y,color = signif))+
  scale_color_manual(name = "Significance:",
                     values = c("#fb8072","#80b1d3","#8dd3c7"))+
  geom_line()+
  theme_bw(base_size = 5) + 
  theme(legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),)

##### Combine plots #####
common_clinic_fill_legend <- get_plot_component(clinic_sl_plot + 
                                                  guides(color = guide_legend(nrow = 1,ncol = 2)) +
                                                  theme(legend.position = "bottom",
                                                        legend.text = element_text(margin = margin(0,5,0,5))),
                                                'guide-box-bottom',return_all = TRUE)

common_clinic_color_legend <- get_plot_component(signif_color_plot + 
                                                   guides(color = guide_legend(nrow = 1,ncol = 3)) +
                                                   theme(legend.position = "bottom",
                                                         legend.text = element_text(margin = margin(0,5,0,5))),
                                                 'guide-box-bottom',return_all = TRUE)

common_clinic_symm_legend <- get_plot_component(group_step_length_symmetry_plot + 
                                                  guides(color = guide_legend(nrow = 1,ncol = 3)) +
                                                  theme(legend.position = "bottom",
                                                        legend.text = element_text(margin = margin(0,5,0,5))),
                                                'guide-box-bottom',return_all = TRUE)

common_cv_legend <- get_plot_component(group_sl_CV_L_plot + 
                                         guides(color = guide_legend(nrow = 1,ncol = 3)) +
                                         theme(legend.position = "bottom",
                                               legend.text = element_text(margin = margin(0,5,0,5))),
                                       'guide-box-bottom',return_all = TRUE)

combined_legends <- plot_grid(NULL,common_cv_legend,NULL,common_clinic_fill_legend,NULL,common_clinic_color_legend,NULL,
                              ncol = 7,
                              rel_widths = c(0.5,1,0.5,0.75,0.15,0.75,0.35))

step_length_cv_plot <- plot_grid(group_sl_CV_L_plot + theme(legend.position = "none"),
                                 NULL,
                                 group_sl_CV_R_plot + theme(legend.position = "none"),
                                 ncol = 3,
                                 nrow = 1,
                                 rel_widths = c(1,0.05,1))

step_length_cv_title <- ggplot() + ggtitle("Step Length Variability") + theme(plot.title = element_text(hjust = 0.5,size = 8))

step_length_cv_w_title_plot <- plot_grid(step_length_cv_title,step_length_cv_plot,
                                         nrow = 2,
                                         rel_heights = c(0.1,1))

step_time_cv_plot <- plot_grid(group_st_CV_L_plot + theme(legend.position = "none"),
                               NULL,
                               group_st_CV_R_plot + theme(legend.position = "none"),
                               ncol = 3,
                               nrow = 1,
                               rel_widths = c(1,0.05,1))

step_time_cv_title <- ggplot() + ggtitle("Step Time Variability") + theme(plot.title = element_text(hjust = 0.5,size = 8))

step_time_cv_w_title_plot <- plot_grid(step_time_cv_title,step_time_cv_plot,
                                       nrow = 2,
                                       rel_heights = c(0.1,1))

combined_clinic_metrics <- plot_grid(group_step_length_symmetry_plot + theme(legend.position = "none"),NULL,clinic_sl_plot + theme(legend.position = "none"),
                                     NULL,NULL,NULL,
                                     group_step_time_symmetry_plot + theme(legend.position = "none"),NULL,clinic_st_plot + theme(legend.position = "none"),
                                     ncol = 3,
                                     nrow = 3,
                                     rel_widths = c(0.5,0.05,1),
                                     rel_heights = c(1,0.05,1))

combined_clinic_legend <- plot_grid(NULL,common_clinic_symm_legend,NULL,common_clinic_fill_legend,NULL,common_clinic_color_legend,NULL,
                                    ncol = 7,
                                    rel_widths = c(0.25,1,0.25,1,0.01,1,0.4))

combined_clinic_metrics_w_legend <- plot_grid(combined_clinic_metrics,
                                              combined_clinic_legend,
                                              nrow = 2,
                                              rel_heights = c(1,0.05))

combined_clinic_var_metrics <- plot_grid(NULL,NULL,NULL,NULL,NULL,
                                         step_length_cv_w_title_plot,NULL,group_step_length_symmetry_plot + theme(legend.position = "none"),NULL,clinic_sl_plot + theme(legend.position = "none"),
                                         NULL,NULL,NULL,NULL,NULL,
                                         step_time_cv_w_title_plot,NULL,group_step_time_symmetry_plot + theme(legend.position = "none"),NULL,clinic_st_plot + theme(legend.position = "none"),
                                         ncol = 5,
                                         nrow = 4,
                                         rel_widths = c(0.75,0.01,0.5,0.01,1.33),
                                         rel_heights = c(0.05,1,0.05,1))

combined_clinic_var_metrics_w_legend <- plot_grid(combined_clinic_var_metrics,
                                                  NULL,
                                                  combined_legends,
                                                  nrow = 3,
                                                  rel_heights = c(1,0.01,0.05))

combined_clinic_var_metrics <- plot_grid(NULL,
                                         step_length_cv_w_title_plot,
                                         NULL,
                                         step_time_cv_w_title_plot,
                                         nrow = 4,
                                         rel_heights = c(0.05,1,0.05,1))

combined_clinic_var_title <- ggplot() + ggtitle("In-clinic Gait Metrics") + theme(plot.title = element_text(hjust = 0.5,size = 8))

combined_clinic_var_metrics_with_title <- plot_grid(combined_clinic_var_title,combined_clinic_var_metrics,
                                                    nrow = 2,
                                                    rel_heights = c(0.05,1))

##### Long-term Rover #####
# Load data
rover_data <- read.csv('Figure5C.csv')

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

# Plots
rover_sl_plot <- ggplot(data = data3 %>% filter(metric == "Stride_Length"),aes(x = Condition, y = mean, color = SubjectID, shape = SubjectID, fill = SubjectID, group = SubjectID))+
  geom_point(position = position_dodge2(width = 0.25), size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),position = position_dodge2(width = 0.25),width = 0.25, linewidth = 0.25, show.legend = FALSE) +
  geom_path(position = position_dodge2(width = 0.25), show.legend = FALSE, linewidth = 0.25) +
  scale_color_manual(name = "",
                     labels = c("Patient 2","Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_x_discrete(labels = c("C" = "cDBS", "RU" = "RU-aDBS", "RD" = "RD-aDBS")) + 
  theme_bw(base_size = 5) +
  xlab("") + 
  ylab("meters") + 
  ggtitle("Stride Length") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

rover_ws_plot <- ggplot(data = data3 %>% filter(metric == "Walking_Speed"),aes(x = Condition, y = mean, color = SubjectID, shape = SubjectID, fill = SubjectID, group = SubjectID))+
  geom_point(position = position_dodge2(width = 0.25), size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),position = position_dodge2(width = 0.25),width = 0.25, linewidth = 0.25, show.legend = FALSE) +
  geom_path(position = position_dodge2(width = 0.25), show.legend = FALSE, linewidth = 0.25) +
  scale_color_manual(name = "",
                     labels = c("Patient 2","Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_x_discrete(labels = c("C" = "cDBS", "RU" = "RU-aDBS", "RD" = "RD-aDBS")) + 
  theme_bw(base_size = 5) +
  xlab("") + 
  ylab("meters/second") + 
  ggtitle("Walking Speed") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

rover_c_plot <- ggplot(data = data3 %>% filter(metric == "Cadence"),aes(x = Condition, y = mean, color = SubjectID, shape = SubjectID, fill = SubjectID, group = SubjectID))+
  geom_point(position = position_dodge2(width = 0.25), size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),position = position_dodge2(width = 0.25),width = 0.25, linewidth = 0.25, show.legend = FALSE) +
  geom_path(position = position_dodge2(width = 0.25), show.legend = FALSE, linewidth = 0.25) +
  geom_signif(data = data.frame(SubjectID = c("P2"),
                                start = c("C"),
                                end = c("RD"),
                                y = c(118),
                                label = c("")),
              aes(y_position = y,xmin = start,xmax = end,annotations = label), color = "#FAA41D", size = 0.25, textsize = 2, vjust = 0.5,tip_length = 0, manual = TRUE)+
  scale_color_manual(name = "",
                     labels = c("Patient 2","Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_x_discrete(labels = c("C" = "cDBS", "RU" = "RU-aDBS", "RD" = "RD-aDBS")) + 
  theme_bw(base_size = 5) +
  xlab("") + 
  ylab("steps/min") + 
  ggtitle("Cadence") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

rover_symm_plot <- ggplot(data = data3 %>% filter(metric == "Symmetry"),aes(x = Condition, y = mean, color = SubjectID, shape = SubjectID, fill = SubjectID, group = SubjectID))+
  geom_point(position = position_dodge2(width = 0.25), size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),position = position_dodge2(width = 0.25),width = 0.25, linewidth = 0.25, show.legend = FALSE) +
  geom_path(position = position_dodge2(width = 0.25), show.legend = FALSE, linewidth = 0.25) +
  scale_color_manual(name = "",
                     labels = c("Patient 2","Patient 3", "Patient 4"),
                     values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_fill_manual(name = "",
                    labels = c("Patient 2", "Patient 3", "Patient 4"),
                    values = c("#FAA41D","#ED2790","#6BBD46"))+
  scale_shape_manual(name = "",
                     labels = c("Patient 2", "Patient 3", "Patient 4"),
                     values = c(22,23,24))+
  scale_x_discrete(labels = c("C" = "cDBS", "RU" = "RU-aDBS", "RD" = "RD-aDBS")) + 
  theme_bw(base_size = 5) +
  xlab("") + 
  ylab("Absolute Asymmetry") +
  ggtitle("Step Length Symmetry") +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 5,margin = margin(0,0,0,0)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

rover_plot_titles <- ggplot() + ggtitle("Home-monitoring Gait Metrics") + theme(plot.title = element_text(hjust = 0.5,size = 8))

common_rover_legend <- get_plot_component(rover_sl_plot + 
                                            guides(color = guide_legend(nrow = 1,ncol = 3)) +
                                            theme(legend.position = "bottom"),
                                          'guide-box-bottom',return_all = TRUE)

combined_rover_metrics <- plot_grid(rover_sl_plot + theme(legend.position = "none"),NULL,rover_symm_plot + theme(legend.position = "none"),
                                    NULL,NULL,NULL,
                                    rover_c_plot + theme(legend.position = "none"),NULL,rover_ws_plot + theme(legend.position = "none"),
                                    ncol = 3,
                                    nrow = 3,
                                    rel_heights = c(1,0.02,1),
                                    rel_widths = c(1,0.02,1))

combined_rover_plots <- plot_grid(combined_rover_metrics,common_rover_legend,
                                  ncol = 1,
                                  nrow = 2,
                                  rel_heights = c(1,0.05))

rover_plots_with_title <- plot_grid(rover_plot_titles,combined_rover_metrics,
                                    nrow = 2,
                                    rel_heights = c(0.05,1))

##### Post-motor Diary #####
# Load data
motor_diary_data <- read.csv('Figure5AB.csv') %>% 
  mutate(Rigidity = as.factor(case_when(str_detect(Rigidity,"Better") ~ "Better",str_detect(Rigidity,"same") ~ "Same",str_detect(Rigidity,"Worse") ~ "Worse")),
         Tremor = as.factor(case_when(str_detect(Tremor,"Better") ~ "Better",str_detect(Tremor,"same") ~ "Same",str_detect(Tremor,"Worse") ~ "Worse")),
         Dyskinesia = as.factor(case_when(str_detect(Dyskinesia,"Better") ~ "Better",str_detect(Dyskinesia,"same") ~ "Same",str_detect(Dyskinesia,"Worse") ~ "Worse")),
         NumFall = factor(NumFall,levels = c("0","1","2-4","5+")),
         NumFreeze = factor(NumFreeze, levels = c("0","1","2-4","5+")))

aggregate_md_data <- motor_diary_data %>% filter(Setting != "") %>% 
  group_by(SubjectID,Setting) %>%
  pivot_longer(cols = c(3:7), names_to = "columns", values_to = "value") %>%
  count(columns, value) %>% 
  mutate(columns = as.factor(columns),
         Setting = factor(case_when(str_detect(Setting,"Clinical") ~ "cDBS",str_detect(Setting,"ramp-up") ~ "RU-aDBS",str_detect(Setting,"ramp-down") ~ "RD-aDBS"),levels = c("cDBS","RU-aDBS","RD-aDBS")))

filt_motor_diary_data <-motor_diary_data %>% filter(Setting != "") %>% 
  mutate(Setting = factor(case_when(str_detect(Setting,"Clinical") ~ "cDBS",str_detect(Setting,"ramp-up") ~ "RU-aDBS",str_detect(Setting,"ramp-down") ~ "RD-aDBS"),levels = c("cDBS","RU-aDBS","RD-aDBS")))

# Plots
cardinal_symptom_plot <- ggplot(aggregate_md_data %>% filter(columns %in% c("Rigidity", "Tremor", "Dyskinesia")),
                                aes(x = Setting, y = n, fill = value)) + 
  geom_bar(position="stack", stat="identity",width = 0.75) + 
  scale_fill_manual(name = "Rating:",
                    # labels = c("Better","Same","Worse"),
                    values = c("#1670B9","#4FA747","#BF2026")) + 
  ylab("Count") + 
  facet_grid(SubjectID ~ columns, switch = "y", labeller = labeller(SubjectID = c("P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4"), columns = c("Dyskinesia" = "Dyskinesia","Rigidity" = "Patient-reported\nStiffness","Tremor" = "Tremor"))) + 
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_text(vjust = 0),
        legend.key.size = unit(0.5,"line"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

fall_freeze_plot <- ggplot(aggregate_md_data %>% filter(columns %in% c("NumFall", "NumFreeze")) %>% mutate(value = factor(value,levels = c("5+","2-4","1","0"))),
                           aes(x = Setting, y = n, fill = value)) + 
  geom_bar(position="stack", stat="identity", width = 0.75) + 
  scale_fill_manual(name = "Amount:",
                    labels = c("5+","2-4", "1","0"),
                    values = c("#BF2026","#F7921E","#4FA747","#1670B9")) + 
  facet_grid(SubjectID~columns, switch = "y", labeller = labeller(columns = c("NumFall" = "Number of Falls","NumFreeze" = "Number of Freezes"),SubjectID = c("P2" = "Patient 2","P3" = "Patient 3","P4" = "Patient 4"))) +
  theme_bw(base_size = 5) + 
  theme(plot.title = element_text(hjust = 0.5,size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_text(vjust = -8),
        legend.key.size = unit(0.5,"line"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 5,vjust=5),
        strip.clip = "off",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"))

combined_md_plots <- plot_grid(fall_freeze_plot,NULL,cardinal_symptom_plot,
                               ncol = 3,
                               rel_widths = c(1,0.01,1.25))

##### Combine all plots #####
clinic_and_rover_gait_metrics <- plot_grid(rover_plots_with_title,NULL,combined_clinic_var_metrics_with_title,
                                           ncol = 3,
                                           rel_widths = c(1,0.05,1))

clinic_and_rover_gait_metrics_w_legend <- plot_grid(clinic_and_rover_gait_metrics,
                                                    common_cv_legend,
                                                    nrow = 2,
                                                    rel_heights = c(1,0.05))

paper_figure <- plot_grid(combined_md_plots,NULL,clinic_and_rover_gait_metrics_w_legend,
                           nrow = 3,
                           rel_heights = c(0.75,0.05,1))