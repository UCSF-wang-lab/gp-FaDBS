library(tidyverse)
library(ggsignif)
library(cowplot)


#### Helper Function ####
outlier_threshold <- function(df,variableName,quantileVal)
{
  iqr_val <- IQR(df[[variableName]],na.rm = TRUE)
  quantile_val <- quantile(df[[variableName]],probs = quantileVal,na.rm = TRUE)
  ifelse(quantileVal>.50, return(quantile_val+1.5*iqr_val), return(quantile_val-1.5*iqr_val))
}

gen_plot <- function(df,side)
{
  out_plot <- ggplot(data = df %>% filter(Side == side)) + 
    geom_density(alpha = 0.5, aes(x = PowerBand, fill = VisitName, after_stat(scaled))) +
    geom_boxplot(data = df %>% filter(VisitName == "Pre optimized DBS"), aes(x = PowerBand, y = -0.30), width = 0.50, outlier.size = 0.25, color = "#8dd3c7" ,inherit.aes = FALSE) +
    geom_boxplot(data = df %>% filter(VisitName == "Post optimized DBS"), aes(x = PowerBand, y = -1.0), width = 0.50, outlier.size = 0.25, color = "#fdb462", inherit.aes = FALSE) + 
    geom_boxplot(data = df %>% filter(VisitName == "aDBS"), aes(x = PowerBand, y = -1.7), width = 0.50, outlier.size = 0.25, color = "#bebada", inherit.aes = FALSE) + 
    scale_fill_manual(labels = c("Pre optimized DBS","Post optimized DBS","aDBS"),
                      values = c("#8dd3c7","#fdb462","#bebada"))+
    scale_y_continuous(breaks = c(-1.7,-1.0,-0.3,0.5),labels = c("","","","Scaled Density")) + 
    xlab("Biomarker Power") + 
    ylab("") + 
    theme_bw(base_size = 5)+
    theme(plot.title = element_text(hjust = 0.5,size = 8),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          legend.key.size = unit(1,"line"),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 5),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.spacing = unit(0.25, "lines"),
          axis.line.x = element_line(linetype = "solid", colour = "black"),
          axis.line.y = element_line(linetype = "solid", colour = "black"))
}

gen_plot2 <- function(df,side)
{
  out_plot <- ggplot(data = df %>% filter(Side == side)) + 
    geom_density(alpha = 0.5, aes(x = PowerBand, fill = VisitName, after_stat(scaled))) +
    geom_boxplot(data = df %>% filter(VisitName == "Post optimized DBS"), aes(x = PowerBand, y = -1.0), width = 0.50, outlier.size = 0.25, color = "#fdb462", inherit.aes = FALSE) + 
    geom_boxplot(data = df %>% filter(VisitName == "aDBS"), aes(x = PowerBand, y = -1.7), width = 0.50, outlier.size = 0.25, color = "#bebada", inherit.aes = FALSE) + 
    scale_fill_manual(labels = c("Post optimized DBS","aDBS"),
                      values = c("#fdb462","#bebada"))+
    scale_y_continuous(breaks = c(-1.7,-1.0,-0.3,0.5),labels = c("","","","Scaled Density")) + 
    xlab("Biomarker Power") + 
    ylab("") + 
    theme_bw(base_size = 5)+
    theme(plot.title = element_text(hjust = 0.5,size = 8),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          legend.key.size = unit(1,"line"),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 5),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.spacing = unit(0.25, "lines"),
          axis.line.x = element_line(linetype = "solid", colour = "black"),
          axis.line.y = element_line(linetype = "solid", colour = "black"))
}

##### Load data ####
data <- read.csv("/Users/USER/Documents/all_patients_biomarker_power.csv")
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

##### Plot #####
P1_left_plot <- gen_plot2(P1,"Left") + ggtitle("Patient 1")
P1_right_plot <- gen_plot2(P1,"Right")
P2_left_plot <- gen_plot(P2,"Left") + ggtitle("Patient 2")
P2_right_plot <- gen_plot(P3,"Right")
P3_left_plot <- gen_plot(P3,"Left") + ggtitle("Patient 3") #+ coord_cartesian(xlim = c(0,600))
P3_right_plot <- gen_plot(P4,"Right") #+ coord_cartesian(xlim = c(0,600))
P4_left_plot <- gen_plot(P4,"Left") + ggtitle("Patient 4") #+ coord_cartesian(xlim = c(0,400))
P5_left_plot <- gen_plot2(P5,"Left") + ggtitle("Patient 5")
P5_right_plot <- gen_plot2(P5,"Right")

##### Combine plots #####
common_legend <- get_plot_component(P2_left_plot + 
                                      guides(fill = guide_legend(nrow = 3,ncol = 1)) +
                                      theme(legend.position = "bottom"),
                                    'guide-box-bottom',return_all = TRUE)

P1_plots <- plot_grid(P1_left_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          NULL,
                          P1_right_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          ncol = 3,
                          nrow = 1,
                          rel_widths = c(1,0.01,1))
P1_plots <- P1_plots + theme(panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 0.5),
                                     plot.margin = unit(c(0.05,0,0.05,0.05),"cm"))

P2_plots <- plot_grid(P2_left_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          NULL,
                          P2_right_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          ncol = 3,
                          nrow = 1,
                          rel_widths = c(1,0.01,1))
P2_plots <- P2_plots + theme(panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 0.5),
                                     plot.margin = unit(c(0.05,0,0.05,0.05),"cm"))

P2_title <- ggdraw() + draw_label("Patient 2",size = 7)

P2_plots <- plot_grid(P2_title,P2_plots,
                          ncol = 1,
                          nrow = 2,
                          rel_heights = c(0.1,1))

P3_plots <- plot_grid(P3_left_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          NULL,
                          P3_right_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          ncol = 3,
                          nrow = 1,
                          rel_widths = c(1,0.01,1))
P3_plots <- P3_plots + theme(panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 0.5),
                                     plot.margin = unit(c(0.05,0,0.05,0.05),"cm"))

P3_title <- ggdraw() + draw_label("Patient 3",size = 7)

P3_plots <- plot_grid(P3_title,P3_plots,
                          ncol = 1,
                          nrow = 2,
                          rel_heights = c(0.1,1))

P4_plots <- plot_grid(P4_left_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                          ncol = 1,
                          nrow = 1)
P4_plots <- P4_plots + theme(panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 0.5),
                                     plot.margin = unit(c(0.05,0,0.05,0.05),"cm"))

P4_title <- ggdraw() + draw_label("Patient 4",size = 7)

P4_plots <- plot_grid(P4_title,P4_plots,
                          ncol = 1,
                          nrow = 2,
                          rel_heights = c(0.1,1))

P5_plots <- plot_grid(P5_left_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                         NULL,
                         P5_right_plot + theme(legend.position = "none",plot.margin = unit(c(0.05,0.25,0.05,0.10),"cm")),
                         ncol = 3,
                         nrow = 1,
                         rel_widths = c(1,0.01,1))
P5_plots <- P5_plots + theme(panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 0.5),
                                   plot.margin = unit(c(0.05,0,0.05,0.05),"cm"))

combined_plots <- plot_grid(P1_left_plot + theme(legend.position = "none") + ylab("Left Hemisphere"), NULL, 
                            P2_left_plot + theme(legend.position = "none") + ylab(""), NULL, 
                            P3_left_plot + theme(legend.position = "none") + ylab(""), NULL, 
                            P5_left_plot + theme(legend.position = "none") + ylab(""), NULL, 
                            P4_left_plot + theme(legend.position = "none") + ylab(""),
                            P1_right_plot + theme(legend.position = "none") + ylab ("Right Hemisphere"), NULL, 
                            P2_right_plot + theme(legend.position = "none") + ylab(""), NULL, 
                            P3_right_plot + theme(legend.position = "none") + ylab(""), NULL, 
                            P5_right_plot + theme(legend.position = "none") + ylab(""), NULL, 
                            common_legend,
                            ncol = 9,
                            nrow = 2,
                            rel_widths = c(1,0.01,1,0.01,1,0.01,1,0.01,1))