library(tidyverse)
library(Cairo)
library(export)
library(cowplot)

# Extract data
accuracy_df = read.csv('/Users/USER/Documents/all_patient_aDBS_accuracy_data.csv')
accuracy_df <- accuracy_df %>% rename(TP = TP..dur.,FN = FN..dur.,FP = FP..dur.,TN = TN..dur.)

P1_data <- accuracy_df %>% filter(SubjectID == "P1") %>% mutate(Side = as_factor(Side))
P2_data <- accuracy_df %>% filter(SubjectID == "P2") %>% mutate(Side = as_factor(Side))
P3_data <- accuracy_df %>% filter(SubjectID == "P3") %>% mutate(Side = as_factor(Side))
P4_data <- accuracy_df %>% filter(SubjectID == "P4") %>% mutate(Side = as_factor(Side))
P5_data <- accuracy_df %>% filter(SubjectID == "P5") %>% mutate(Side = as_factor(Side))

# Create plots
P1_accuracy_plot <- ggplot(data = P1_data, aes(x = Side, y = Accuracy, fill = Side)) +
  geom_col() +
  geom_text(data = P1_data %>% mutate(Accuracy = round(Accuracy,digits = 3)),
            aes(label = paste0(Accuracy*100,"%")), vjust = -0.5, size = 2.5, color = 'black') + 
  geom_segment(aes(x = 0.525,xend = 1.475, y = 0.3495, yend = 0.3495),linetype = 2)+
  geom_segment(aes(x = 1.525,xend = 2.475, y = 0.3472, yend = 0.3472),linetype = 2)+
  scale_fill_manual(name = "Hemisphere:",
                    labels = c("Left","Right"),
                    values = c("#7261AA","#62B3C4"))+
  coord_cartesian(ylim = c(0,0.71)) + 
  scale_y_continuous(labels = scales::percent)+
  ggtitle('Patient 1') + 
  theme_bw(base_size = 8)+
  theme(text = element_text(family = "Helvetica"),
        plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank())

P2_accuracy_plot <- ggplot(data = P2_data, aes(x = Side, y = Accuracy, fill = Side)) +
  geom_col() +
  geom_text(data = P2_data %>% mutate(Accuracy = round(Accuracy,digits = 3)),
            aes(label = paste0(Accuracy*100,"%")), vjust = -0.5, size = 2.5, color = 'black') + 
  geom_segment(aes(x = 0.525,xend = 1.475, y = 0.3893, yend = 0.3893),linetype = 2)+
  geom_segment(aes(x = 1.525,xend = 2.475, y = 0.4194, yend = 0.4194),linetype = 2)+
  scale_fill_manual(name = "Hemisphere:",
                    labels = c("Left","Right"),
                    values = c("#7261AA","#62B3C4"))+
  coord_cartesian(ylim = c(0,0.71)) + 
  scale_y_continuous(labels = scales::percent)+
  ggtitle('Patient 2') + 
  theme_bw(base_size = 8)+
  theme(text = element_text(family = "Helvetica"),
        plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank())

P3_accuracy_plot <- ggplot(data = P3_data, aes(x = Side, y = Accuracy, fill = Side)) +
  geom_col() +
  geom_text(data = P3_data %>% mutate(Accuracy = round(Accuracy,digits = 3)),
            aes(label = paste0(Accuracy*100,"%")), vjust = -0.5, size = 2.5, color = 'black') + 
  geom_segment(aes(x = 0.525,xend = 1.475, y = 0.3382, yend = 0.3382),linetype = 2)+
  geom_segment(aes(x = 1.525,xend = 2.475, y = 0.3568, yend = 0.3568),linetype = 2)+
  scale_fill_manual(name = "Hemisphere:",
                    labels = c("Left","Right"),
                    values = c("#7261AA","#62B3C4"))+
  coord_cartesian(ylim = c(0,0.71)) + 
  scale_y_continuous(labels = scales::percent)+
  ggtitle('Patient 3') + 
  theme_bw(base_size = 8)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank())

P4_accuracy_plot <- ggplot(data = P4_data, aes(x = Side, y = Accuracy, fill = Side)) +
  geom_col() +
  geom_text(data = P4_data %>% mutate(Accuracy = round(Accuracy,digits = 3)),
            aes(label = paste0(Accuracy*100,"%")), vjust = -0.5, size = 2.5, color = 'black') + 
  geom_segment(aes(x = 0.525,xend = 1.475, y = 0.4023, yend = 0.4023),linetype = 2)+
  scale_fill_manual(name = "Hemisphere:",
                    labels = c("Left","Right"),
                    values = c("#7261AA","#62B3C4"))+
  coord_cartesian(ylim = c(0,0.75)) + 
  scale_y_continuous(labels = scales::percent)+
  ggtitle('Patient 4') + 
  theme_bw(base_size = 8)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank())

P5_accuracy_plot <- ggplot(data = P5_data, aes(x = Side, y = Accuracy, fill = Side)) +
  geom_col() +
  geom_text(data = P5_data %>% mutate(Accuracy = round(Accuracy,digits = 3)),
            aes(label = paste0(Accuracy*100,"%")), vjust = -0.5, size = 2.5, color = 'black') + 
  geom_segment(aes(x = 0.525,xend = 1.475, y = 0.3667, yend = 0.3667),linetype = 2)+
  geom_segment(aes(x = 1.525,xend = 2.475, y = 0.4665, yend = 0.4665),linetype = 2)+
  scale_fill_manual(name = "Hemisphere:",
                    labels = c("Left","Right"),
                    values = c("#7261AA","#62B3C4"))+
  coord_cartesian(ylim = c(0,0.71)) + 
  scale_y_continuous(labels = scales::percent)+
  ggtitle('Patient 5') + 
  theme_bw(base_size = 8)+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.border = element_rect(linetype = "solid",color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank())

dashed_line_plot <- ggplot(data = data.frame(x = 1:10,y = (1:10)^2),
                           aes(x = x,y = y,linetype = "Dashed Line"))+
  geom_line(size = 1) +
  scale_linetype_manual(values = c("Dashed Line" = "dashed"),
                        labels = c("")) +
  labs(linetype = "Percent of gait cycle time\nin contralateral swing:") +
  theme_bw(base_size = 8) + 
  theme(legend.key.height = unit(1,"cm"),
        legend.key.width = unit(2.5,"cm"),
        legend.text = element_text(size = 5),
        legend.box.background = element_rect(colour = "black"))

# Combine plots
common_legend <- get_plot_component(dashed_line_plot + 
                                      guides(linetype = guide_legend(nrow = 1,ncol = 1)),
                                    'guide-box-right',return_all = TRUE)

combined_accuracy_grid_plot <- plot_grid(P1_accuracy_plot + theme(legend.position = "none"), NULL,
                                         P2_accuracy_plot + theme(legend.position = "none"),
                                         P3_accuracy_plot + theme(legend.position = "none"), NULL,
                                         P4_accuracy_plot + theme(legend.position = "none"),
                                         P5_accuracy_plot + theme(legend.position = "none"), NULL,
                                         NULL,
                                         ncol = 3,
                                         nrow = 3,
                                         rel_widths = c(1,0.01,1))
