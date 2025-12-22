library(tidyverse)
library(cowplot)
library(ggsignif)
library(rstatix)
library(ggpattern)

#### Grab and filter data ####
data <- read.csv('/Users/USER/Documents/gait_cycle_spectrogram_data.csv')

data_tile <- data %>% filter(FreqBandName %in% c("Theta","Alpha","Low Beta","High Beta","Low Gamma")) %>% 
  mutate(FreqBandName = as_factor(FreqBandName),FreqBandName = fct_relevel(FreqBandName,c("Theta","Alpha","Low Beta","High Beta","Low Gamma"))) %>% 
  mutate(Key = if_else(str_detect(Key,"Key0|Key1"),"Pallidum",if_else(str_detect(Key,"Key2"),"M1","PM"))) %>% 
  mutate(Key = as_factor(Key)) %>% 
  mutate(Key = fct_relevel(Key,c("Pallidum","M1","PM"))) %>% 
  mutate(Hemisphere = as_factor(Hemisphere))

#### Stats ####
SW_normality_test <- data_tile %>% group_by(SubjectID,FreqBandName) %>% shapiro_test(Power)

WR_test <- data_tile %>% group_by(SubjectID,Hemisphere,Key,FreqBandName) %>% wilcox_test(Power ~ GaitPhase, detailed = TRUE) %>% 
  mutate(signif = if_else(p<0.05,if_else(p<0.001,if_else(p<0.0001,"***","**"),"*"),NA)) %>% 
  mutate(plotPattern = if_else(Hemisphere == "Left",if_else(estimate>0,"s","n"),if_else(estimate>0,"n","s")))

WR_test$plotPattern[WR_test$p>=0.15] = "n"
WR_test <- WR_test %>% mutate(plotPattern = as_factor(plotPattern))

P4_right_blank <- setNames(data.frame(matrix(ncol = ncol(WR_test), nrow = 0)), colnames(WR_test))
P4_right_blank[1,"SubjectID"] <- "P4"
P4_right_blank[1,"Hemisphere"] <- "Right"
P4_right_blank[1,"Key"] <- "Pallidum"
P4_right_blank[1,"FreqBandName"] <- "Theta"

WR_test <- bind_rows(WR_test,P4_right_blank) %>% 
  mutate(Key = as_factor(Key)) %>% 
  mutate(Key = fct_relevel(Key,c("Pallidum","M1","PM"))) %>% 
  mutate(FreqBandName = as_factor(FreqBandName)) %>% 
  mutate(FreqBandName = fct_relevel(FreqBandName,c("Theta","Alpha","Low Beta","High Beta","Low Gamma"))) 

#### Plots ####
heatmap_P1 <- ggplot(WR_test %>% filter(SubjectID == "P1"),aes(x = Key,y = FreqBandName,fill = p, pattern = plotPattern)) + 
  geom_tile_pattern(color = "black",
                    pattern_fill = "black",
                    pattern_angle = 45,
                    pattern_density = 0.1) +
  scale_pattern_manual(values = c(n = "none",s = "stripe"),labels = c("Contralateral Leg Swing","Ipsilateral Leg Swing")) +
  facet_grid(~Hemisphere) + 
  scale_y_discrete(labels = c(expression(theta),expression(alpha),expression(beta),expression(beta),expression(gamma))) +
  binned_scale("fill",
               "foo",
               ggplot2:::pal_binned(scales::manual_pal(c("#FDE725FF","#22A884FF","#414487FF"))),
               breaks = c(0.0499,0.1499,0.2),
               labels = c("<0.05","<0.15",">0.15"),
               limits = c(0,0.2)) + 
  guides(fill = guide_legend(override.aes = list(pattern = "none"),title = "p-value",label.vjust = 0.5),
         pattern = guide_legend(override.aes = list(fill = "white"),title = "")) +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        strip.text.y.right = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(5,"pt"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(10,'pt'))

heatmap_P2 <- ggplot(WR_test %>% filter(SubjectID == "P2"),aes(x = Key,y = FreqBandName,fill = p, pattern = plotPattern)) + 
  geom_tile_pattern(color = "black",
                    pattern_fill = "black",
                    pattern_angle = 45,
                    pattern_density = 0.1) +
  scale_pattern_manual(values = c(n = "none",s = "stripe"),labels = c("Contralateral Leg Swing","Ipsilateral Leg Swing")) +
  facet_grid(~Hemisphere) + 
  scale_y_discrete(labels = c(expression(theta),expression(alpha),expression(beta),expression(beta),expression(gamma))) +
  binned_scale("fill",
               "foo",
               ggplot2:::pal_binned(scales::manual_pal(c("#FDE725FF","#22A884FF","#414487FF"))),
               breaks = c(0.0499,0.1499,0.2),
               labels = c("<0.05","<0.15",">0.15"),
               limits = c(0,0.2)) + 
  guides(fill = guide_legend(override.aes = list(pattern = "none"),title = "p-value",label.vjust = 0.5),
         pattern = guide_legend(override.aes = list(fill = "white"),title = "")) +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        strip.text.y.right = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(5,"pt"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(10,'pt'))

heatmap_P3 <- ggplot(WR_test %>% filter(SubjectID == "P3"),aes(x = Key,y = FreqBandName,fill = p, pattern = plotPattern)) + 
  geom_tile_pattern(color = "black",
                    pattern_fill = "black",
                    pattern_angle = 45,
                    pattern_density = 0.1) +
  scale_pattern_manual(values = c(n = "none",s = "stripe"),labels = c("Contralateral Leg Swing","Ipsilateral Leg Swing")) +
  facet_grid(~Hemisphere) + 
  scale_y_discrete(labels = c(expression(theta),expression(alpha),expression(beta),expression(beta),expression(gamma))) +
  binned_scale("fill",
               "foo",
               ggplot2:::pal_binned(scales::manual_pal(c("#FDE725FF","#22A884FF","#414487FF"))),
               breaks = c(0.0499,0.1499,0.2),
               labels = c("<0.05","<0.15",">0.15"),
               limits = c(0,0.2)) + 
  guides(fill = guide_legend(override.aes = list(pattern = "none"),title = "p-value",label.vjust = 0.5),
         pattern = guide_legend(override.aes = list(fill = "white"),title = "")) +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        strip.text.y.right = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(5,"pt"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(10,'pt'))

heatmap_P4 <- ggplot(WR_test %>% filter(SubjectID == "P4"),aes(x = Key,y = FreqBandName,fill = p, pattern = plotPattern)) + 
  geom_tile_pattern(color = "black",
                    pattern_fill = "black",
                    pattern_angle = 45,
                    pattern_density = 0.1) +
  scale_pattern_manual(values = c(n = "none",s = "stripe"),labels = c("Contralateral Leg Swing","Ipsilateral Leg Swing")) +
  facet_grid(~Hemisphere) + 
  scale_y_discrete(labels = c(expression(theta),expression(alpha),expression(beta),expression(beta),expression(gamma))) +
  binned_scale("fill",
               "foo",
               ggplot2:::pal_binned(scales::manual_pal(c("#FDE725FF","#22A884FF","#414487FF"))),
               breaks = c(0.0499,0.1499,0.2),
               labels = c("<0.05","<0.15",">0.15"),
               limits = c(0,0.2)) + 
  guides(fill = guide_legend(override.aes = list(pattern = "none"),title = "p-value",label.vjust = 0.5),
         pattern = guide_legend(override.aes = list(fill = "white"),title = "")) +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        strip.text.y.right = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(5,"pt"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(10,'pt'))

heatmap_P5 <- ggplot(WR_test %>% filter(SubjectID == "P5") %>% mutate(Key = recode(Key,"M1" = "S1","PM" = "M1")),
                        aes(x = Key,y = FreqBandName,fill = p, pattern = plotPattern)) + 
  geom_tile_pattern(color = "black",
                    pattern_fill = "black",
                    pattern_angle = 45,
                    pattern_density = 0.1) +
  scale_pattern_manual(values = c(n = "none",s = "stripe"),labels = c("Contralateral Leg Swing","Ipsilateral Leg Swing")) +
  facet_grid(~Hemisphere) + 
  scale_y_discrete(labels = c(expression(theta),expression(alpha),expression(beta),expression(beta),expression(gamma))) +
  binned_scale("fill",
               "foo",
               ggplot2:::pal_binned(scales::manual_pal(c("#FDE725FF","#22A884FF","#414487FF"))),
               breaks = c(0.0499,0.1499,0.2),
               labels = c("<0.05","<0.15",">0.15"),
               limits = c(0,0.2)) + 
  guides(fill = guide_legend(override.aes = list(pattern = "none"),title = "p-value",label.vjust = 0.5),
         pattern = guide_legend(override.aes = list(fill = "white"),title = "")) +
  theme_bw(base_size = 5)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 5),
        strip.text.y.right = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(5,"pt"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(10,'pt'))

p_value_legend <- get_plot_component(heatmap_P1 + 
                                       guides(fill = guide_legend(override.aes = list(pattern = "none"),nrow = 1,ncol = 3,title = "p-value", title.position = "top", title.hjust = 0.5),
                                              pattern = "none") +
                                       theme(legend.position = "bottom"),
                                     'guide-box-bottom',return_all = TRUE)

pattern_legend <- get_plot_component(heatmap_P1 + 
                                       guides(fill = "none",
                                              pattern = guide_legend(override.aes = list(fill = "white", color = "black"),nrow = 1, ncol = 2, title = "Direction",title.position = "top", title.hjust = 0.5)) +
                                       theme(legend.position = "bottom",
                                             legend.key.spacing.x = unit(0.5,'cm')),
                                     'guide-box-bottom',return_all = TRUE)

common_legend <- plot_grid(NULL,
                           p_value_legend,
                           NULL,
                           pattern_legend,
                           NULL,
                           ncol = 5,
                           nrow = 1,
                           rel_widths = c(0.6,0.5,0.2,0.5,0.6))

heatmap_combined <- plot_grid(heatmap_P1 + theme(legend.position = "none"), 
                              NULL,
                              heatmap_P3 + theme(legend.position = "none"),
                              NULL,
                              heatmap_P3 + theme(legend.position = "none"),
                              NULL,
                              heatmap_P4 + theme(legend.position = "none"),
                              NULL,
                              heatmap_P5 + theme(legend.position = "none"),
                              ncol = 9,
                              nrow = 1,
                              rel_widths = c(1,0.01,1,0.01,1,0.01,1,0.01,1))

heatmap_combined_with_legend <- plot_grid(heatmap_combined,
                                          NULL,
                                          common_legend,
                                          NULL,
                                          ncol = 1,
                                          nrow = 4,
                                          rel_heights = c(1,0.05,0.05,0.1))