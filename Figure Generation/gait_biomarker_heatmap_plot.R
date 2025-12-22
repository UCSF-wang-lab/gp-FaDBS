library(tidyverse)
library(cowplot)

##### Read in data and filter #####
data <- read.csv('/Users/USER/Documents/search_data.csv')
data_filt <- data %>% filter(Hemisphere == "left",Bitshift == X, FFT_Int == Y)

FreqBins = seq(from = 0, to = 60, by = 500/Y)

data_filt <- data_filt %>% 
  mutate(FreqBandSplit = str_split(FreqBand1,'--', simplify = TRUE)) %>% 
  mutate(Freq1 = FreqBandSplit[,1], Freq2 = FreqBandSplit[,2]) %>% 
  mutate(Freq1 = str_remove(Freq1,"key1_"), Freq2 = str_remove(Freq2,"key1_")) %>% 
  mutate(Freq1 = str_replace(Freq1,"_","."), Freq2 = str_replace(Freq2,"_",".")) %>% 
  mutate(Freq1 = as.numeric(Freq1), Freq2 = as.numeric(Freq2)) %>% 
  select(-FreqBandSplit) %>% 
  mutate(Freq1Bin = findInterval(Freq1, vec = FreqBins), Freq2Bin = findInterval(Freq2, vec = FreqBins))

FreqBin1_max <- data_filt$Freq1Bin[which.max(data_filt$GaitPhaseAccuracy)]
FreqBin2_max <- data_filt$Freq2Bin[which.max(data_filt$GaitPhaseAccuracy)]


##### Plots #####
biomarker_heatmap <- ggplot()+
  geom_tile(data = data_filt,aes(x = Freq1Bin, y = Freq2Bin, fill = GaitPhaseAccuracy)) +
  geom_tile(data = data_filt[data_filt$Freq1Bin == FreqBin1_max & data_filt$Freq2Bin == FreqBin2_max,], aes(x = Freq1Bin, y = Freq2Bin,fill = GaitPhaseAccuracy), color = "red", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "grey75", linewidth = 0.5, linetype = "dashed") + 
  ylab("Frequency Band Stop (Hz)") + 
  xlab("Frequency Band Start (Hz)") + 
  scale_fill_viridis_c(na.value = "white", name = "Predicted Accuracy")+
  scale_x_continuous(breaks = c(4, 8, 17, 22, 30),
                     labels = as.character(round(FreqBins[c(4, 8, 17, 22, 30)],2))) +
  scale_y_continuous(breaks = c(4, 8, 17, 22, 30),
                     labels = as.character(round(FreqBins[c(4, 8, 17, 22, 30)],2))) +
  coord_cartesian(xlim = c(4, 31)) + 
  theme_bw(base_size = 8) + 
  theme(axis.text = element_text(size = 5),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Helvetica"))

legend <- get_plot_component(biomarker_heatmap + 
                          theme(legend.direction = "horizontal",
                                legend.key.height = unit(0.15,"cm"),
                                legend.key.width = unit(0.5,"cm"),
                                legend.text = element_text(size = 4, family = "Helvetica"),
                                legend.title = element_text(size = 5,vjust = 1,hjust = 0.75,margin = margin(r = 0.25, unit = 'cm'),family = "Helvetica")),
                          'guide-box-right',return_all = TRUE)

plots_arranged <- plot_grid(biomarker_heatmap + theme(legend.position = "none"),
                         legend,
                         ncol = 1,
                         nrow = 2,
                         rel_heights = c(1.1,0.1))