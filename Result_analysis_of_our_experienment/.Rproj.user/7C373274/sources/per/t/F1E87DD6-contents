#9.绘制task2中，五折交叉验证实验中的折线图
#
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggthemes)

file_path <- "training_result_plot.xlsx"
IBD_data <- read_excel(file_path, sheet = "IBD")
GC_data <- read_excel(file_path, sheet = "GC")
ESRD_data <- read_excel(file_path, sheet = "ESRD")

method_colors <- c("NMSE" = "#E34D36", "ALL_SCC" = "#4FB9D5", "RMSE" = "#029f85", 
                   "TOP50_SCC" = "#455882", "dash2" = "#EF9C80", "MAAPE" = "#8691AC", "dash4" = "#95CFC1")

group1_colors <- c(method_colors["TOP50_SCC"], method_colors["ALL_SCC"], method_colors["MAAPE"])

group2_colors <- c(method_colors["NMSE"], method_colors["RMSE"])

IBD_data_long <- gather(IBD_data, key = "Indicator", value = "Value", -Method)
GC_data_long <- gather(GC_data, key = "Indicator", value = "Value", -Method)
ESRD_data_long <- gather(ESRD_data, key = "Indicator", value = "Value", -Method)

IBD_data_long$Method <- factor(IBD_data_long$Method, levels = c("MelonnPan", "ENVIM", "MMINP", "MiMeNet", "BiomeNED", "mNODE", "LOCATE"))
GC_data_long$Method <- factor(GC_data_long$Method, levels = c("MelonnPan", "ENVIM", "MMINP", "MiMeNet", "BiomeNED", "mNODE", "LOCATE"))
ESRD_data_long$Method <- factor(ESRD_data_long$Method, levels = c("MelonnPan", "ENVIM", "MMINP", "MiMeNet", "BiomeNED", "mNODE", "LOCATE"))
#==============================IBD数据集========================================
# 提取第一组指标
group1_indicators <- c("TOP50_SCC", "ALL_SCC", "MAAPE")
group1_data <- filter(IBD_data_long, Indicator %in% group1_indicators)

# 提取第二组指标
group2_indicators <- c("NMSE", "RMSE")
group2_data <- filter(IBD_data_long, Indicator %in% group2_indicators)

coeff1 <- 10

group2_data_scaled <- group2_data %>%
  mutate(Scaled_Value = Value / coeff1)

IBD_plot <- ggplot(data = IBD_data_long,aes(x = Method)) +
  geom_bar(data = group1_data, aes(y = Value, fill = Indicator), stat = "identity", position = "dodge",width = 0.5,alpha=0.8) + 
  geom_line(data = group2_data_scaled, aes(x = Method, y = Scaled_Value, group = Indicator, color=Indicator),size=1) +
  labs(title = "IBD Dataset", x = "Method", y = "Value") +
  scale_y_continuous(name = "Y1 Value(bar)", sec.axis = sec_axis(~.*coeff1, name = "Y2 Value(solid)")) +
  scale_color_manual(values = group2_colors) + # Specify colors for lines
  scale_fill_manual(values = group1_colors) +
  theme_stata() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),  # Adjust the size of the first y-axis title
        axis.title.y.right = element_text(size = 14)) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))
#==============================GC数据集========================================
# 提取第一组指标
group1_indicators <- c("TOP50_SCC", "ALL_SCC", "MAAPE")
group1_data <- filter(GC_data_long, Indicator %in% group1_indicators)

# 提取第二组指标
group2_indicators <- c("NMSE", "RMSE")
group2_data <- filter(GC_data_long, Indicator %in% group2_indicators)

coeff <- 2

# Create a new data frame for the second group of indicators with a separate column for the scaled values
group2_data_scaled <- group2_data %>%
  mutate(Scaled_Value = Value / coeff)

# Create the dual y-axis plot
GC_plot <- ggplot(data = IBD_data_long,aes(x = Method)) +
  geom_bar(data = group1_data, aes(y = Value, fill = Indicator), stat = "identity", position = "dodge",width = 0.5,alpha=0.8) + 
  geom_line(data = group2_data_scaled, aes(x = Method, y = Scaled_Value, group = Indicator, color=Indicator),size=1) +
  labs(title = "GC Dataset", x = "Method", y = "Value") +
  scale_y_continuous(name = "Y1 Value(bar)", sec.axis = sec_axis(~.*coeff, name = "Y2 Value(solid)")) +
  scale_color_manual(values = group2_colors) + # Specify colors for lines
  scale_fill_manual(values = group1_colors) +
  theme_stata() +
  theme(legend.position = "right",
        legend.key.height = unit(1.2, "line"),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),  # Adjust the size of the first y-axis title
        axis.title.y.right = element_text(size = 14),
        legend.spacing = unit(0.8, "lines")) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))+
  guides(color = guide_legend(title = "Line Legend", order = 1),
         fill = guide_legend(title="Bar Legend",order = 2),title.position = "top", title.hjust = 0.5)
#==============================ESRD数据集========================================
# 提取第一组指标
group1_indicators <- c("TOP50_SCC", "ALL_SCC", "MAAPE")
group1_data <- filter(ESRD_data_long, Indicator %in% group1_indicators)

# 提取第二组指标
group2_indicators <- c("NMSE", "RMSE")
group2_data <- filter(ESRD_data_long, Indicator %in% group2_indicators)

coeff <- 2

# Create a new data frame for the second group of indicators with a separate column for the scaled values
group2_data_scaled <- group2_data %>%
  mutate(Scaled_Value = Value / coeff)

# Create the dual y-axis plot
ESRD_plot <- ggplot(data = IBD_data_long,aes(x = Method)) +
  geom_bar(data = group1_data, aes(y = Value, fill = Indicator), stat = "identity", position = "dodge",width = 0.5,alpha=0.8) + 
  geom_line(data = group2_data_scaled, aes(x = Method, y = Scaled_Value, group = Indicator, color=Indicator),size=1) +
  labs(title = "ESRD Dataset", x = "Method", y = "Value") +
  scale_y_continuous(name = "Y1 Value(bar)", sec.axis = sec_axis(~.*coeff, name = "Y2 Value(solid)")) +
  scale_color_manual(values = group2_colors) + # Specify colors for lines
  scale_fill_manual(values = group1_colors) +
  theme_stata() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),  # Adjust the size of the first y-axis title
        axis.title.y.right = element_text(size = 14)) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))

combined_plots <- plot_grid(IBD_plot, GC_plot, ESRD_plot,  nrow = 3)
ggsave("Final_Line_Plots.jpg", combined_plots, width = 14, height = 10, dpi = 300)
ggsave("Final_Line_Plots.svg", combined_plots, width = 14, height = 10, dpi = 300)