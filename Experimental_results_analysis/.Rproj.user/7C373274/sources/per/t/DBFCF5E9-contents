#7.Task1中关联数量柱状图及微生物-代谢物关联热图的绘制
#
library(tidyverse)
library(pheatmap)
library(ggthemes)
library(ggplot2)
library(ggpubr)

#================================================================================
data <- read.csv("Result_corr/associations_heatmap.csv", header = TRUE)
df <- data[,-c(25)]
bk <- c(seq(-1, 1,by = 1))
rownames(df) <- data$Metabolite

p <- pheatmap(df,
              cluster_cols = T, cluster_rows = T, scale = "none",
              treeheight_col = 20, treeheight_row = 20,
              display_numbers = F,
              legend = F,
              legend_labels = c("Exist"),
              border_color = "black",
              color = colorRampPalette(c("grey","white", "#E34D36"))(length(bk)),
              fontsize = 6.5, 
              fontfamily = "Arial",
              main = "Associations Between Differential Metabolites and Representative Microbes")
ggsave("Result_corr/heatmap_mm.png", plot = p, width = 6, height = 3, dpi = 300)
ggsave("Result_corr/heatmap_mm.svg", plot = p, width = 6, height = 3, dpi = 300)
#================================================================================
data <- read.csv("Result_corr/transposed_data-24.csv", header = TRUE)
df <- data[,-c(25,26)]
bk <- c(seq(-1, 2,by = 1))
rownames(df) <- data$Metabolite
rownames(df)
annotation_row = data.frame(Methods = factor(data$Methods))
rownames(annotation_row) <- data$Metabolite

mycolors <- c("#66C2A5","#FC8D62","#BEB8DC","#4FB9D5","#F6CAE5","#8691AC","#BB9727")
names(mycolors) <- unique(annotation_row$Methods)
mycolors <- list(Methods = mycolors)

p <- pheatmap(df,
               cluster_cols = F, cluster_rows = F, scale = "none",
               treeheight_col = 0, treeheight_row = 0,
               display_numbers = F,
               annotation_row = annotation_row,
               annotation_colors = mycolors,
               legend = F,
               legend_labels = c("Not exist", "Exist"),
               border_color = "black",
               color = colorRampPalette(c("grey","white", "#E34D36","#DFE9F4"))(length(bk)),
               fontsize = 10, 
               fontfamily = "Arial",
               main = "Associations Between 24 Representative Microbes and 6 Differential Metabolites")
ggsave("Result_corr/heatmap.png", plot = p, width = 9, height = 8, dpi = 300)
ggsave("Result_corr/heatmap.svg", plot = p, width = 9, height = 8, dpi = 300)
#=============================================================================
library(ggplot2)
library(tidyr)

# 读取CSV数据
data <- read.csv("Result_corr/duidie.csv", header = TRUE)

# 创建一个命名的颜色向量
colors <- c("CHO" = "#66C2A5",
            "PYR" = "#FC8D62",
            "NIC" = "#BEB8DC",
            "CRE" = "#4FB9D5",
            "DPA" = "#F6CAE5",
            "PNT" = "#8691AC")

# 将数据变为长格式
data_long <- pivot_longer(data, cols = -Methods, names_to = "Metabolite", values_to = "Count")
data_long$Methods <- factor(data_long$Methods, levels = c("Spearman", "O2PLS", "SparseCCA", "MIMOSA2", "BiomeNED", "ENVIM", "mNODE"))

data_filtered <- data_long[data_long$Count > 0, ]

# 绘制堆叠柱状图，并在每个柱形上添加标签
stacked_barplot <-ggplot(data_filtered, aes(x = Methods, y = Count, fill = Metabolite, label = Count)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = colors) +
  labs(title = "",
       x = "",
       y = "Metabolite-Microbe Associations Count") +
  theme_stata()+
  theme(legend.position = "right",
        # legend.text = element_text(size = 10),
        # axis.text.x = element_text(angle = 0, hjust = 1),
        legend.key.size = unit(1, "lines"),
        plot.title = element_text(color = "black", size = 15, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +  # 旋转X轴标签
  geom_text(position = position_stack(vjust = 0.5), size = 4, color = "black")+
  # guides(fill = guide_legend(nrow = 2)) +  
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))
# 保存为PNG
ggsave("Result_corr/stacked_barplot.png", plot = stacked_barplot, width = 12, height = 6, dpi = 300)
# 保存为SVG
ggsave("Result_corr/stacked_barplot.svg", plot = stacked_barplot, width = 12, height = 6, dpi = 300)
