#6.绘制典型代谢物（6个）的对比柱形图
#
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggthemes)
library(gridExtra)
library(ggsignif)
#=============获取典型代谢物在正负样本中的丰度数据=================
met11 <-read.csv("Difference_test/met_new.csv",row.names = 1, check.names = FALSE)
met11 <-t(met11)
met11 <- as.data.frame(met11)  # 将其转换回数据框形式
met11_select <- met11[, c("Sub_Class", "nicotinate","pyridoxamine","pentadecanoate","cholate", "creatinine","docosapentaenoic acid")]
write.csv(met11_select,"Difference_test/met11_select.csv")

#=============绘制差异检验图=================
# 读取数据
data <- read.csv('Difference_test/met11_select.csv', row.names = 1, check.names = FALSE)
colnames(data)[1] <- "Label"

# 检查数据结构
labels <- unique(data$Label)

# 代谢物列表
metabolites <- colnames(data)[-1]
boxplot_list <- list()
letters <- c("A.", "B.", "C.", "D.", "E.", "F.")

# 添加 P 值列表，从significant_metabolites_with_pvalues文件获取
KWP <- c("3.44E-08", "1.59E-05", "4.72E-05", "1.22E-06", "9.10E-05", "7.83E-03")
pairwise_pvalues <- list(
  Control_CD = c(1.64E-09, 6.85E-07, 5.92E-07, 7.56E-08, 5.09E-06, 0.003819272),
  Control_UC = c(0.000214374, 0.009654401, 0.039654913, 0.00156633,0.022917815,0.072438219)
)

# 将显著性标记转换为符号
get_signif_label <- function(p) {
  if (p < 0.001) {
    return("****")
  } else if (p < 0.01) {
    return("***")
  } else if (p < 0.05) {
    return("**")
  } else if (p < 0.1) {
    return("*")
  } else {
    return("")
  }
}

for (i in 1:length(metabolites)) {
  metabolite <- metabolites[i]

  boxplot <- ggplot(data, aes(x = Label, y = .data[[metabolite]])) +
    geom_boxplot(color = "black", fill = NA, width = 0.5) +
    geom_jitter(aes(color = Label), width = 0.2, alpha = 0.5) +
    theme_base() +
    labs(title = paste(letters[i], metabolite),  # 根据当前代谢物名称设置标题
         x = paste("Label (KW p=", KWP[i], ")"), y = "Relative Abundance") +
    theme(legend.position = "none",
          plot.title = element_text(size = 13, hjust = 0.5),  # 设置标题字号为 13
          axis.text = element_text(size = 11),  # 设置轴标签字号为 11
          axis.title = element_text(size = 11)) +
    scale_color_manual(values = c("Control" = "#029f85", "CD" = "#EF9C80", "UC" = "#E34D36")) + 
    coord_cartesian(ylim = c(min(data[[metabolite]]) - 1, max(data[[metabolite]]) + 1)) +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA))
  built_plot <- ggplot_build(boxplot)
  y_range <- built_plot$layout$panel_scales_y[[1]]$range$range
  
  significance_height <- max(y_range) + 0.8 # 例如，加一点偏移

  boxplot <- boxplot +
    annotate("text", x = 1.5, y = significance_height, 
             label = get_signif_label(pairwise_pvalues$Control_CD[i]), size = 5) +
    annotate("text", x = 2.5, y = significance_height-0.5, 
             label = get_signif_label(pairwise_pvalues$Control_UC[i]), size = 5)

  for (label in labels) {
    subset_data <- subset(data, Label == label)
    boxplot <- boxplot + geom_boxplot(data = subset_data, aes(group = Label), color = "black", fill = NA, width = 0.5, position = position_dodge(width = 0.75))
  }
  boxplot <- boxplot + scale_x_discrete(breaks = unique(data$Label))
  boxplot_list[[metabolite]] <- ggplotGrob(boxplot)
}

combined_plot <- grid.arrange(
  grobs = boxplot_list,
  ncol = 3,
  nrow = 2
)

# 保存图像
ggsave("Difference_test/combined_boxplot.png", plot = combined_plot, width = 10, height = 5, units = "in", dpi = 300)
ggsave("Difference_test/combined_boxplot.svg", plot = combined_plot, width = 10, height = 5, units = "in", dpi = 300)