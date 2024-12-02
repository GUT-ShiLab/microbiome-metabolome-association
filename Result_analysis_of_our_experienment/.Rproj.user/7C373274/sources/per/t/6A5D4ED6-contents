#4.计算微生物-代谢物spearman相关性、Sparse CCA、O2PLS微生物/代谢物载荷
library(ggplot2)
library(reshape2)
metab_clr <- read.csv("correlation/IBD2019_ann/metabolome_clr.csv", row.names = 1, check.names = FALSE)
micro_clr <- read.csv("correlation/IBD2019_ann/microbiome_clr.csv", row.names = 1)

#----------------------------计算Spearman相关性-----------------------------------
# 计算相关系数
correlation <- cor(t(metab_clr), t(micro_clr), method = "spearman")
p_values <- matrix(NA, ncol = ncol(correlation), nrow = nrow(correlation))
# 计算每对变量之间的p值
for (i in 1:nrow(correlation)) {
  for (j in 1:ncol(correlation)) {
    metab_i <- as.numeric(metab_clr[i, ])
    micro_j <- as.numeric(micro_clr[j, ])
    p_values[i, j] <- cor.test(metab_i, micro_j, method = "spearman")$p.value
  }
}

# 使用p.adjust函数对p值进行多重检验校正
q_values <- p.adjust(p_values, method = "fdr")
significant <- (q_values < 0.25)

# 创建包含所有结果的数据框
correlation_data <- data.frame(Row_Variable = character(),
                               Col_Variable = character(),
                               Correlation = numeric(),
                               P_Value = numeric(),
                               Q_Value = numeric(),
                               Significant = logical(),
                               stringsAsFactors = FALSE)
k <- 1
for (i in 1:nrow(correlation)) {
  for (j in 1:ncol(correlation)) {
    # 获取代谢物和微生物的名称
    metab_name <- rownames(correlation)[i]
    micro_name <- colnames(correlation)[j]
    
    # 将相关性得分、p值、q值和显著性结果添加到数据框中
    correlation_data <- rbind(correlation_data, data.frame(Row_Variable = metab_name,
                                                           Col_Variable = micro_name,
                                                           Correlation = correlation[i, j],
                                                           P_Value = p_values[i, j],
                                                           Q_Value = q_values[k],
                                                           Significant = significant[k]))
    k <- k + 1
  }
}

write.csv(correlation, file = "correlation/IBD2019_ann/spearman_correlation_scores.csv", row.names = FALSE)
# 导出结果为CSV文件
write.csv(correlation_data, file = "correlation/IBD2019_ann/spearman_correlation_scores_with_pq_values.csv", row.names = FALSE)

#----------------------------Sparse CCA-----------------------------------
library("PMA")
set.seed(3189)
x <- t(metab_clr)
z <- t(micro_clr)
out <- CCA(x,z,typex="standard",typez="standard",K=3)
print(out,verbose=TRUE) # To get less output, just print(out)
# Or can use CCA.permute to choose optimal parameter values
perm.out <- CCA.permute(x,z,typex="standard",typez="standard",nperms=10)
print(perm.out)
plot(perm.out)
out1 <- CCA(x,z,typex="standard",typez="standard",K=1,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
           v=perm.out$v.init)
# print(out1)
micro_names <- out1$znames
micro_weights <- out1$v
metab_names <- out1$xnames
metab_weights <- out1$u

# 创建数据框
mic_component_data <- data.frame(Microorganism = micro_names,
                             Microorganism_Weight = micro_weights)

met_component_data <- data.frame(Metabolite = metab_names,
                             Metabolite_Weight = metab_weights)
# 将数据保存为CSV文件
write.csv(mic_component_data, file = "correlation/IBD2019_ann/mic_CCA_weights.csv", row.names = FALSE)
write.csv(met_component_data, file = "correlation/IBD2019_ann/met_CCA_weights.csv", row.names = FALSE)

#----------------------------O2-PLS-----------------------------------
library("OmicsPLS")
library(magrittr) # needs to be run every time you start R and want to use %>%
library(ggplot2)

tax = t(micro_clr)
met = t(metab_clr)

set.seed(123)
crossval_o2m(tax, met, 2:5,1:3,1:3,nr_folds = 10) #10折交叉验证
modelfit<-o2m(tax, met, 5, 2, 2)  #基于交叉验证结果确定成分数目参数
print (modelfit)
#----------------选择前50个微生物-------------------
xj<- loadings(modelfit, "Xjoint", 1:2) %>% abs %>% rowSums
x_df <- data.frame(Feature = colnames(tax), Weight = xj)
write.csv(x_df, file = "correlation/IBD2019_ann/O2pls_mic_loadings.csv", row.names = FALSE)
xj[-(order(xj,decreasing=T)[1:50])] = 0
xj <- sign(xj)
print(xj[xj==1])
selected_x <- xj[xj == 1]
# 创建一个数据框，将选择的元素放入其中
selectedx_df <- data.frame(selected_x)
# 导出数据框到 CSV 文件
write.csv(selectedx_df, file = "correlation/IBD2019_ann/O2pls_selected_x.csv", row.names = TRUE)

#----------------选择前50个代谢物物-------------------
yj<- loadings(modelfit, "Yjoint", 1:2) %>% abs %>% rowSums
y_df <- data.frame(Feature = colnames(met), Weight = yj)
write.csv(y_df, file = "correlation/IBD2019_ann/O2pls_met_loadings.csv", row.names = FALSE)
yj[-(order(yj,decreasing=T)[1:50])] = 0
yj <- sign(yj)
print (yj[yj==1])
selected_y <- yj[yj == 1]
# 创建一个数据框，将选择的元素放入其中
selectedy_df <- data.frame(selected_y)
# 导出数据框到 CSV 文件
write.csv(selectedy_df, file = "correlation/IBD2019_ann/O2pls_selected_y.csv", row.names = TRUE)