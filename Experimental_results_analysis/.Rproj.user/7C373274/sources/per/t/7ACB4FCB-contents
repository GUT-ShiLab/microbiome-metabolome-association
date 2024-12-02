#8.计算task2： 五折交叉验证实验中各方法综合指标的排名
#
# 轮流读取三个疾病数据集并转换为矩阵
data_IBD <- read.csv("friedman_test/friedman_test_ESRD.csv", header = TRUE, row.names = 1, check.names = FALSE)
data_IBD <- as.matrix(data_IBD)

# 确定指标的方向
metrics_smaller_is_better <- c("RMSE", "NMSE", "MAAPE")
metrics_larger_is_better <- c("TOP50_SCC", "ALL_SCC", "Ratio")

# 创建一个空的矩阵来存储处理后的数据
processed_data <- matrix(nrow = nrow(data_IBD), ncol = ncol(data_IBD))

# 处理数据，根据指标方向取反或直接使用
for (i in 1:nrow(data_IBD)) {
  metric <- rownames(data_IBD)[i]
  if (metric %in% metrics_smaller_is_better) {
    # 越小越好，直接排序
    processed_data[i, ] <- rank(data_IBD[i, ])
  } else if (metric %in% metrics_larger_is_better) {
    # 越大越好，取负值后排序
    processed_data[i, ] <- rank(-data_IBD[i, ])
  }
}

# 执行Friedman检验
result_friedman <- friedman.test(t(processed_data))

# 计算每个方法的平均排名
avg_ranks_IBD <- colMeans(processed_data)

# 打印结果
print("ESRD数据的平均排名：")
print(avg_ranks_IBD)

print("Friedman检验结果：")
print(result_friedman)

#=================计算三个数据集上综合排序后的friedman_test===========

dataset_combine <- read.csv("friedman_test/friedman_test_rank.csv", header = TRUE, row.names = 1, check.names = FALSE)
dataset_combine <- as.matrix(dataset_combine)
# 对前三个指标直接排名，对后三个指标取负值后排名
# 我们使用ifelse函数来根据条件应用不同的排名逻辑
rankings <- t(apply(dataset_combine, 1, function(x) rank(x)))
result_friedman <- friedman.test(rankings)

avg_ranks <- colMeans(rankings)

print("综合排名：")
print(avg_ranks)

print("Friedman检验结果：")
print(result_friedman)
