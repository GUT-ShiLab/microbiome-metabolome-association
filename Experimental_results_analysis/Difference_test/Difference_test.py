from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.feature_selection import SelectKBest, f_classif
import pandas as pd
from scipy.stats import kruskal, wilcoxon
from scipy.stats import mannwhitneyu
import numpy as np
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import csv

# 读取数据
data = pd.read_csv('met_new.csv').T
data.columns = data.iloc[0]
data = data[1:]

# 提取特征和标签
y = data['Class']
y1= data['Sub_Class']
X = data.drop(columns=['Sub_Class', 'Class'])
# 将数据中的非数值型值转换为数值型值
X = X.apply(pd.to_numeric, errors='coerce')

# 创建一个空列表，用于存储每个特征的 Kruskal-Wallis 检验结果
p_values_kw = []
# 对每个特征执行 Kruskal-Wallis 检验
for column in X.columns:
    # 从数据框架中提取当前特征的值，并按照不同的标签进行分组
    groups = [X[column][y == label] for label in y.unique()]

    # 执行 Kruskal-Wallis 检验
    _, p_value = kruskal(*groups)

    # 将检验结果添加到列表中
    p_values_kw.append(p_value)

# FDR 校正
significant_features_kw = np.array(p_values_kw) < 0.05
significant_metabolites_kw = X.columns[significant_features_kw]

print("Number of significant metabolites after Kruskal-Wallis test:", len(significant_metabolites_kw))

# 创建一个空列表，用于存储每个特征的  Mann-Whitney U 秩和检验结果
p_values_ww = {}

# 初始化显著代谢物列表为 KW 检验的结果
significant_metabolites_ww = list(significant_metabolites_kw)

# 对每个特征执行 Wilcoxon 秩和检验
for column in significant_metabolites_kw:
    # 对每两个不同的子分类标签进行比较
    unique_values = y1.unique()

    # 初始化存储当前特征检验结果的列表
    p_values_ww[column] = []

    for i in range(len(unique_values)):
        for j in range(i + 1, len(unique_values)):
            group_1 = X[column][y1 == unique_values[i]]
            group_2 = X[column][y1 == unique_values[j]]
            # 执行 Wilcoxon 秩和检验
            _, p_value = mannwhitneyu(group_1, group_2)
            # 将检验结果添加到列表中
            p_values_ww[column].append(p_value)

    # 如果所有的 p-value 都不显著，则从显著代谢物列表中移除该特征
    if (p_values_ww[column][1] > 0.01 and p_values_ww[column][2] > 0.01):
        significant_metabolites_ww.remove(column)


print("Number of significant metabolites after  Mann-Whitney U test:", len(significant_metabolites_ww))
significant_metabolites_with_pvalues = []

for metabolite in significant_metabolites_ww:
    # 获取当前代谢物的 KW 检验的 p-value
    kw_pvalue = p_values_kw[X.columns.get_loc(metabolite)]

    # 获取当前代谢物的 WW 检验的 p-values
    ww_pvalues = p_values_ww[metabolite]

    # 添加到列表中
    significant_metabolites_with_pvalues.append([metabolite, kw_pvalue] + ww_pvalues)


# 将显著代谢物及其P-value值写入 CSV 文件
def write_list_to_csv(filename, data):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        # 写入表头
        writer.writerow(["Metabolite", "KW_pvalue", "WW_pvalue_1", "WW_pvalue_2", "WW_pvalue_3"])
        # 写入数据
        for item in data:
            writer.writerow(item)

write_list_to_csv("significant_metabolites_with_pvalues.csv", significant_metabolites_with_pvalues)
