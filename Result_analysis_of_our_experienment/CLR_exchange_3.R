#3.预处理task1的数据，执行clr变换
#
library("compositions")
library(ggplot2)
library(reshape2)
metab <- read.csv("correlation/IBD2019_ann/metabolome.csv", row.names = 1, check.names = FALSE)
micro <- read.csv("correlation/IBD2019_ann/microbiome.csv", row.names = 1)
metab[metab == 0] <- 1
# 将微生物数据集中的每个样本的0值替换为当前样本的最小非零值的一半
for (i in 1:ncol(micro)) {
  min_nonzero <- min(micro[micro[,i] != 0, i])
  micro[micro[,i] == 0, i] <- min_nonzero / 2
}
metab_clr <- apply(metab,2,clr)
micro_clr <- apply(micro,2,clr)

write.csv(metab_clr,"correlation/IBD2019_ann/metabolome_clr.csv",row.names = TRUE)
write.csv(micro_clr,"correlation/IBD2019_ann/microbiome_clr.csv",row.names = TRUE)
# write.csv(t(metab_clr),"correlation/IBD2019_ann/microbiome_clr222.csv",row.names = TRUE)

