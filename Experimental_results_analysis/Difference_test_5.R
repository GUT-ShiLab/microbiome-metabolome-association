#5.处理得到符合检验要求的数据框
#
mic <- read.csv("Difference_test/microbiome.csv",row.names=1)
met <- read.csv("Difference_test/metabolome.csv",row.names=1, check.names = FALSE)
y <- read.csv("Difference_test/label.csv",row.names=1)
met_transposed <- t(met)

# 匹配y的行名和met转置后的行名，并传递标签值给转置后的met的Sub_Class列
for (i in 1:nrow(y)) {
  sample_name <- rownames(y)[i]
  label <- y$Study.Group[i]
  if (sample_name %in% rownames(met_transposed)) {
    met_transposed[sample_name,"Sub_Class"] <- label
  }
}

met_transposed[,"Class"] <- ifelse(met_transposed[,"Sub_Class"] == "UC" | met_transposed[,"Sub_Class"] == "CD", "IBD", "Control")

rownames(met_transposed) <- gsub(".*\\.", "", rownames(met_transposed))
write.csv(t(met_transposed),"Difference_test/met_new.csv",row.names = TRUE)
#==================KW检验以及MW检验代码参见Python文件Difference_test.py=======================


