#2.选择IBD疾病数据集中的PRISM队列作为Task1的研究数据
#
library(readr)
library(dplyr)

genera <- read_csv('data_process/IBD2019/genera_filtered.csv', show_col_types = FALSE) %>% as.data.frame()
mtb <- read_csv('data_process/IBD2019/mtb_filtered.csv', show_col_types = FALSE) %>% as.data.frame()

rownames(mtb) <- mtb[[1]]
mtb <- mtb[,-1]
rownames(genera) <- genera[[1]]
genera <- genera[,-1]

training_genera <- genera[!grepl('Validation', rownames(genera)), ]
testing_genera <- genera[grepl('Validation', rownames(genera)), ]

training_mtb <- mtb[!grepl('Validation', rownames(mtb)), ]
testing_mtb <- mtb[grepl('Validation', rownames(mtb)), ]

training_genera_transposed <- t(training_genera)%>% as.data.frame()
testing_genera_transposed <- t(testing_genera)%>% as.data.frame()

training_mtb_transposed <- t(training_mtb)%>% as.data.frame()
testing_mtb_transposed <- t(testing_mtb)%>% as.data.frame()

write.csv(training_genera_transposed, 'correlation/IBD2019_ann/microbiome.csv',row.names=TRUE)
write.csv(training_mtb_transposed, 'correlation/IBD2019_ann/metabolome.csv', row.names=TRUE)




