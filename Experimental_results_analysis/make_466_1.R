#1.预处理论文task-1的数据（微生物数据和代谢物数据）
#
matab_ann <- read.csv("data_process/IBD2019/metabolome_annotation.csv")
genera <- read.table("data_process/IBD2019/genera.tsv",header = TRUE,sep = "\t")
mtb <- read.table("data_process/IBD2019/mtb.tsv",header = TRUE,sep = "\t",check.names = FALSE)
#======================提取代谢物注释，合并相同代谢物数据=======================
columns_to_remove <- which(substr(colnames(mtb), nchar(colnames(mtb)) - 1, nchar(colnames(mtb))) == "NA")
mtb_filtered <- mtb[, -columns_to_remove]
rownames(mtb_filtered) <- mtb_filtered[,1]
mtb_filtered <- mtb_filtered[, -1]

capitalize_first_letter <- function(string) {
  first_letter <- regmatches(string, regexpr("[A-Za-z]", string))
  if (length(first_letter) > 0) {
    string <- gsub(first_letter, toupper(first_letter), string)
  }
  return(string)
}

metabolites_names <- gsub(".*: ", "", colnames(mtb_filtered))
unique_metabolites <- unique(metabolites_names)
combined_data <- data.frame(matrix(ncol = length(unique_metabolites), nrow = nrow(mtb_filtered)))
colnames(combined_data) <- unique_metabolites
for (metabolite in unique_metabolites) {
  matching_columns <- metabolites_names == metabolite

  if (sum(matching_columns) == 1) {
    col_name <- gsub(".*: ", "", metabolite)
    combined_data[, col_name] <- mtb_filtered[, matching_columns]
  } else {

    col_name <- gsub(".*: ", "", metabolite)
    combined_data[, col_name] <- rowSums(mtb_filtered[, matching_columns])
  }
}
rownames(combined_data) <- rownames(mtb_filtered)
column_names <- colnames(combined_data)

column_names_df <- data.frame(Column_Names = column_names)

write.csv(column_names_df, file = "data_process/IBD2019/compound_names.csv", row.names = FALSE)
write.csv(combined_data, file = "data_process/IBD2019/mtb_filtered.csv", row.names = TRUE)
#=============================针对特定代谢物继续进行合并========================
merge_list <- list(
  "lithocholate" = "lithocholic acid",
  "pantothete" = "pantothenate",
  "azelate" = "azelaic acid",
  "N-acetylglutamic acid" = "N-acetylglutamate",
  "7-ketodeoxycholate" = "ketodeoxycholate"
)

for (old_metabolite in names(merge_list)) {
  new_metabolite <- merge_list[[old_metabolite]]
  
  if (old_metabolite %in% colnames(combined_data) && new_metabolite %in% colnames(combined_data)) {
    combined_data[, new_metabolite] <- combined_data[, new_metabolite] + combined_data[, old_metabolite]
    combined_data <- combined_data[, !colnames(combined_data) == old_metabolite, drop = FALSE]
  }else{
   print(paste(old_metabolite,seq="no comluns")) 
  }
}
rownames(combined_data) <- rownames(mtb_filtered)

write.csv(combined_data, file = "data_process/IBD2019/mtb_filtered.csv", row.names = TRUE)
#====================对微生物特征进行处理，合并相同微生物数据===================
rownames(genera) <- genera[,1]
genera <- genera[, -1]
genera_names <- gsub(".*g__([^.]*)\\.?.*", "g__\\1", colnames(genera))
colnames(genera) <- genera_names

extra_letters <- grep("^g__[A-Za-z0-9.]+_[A-Z]$", colnames(genera))

for (feature in colnames(genera)[extra_letters]) {
  target_feature <- gsub("_[A-Z]$", "", feature)
  
  if (target_feature %in% colnames(genera)) {
    genera[, target_feature] <- genera[, target_feature] + genera[, feature]
  } else {
    genera[, target_feature] <- genera[, feature]
  }
}

genera <- genera[, -extra_letters]
row_sums <- rowSums(genera)
print(row_sums)
#==================剔除流行度和平均丰度较低的微生物和代谢物特征=================
average_abundance <- colMeans(genera)
selected_features <- average_abundance > 0.0001
selected_features <- selected_features & colMeans(genera>0) >= 0.1
genera_filtered <- genera[, selected_features]
write.csv(genera_filtered, file = "data_process/IBD2019/genera_filtered.csv", row.names = TRUE)