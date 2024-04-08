library(stringr)
options(max.print=1000000)
load('RES.RData')
baby = read.csv('baby.csv')
baby_lists <- as.character(baby[,1])
col_list <- colnames(read.csv('sample.csv'))[-1]

data_baby <- data_mat[which(1-is.na(match(rownames(data_mat), baby_lists)) == 1), col_list]

result = data.frame()

for (row_length in 1:dim(data_baby)[1]){
  row_matrix <- data.frame(rownames(data_baby)[row_length])
  for (Letter in LETTERS){
    extract_col_list <- col_list[str_detect(col_list, Letter)]
    row_matrix <- cbind(row_matrix, min(data_baby[row_length,extract_col_list], na.rm=TRUE))
  }
  result <- rbind(result, row_matrix)
}
colnames(result) <- c('ID', LETTERS)
write.csv(result, file='result_1.csv', quote=FALSE)
