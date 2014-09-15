CC2cluster <- function(filename){
  dataFrame = read.table(file(filename),header = FALSE,row.names = 1,fill = TRUE)
  
  rownamesCC_raw = c(rownames(dataFrame))
  dm_raw = as.matrix(dataFrame)  
  dm = dm_raw[str_detect(rownamesCC_raw,"_"),]
  rownamesCC = rownamesCC_raw[str_detect(rownamesCC_raw,"_")]
  num_row = length(rownamesCC)
  
#  quantile95 = matrix(data = rep(0,num_row*2),nrow = num_row,ncol = 2)
  quantile95 = rep(0,num_row)
  for(px in 1:num_row){
    quantile95[px] = quantile(dm[px,],0.95,na.rm =TRUE)
 #   quantile95[px,2] = quantile(dm[px,],0.75,na.rm =TRUE)
  }
  # set the inital point 

  
  kmeansR = kmeans(quantile95,centers = c(max(quantile95),min(quantile95)) )
  
  group1 = rownamesCC[kmeansR$cluster == 1 ]
  group2 = rownamesCC[kmeansR$cluster == 2 ]
  
  if(kmeansR$centers[1] > kmeansR$centers[2]){
    cat(c(group1,"\n")) }else{
      cat(c(group2,"\n"))
    } 
}

library("stringr", lib.loc="/usr/lib/R/site-library")
sys_argv = commandArgs()
operation_folder= sys_argv[-c(1:5)]

setwd(operation_folder)
cpd_files = list.files(path = ".",pattern = "*.lst")
numFile = length(cpd_files)

for(filei in 1:numFile){
#cat(c(cpd_files[filei],"\n"))
CC2cluster(cpd_files[filei])
}