library(Matrix)
library(readr)
library(dplyr)

# read prefiltered Matrix
cl_indata <- Matrix::readMM("/mnt/workspace_stud/stud10/output/esophagus_mucosa/cl_indata_esophagus_mucosa.mtx")

# read peaks.tsv file 
peaks <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/peaks.tsv', 
                  col.names = 'peaks', sep = ',', row.names = NULL, header = FALSE)
head(peaks)

cluster <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/cluster.tsv', 
                    col.names = 'cluster', sep = ',', row.names = NULL, header = FALSE)
head(cluster)

# set col and rown names of cl_indata
row.names(cl_indata) <- peaks$peaks
colnames(cl_indata) <- cluster$cluster
head(cl_indata)

# create 2 variables to check highest value of cluster as well as amount of clusters.
max <- max(unique(cluster$cluster))
len <- length(unique(cluster$cluster))

# create variables with the name v<number> and assign them with 0. 
vals <- vector("list", length = len)
assign("vals", rep(0, len))
n    <- length(vals)
lhs  <- paste("v",    0:max,     sep="")
rhs  <- paste("vals[",0:max,"]", sep="")
eq   <- paste(paste(lhs, rhs, sep="<-"), collapse=";")
eval(parse(text=eq))

a <- matrix(nrow=1 ,ncol= 19, byrow = TRUE)
for(r in 1:nrow(cl_indata)) {
  for(c in 1:ncol(cl_indata)) {
    if (cl_indata[r, c] == 1){
      a[r, 1] = rownames(cl_indata)[r]
      if(colnames(cl_indata)[c] == "0"){
        v0 <- v0 + 1
      } else if (colnames(cl_indata)[c] == '1'){
        v1 <- v1 + 1
      } else if (colnames(cl_indata)[c] == '2'){
        v2 <- v2 + 1
      } else if (colnames(cl_indata)[c] == '3'){
        v3 <- v3 + 1
      } else if (colnames(cl_indata)[c] == '4'){
        v4 <- v4 + 1
      } else if (colnames(cl_indata)[c] == '5'){
        v5 <- v5 + 1
      } else if (colnames(cl_indata)[c] == '6'){
        v6 <- v6 + 1
      } else if (colnames(cl_indata)[c] == '7'){
        v7 <- v7 + 1
      } else if (colnames(cl_indata)[c] == '8'){
        v8 <- v8 + 1
      } else if (colnames(cl_indata)[c] == '9'){
        v9 <- v9 + 1
      } else if (colnames(cl_indata)[c] == '10'){
        v10 <- v10 + 1
      } else if (colnames(cl_indata)[c] == '11'){
        v11 <- v11 + 1
      } else if (colnames(cl_indata)[c] == '12'){
        v12 <- v12 + 1
      } else if (colnames(cl_indata)[c] == '13'){
        v13 <- v13 + 1
      } else if (colnames(cl_indata)[c] == '14'){
        v14 <- v14 + 1
      } else if (colnames(cl_indata)[c] == '15'){
        v15 <- v15 + 1
      } else if (colnames(cl_indata)[c] == '16'){
        v16 <- v16 + 1
      } else if (colnames(cl_indata)[c] == '17'){
        v17 <- v17 + 1
      }
      
    }  
  }
  a <- suppressWarnings(rbind(a, c(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17))) 
  # set all variables to 0 again
  eval(parse(text=eq))
}


#create variable to find number of rows and substract 1 
x <- nrow(a) - 1
# Move Column one row down
a[, 1]  <- matrix(c(NA, a[1:x, 1]))
# Change colnames and rownames 
row.names(a) <- a[, 1]
colnames(a) <- c("Peak", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                 "13", "14", "15", "16", "17")

# delete row with NA and column Peak
a <- a[complete.cases(a), ]
a <- a[,-1]

# convert matrix to data frame. 
a <- as.data.frame(a)

# store data frame in dir /stud10. 
write.csv(a, '/mnt/workspace_stud/stud10/cluster_peaks.csv')