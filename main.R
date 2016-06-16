## Script File for MSR
## Satish Kumar & Vijay Kumar
## MTech SERC, IISc
##

###################################################
### Directory path and Library 
###################################################
setwd('/home/satish/Desktop/myR/camera_estimation') # Working directory
library(msr)
library("R.matlab")
source('~/Desktop/myR/msc_qm.R')
source('~/Desktop/myR/plot.msc.R')
###################################################
### Input parameters 
###################################################
datafile = "camera_estimation.csv" # Data File in CSV Formate
yidx = 1 # collumn index for output

# MS Complex parameter
knn = 50 # knn , generally 5d
#nLevels = 4 # For hirarchical MS Complex
pLevel = 0.1 # Exact percistance value
lidx = 1 # Level choosen for ploting, 1 for single percistance

# Ploting flags
flag1 = 0 # for ploting 2d complex
flag2 = 0 # for ploting persistence plot
flag3 = 1 # For visualization of HD data using regression curve

# Model
model1 = 0 # linear
model2 = 1 # Quadratic

###################################################
### Read the data and Normalize 
###################################################
dataset = read.csv(datafile)
#dataset = fourpeaks(1000)
N = nrow(dataset)
P = ncol(dataset) - 1

X <- matrix(,nrow=N,ncol=P)
Y <- vector("double", N)

if (yidx == 1){
  for (i in 2:(P+1)){
    X[,i-1] <- (dataset[,i] - min(dataset[,i]))/(max(dataset[,i])-min(dataset[,i])) 
  }    
} else if(yidx == P+1){
  for (i in 1:P){
   X[,i] <- (dataset[,i] - min(dataset[,i]))/(max(dataset[,i])-min(dataset[,i]))
  } 
} else{
  print("Error in data!!!")
}

Y <- (dataset[,yidx] - min(dataset[,yidx]))/(max(dataset[,yidx])-min(dataset[,yidx])) 


###################################################
### Create Morse-Smale complex [With SVM]#########
###################################################
ms <- msc.nn.svm(y = Y, x = X, pLevel = pLevel, knn = knn, cost = 1, precompute = F)
#ms <- msc.nn.svm(y = Y, x = X, nLevels = nLevels, knn = knn, cost = 1, precompute = T)


if ((P == 2) && (flag1 == 1)){
##################################################
### plot 2D surface of morse smale complex
###################################################
par(mfcol = c(1, 2))
fancy <- c('#0066CC', '#CCCC00', '#D22905')
colormap = colorRamp(fancy, interpolate = "linear", bias = 0.5 )
colors <- rgb(colormap(Y), maxColorValue = 255)
par(mar = c(8, 5, 5, 4))
plot(X[ , 1], X[ , 2], type = "p", xlab = "", ylab = "", col = colors, pch = 19, cex.axis = 2.5 )
pal <- brewer.pal(9, "Set1")
pal <- colorRampPalette(pal)
pal <- pal(length(ms$level[[lidx]]$mins)) 
colors <- pal[ms$level[[lidx]]$partition]
par(mar=c(8,5,5,4))
plot(X[ , 1], X[ , 2], col = colors, type = "p", xlab = "", ylab = "", pch = 19, cex.axis = 2.5)
}

if ((lidx != 1) && (flag2 == 1)){
###################################################
###  plot persistencies
###################################################
np <- length(ms$persistence)
ms$persistence[np] <- 1
par(mar = c(6, 4, 4, 4) + 3)
plot(ms$persistence, np:1 + 1, xlab = "Persistence Percentage", ylab = "Extrema",
     cex = 2, cex.lab = 2, cex.axis = 2, type = "p", pch = 19)
}

if (flag3 == 1){
###################################################
### Visualize Morse-Smale complex at level lidx
###################################################
ms$predictLevel <- lidx
plot(ms,drawStdDev=T,axesOn=TRUE,span=0.9,plot=TRUE)
}

if (model1 == 1){
  ###################################################
  ### Linear model for regression at level lidx
  ###################################################
  ms$predictLevel <- lidx
  msr.lm <- msc.lm(ms)
  #p1_lm = predict(msr.lm, newdata = dataset)
  #m1_lm = mean((Y-p1_lm)^2)
  #cat("Linear Model:",sqrt(m1_lm),"\n")
}

if (model2 == 1){
  ###################################################
  ### Quadratic model for regression at level lidx
  ###################################################
  msr.qm <- msc_qm(ms,lidx)

###################################################
### predict partition on test data 
###################################################
ms$predictLevel <- lidx
prob <- predict(ms, newdata = X)

###################################################
### Quadratic Model Evaluation 
###################################################
predicted_y <- vector("double", length(N))
predicted_partition <- vector("double", length(N))
error <- vector("double", length(N))

for (i in 1:N){
  partition_test = which( prob == max(prob[i,]), arr.ind = TRUE)
  predicted_partition[i] = partition_test[,2]
  label = partition_test[,2]
  #print(label)
  
  predicted_y[i] = msr.qm[label,1]
  # Linear Terms
  for (k in 1:P){
    predicted_y[i] = predicted_y[i] + (msr.qm[label,k+1]*X[i,k])
  }
  
  k=P+2
  # Quadratic terms
  for (i1 in 1:P){
    for (j1 in i1:P){
      predicted_y[i] = predicted_y[i] + (msr.qm[label,k]*X[i,i1]*X[i,j1])
      k = k+1
    }
  }
  error[i] = predicted_y[i] - Y[i]
}

RMSE = sqrt(sum(error^2)/N)
cat("Quadratic Model:",RMSE)
}
###################################################
### Saving data in csv file for visualization
###################################################
csvdata1 = cbind(ms$y,ms$x[,1],ms$x[,2],ms$level[[lidx]]$partition)
file1 = "train.csv"
writeLines("y,x1,x2,p",file1)
write.table(csvdata1,file1,row.names=F,col.names=F,sep=",",append=T)

if (model2 == 1){
csvdata2 = cbind(predicted_y,X[,1],X[,2],predicted_partition)
file2 = "test.csv"
writeLines("y,x1,x2,p",file2)
write.table(csvdata2,file2,row.names=F,col.names=F,sep=",",append=T)
}