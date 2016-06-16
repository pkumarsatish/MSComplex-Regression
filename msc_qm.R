# Function file for quadratic model
# Inputs : ms, lidx

msc_qm <- function (ms, lidx)
{
  qm <- vector("list", length(ms$level[[lidx]]$partitionSize))
  ncoeff = 1 + P + (P*(P+1)/2)
  coeff <- matrix(,nrow=length(ms$level[[lidx]]$partitionSize),ncol=ncoeff)

  for (i in 1:length(ms$level[[lidx]]$partitionSize)){
    
    no_points_incell = ms$level[[lidx]]$partitionSize[i]
    y_local <- vector(mode="double", length=no_points_incell)
    x_local <- matrix(, nrow = no_points_incell, ncol = P)
    
    # Getting xlocal and y lcoal
    k = 1;
    for (j in 1:N){
      if (ms$level[[lidx]]$partition[j] == i){
        x_local[k,] = ms$x[j,]
        y_local[k] = ms$y[j]
        k = k+1
      }
    }
    
    myxin <- matrix(,nrow=length(y_local),ncol=ncoeff-1)
    myxin = x_local;
    # Forming Quadratic Structure for model
    for (i1 in 1:P){
      for (j1 in i1:P){
        myxin <- cbind(myxin,x_local[,i1]*x_local[,j1])
      }
    }
    
    #print(cbind(ncoeff,no_points_incell))
    qm[[i]] <- lm(y_local ~ myxin) # Enter linear regre model
    #print(summary(lm[[i]]))
    
    if (sum(is.na(qm[[i]]$coefficients)) != 0){
      qm[[i]]$coefficients[is.na(qm[[i]]$coefficients)] <- 0
      cat("WARNING... Segment #",i,"Least square has multiple solution!!\n")
    }
    
    k = 1
    coeff[i,1]=qm[[i]]$coefficients[k] # constant term
    k = k+1
    
    # Linear Terms
    for (k in 2:(P+1)){
      coeff[i,k]=qm[[i]]$coefficients[k]
    }
    
    # Quadratic terms
    for (k in (P+2):ncoeff){
      #print(cbind(i,k))
      coeff[i,k]=qm[[i]]$coefficients[k]
    }    
  }
  
  coeff
}