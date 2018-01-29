#' ---
#' title: ""
#' author: ""
#' date: ""
#' ---
# ACPA method
# An approximation for project completion time
# Available at https://github.com/tommlogan/acpa

# Import Libraries ---------------------------
library(data.table)
source('calc_clark.R')

# Define main function ---------------------------
AcpaApproximation <- function(x, path.properties, path.correlations){
  # INPUTS:
    # x is the quantile
    # path.properties is a data.table of mean and standard deviations
    # path.correlations is a data.table of correlations between the paths
  # assign k
  k <- dim(path.properties)[1]
  # calculate the g value
  path.properties[, 'g'] <- path.properties[, CalcLittleg(x,mean,stdev)]
  
  # sort the paths by their g value
  path.properties <- path.properties[order(-g)]
  # sort the correlation matrix
  path.correlations.sorted <- copy(path.correlations)
  path.order <- path.properties$path.number
  j <- 0
  for (i in path.order){
    j <- j + 1
    path.correlations.sorted[j,] <- path.correlations[i,][path.order]
  }
  
  # calculate the exponent
  alpha.values <- CalculateCorrelations(k, path.properties,path.correlations.sorted)
  
  # calculate the dependency with the gombel copula
  zk <- NestedCopula(k,path.properties,alpha.values)
  
  return(exp(-zk))
}


# Define sub functions ---------------------------
CalcLittleg <- function(t, path.mean, path.stdev){
  
  # calculate the lower probability of the normcdf
  prob.less <- pnorm(t,path.mean,path.stdev)
  
  # take the negative log
  g <- -log(prob.less)
  
  return(g)
}


SelectExponents <- function(path.properties,path.correlations){
  # Calculate the exponent which the dependencies are raised to
  
  # determine the path order
  path.order <- path.properties[,path.number]
  
  # calculate the exponent values from the correlations
  
  exponents <- CalculateExponents(path.correlations)
  
  # select the values
  p.values <- c(exponents[path.order[1],path.order[2]],
                exponents[path.order[1],path.order[3]])
  
  return (p.values)
}


CalculateExponents <- function(path.correlations){
  # calculate p matrix

  # calculate p
  p <- 1./(1-asin(path.correlations)*(2/pi));
  
  # set diag to zero
  diag(p) <- 0;

  return(p)
}



NestedCopula <- function(k,path.properties,alpha.values){
  # this function calculates the nested gombel copula
  #
  # init z list
  z <- vector("numeric", k) 
  # init the first value of the list
  z[1] <- path.properties[1,g]
  # loop through this and calculate each value for z
  for (i in 2:k){
    z[i] <- (z[i-1]^alpha.values[i-1] + path.properties[i,g]^alpha.values[i-1]) ^ (1/alpha.values[i-1])
  }
  
  return(z[k])
  
}