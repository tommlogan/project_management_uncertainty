main <- function(){
  # calculate the affect of uncertainty on project completion time estimation and crashing strategies
  # 1. randomly generate projects
  # 2. simulate crashing rules
  # 3. calculate probability of on-time completion
  # 
  # Args:
  #   activity.num: number of activities in the project
  #   link.prob: likelihood of an activity being in a path - this is the density of the upper triangular matrix
  #   
  #   
  # Returns:
  #   
  
  # User defined variables
  activity.num <- 8
  link.prob <- 0.5
  set.seed(15673)
  
  # import libraries
  library(data.table)
  
  # import functions
  source('code/generate_projects_B.R')
  
  
  # Generate random projects  
  projects <- GenerateProjects_B(activity.num, link.prob)
  path.matrix <- projects$path.matrix
  activities <- projects$activities
  
  # Crash strategies 
  #- this code will eventually modify the activity means and variances
  
  
  # Calculate path properties
  project <- CalculatePathProperties(path.matrix, activities)
  
  # Determine probability of on-time completion
  
  
}


CalculatePathProperties <- function(path.matrix, activities){
  # calculate the properties of each path:
  #   mean
  #   variance
  #
  # Return:
  #   data.table of path properties: mean and variance
  
  # number of paths in project
  path.num <- nrow(path.matrix)

  
  # project properties
  path.properties <- data.table(
    mean = lapply(seq(path.num), function(i) sum(activities$means[which(path.matrix[i,]==1)])),
    variance = lapply(seq(path.num), function(i) sum(activities$variances[which(path.matrix[i,]==1)]))
  )
  
  # covariance
  cov.matrix <- diag(path.properties$variance)
  for (i in seq(1,path.num-1)){
    covs <- lapply(seq(i+1,path.num), function(j) {
      activities.shared <- which(colSums(path.matrix[c(i,j),])==2)
      sum(activities$variances[activities.shared])
    })
    cov.matrix[i,(i+1):path.num] <- covs[[1]]
  }
  
  # reflect the covaraince matrix
  cov.matrix[lower.tri(cov.matrix)] <- t(cov.matrix)[lower.tri(cov.matrix)]
  
  # return
  project = list(properties = path.properties, covariance = cov.matrix)
  return(project)
}


CaculateCompletionProbability <- function(path.properties)