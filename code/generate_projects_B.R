GenerateProjects_B <- function(activity.num, path.num, link.prob){
  # generate uTri binary network
  # draw activity mean and standard deviation times
  # 
  # Args:
  #   activity.num: number of activities in the project
  #   path.num: number of paths in the project
  #   link.prob: likelihood of an activity being in a path
  #   
  # Returns:
  #   
  
  # import libraries
  
  
  # generate an incidence matrix. A "1" indicates that an activity is connected
  incidence.matrix <- matrix(rbinom(activity.num * activity.num, 1, link.prob), activity.num, activity.num)
  incidence.matrix[lower.tri(incidence.matrix, T)] <- NA
  
  
  
  # assign activity mean and variance
  activity.means <- runif(path.num, 0, 1)
  activity.variances <- rexp(path.num, 4)
  
  return(list(path.matrix = path.matrix, means = activity.means, var = activity.variances))
  
}
