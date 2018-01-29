GenerateProjects_A <- function(activity.num, path.num, link.prob){
# generate matrix of paths and activities
# draw activity mean and standard deviation times
# 
# Args:
#   activity.num: number of activities in the project
#   path.num: number of paths in the project
#   link.prob: likelihood of an activity being in a path
#   
# Returns:
#   path matrix indicating which activities are within a path
#   list of activity mean and variances
 
  # import libraries
  
  
  # generate a path matrix. Each row is a path, each column is an activity. A "1" indicates that that activity is in the path.
  path.matrix <- matrix(rbinom(activity.num * path.num, 1, link.prob), path.num, activity.num)
  
  # remove reduntant/dominated paths, i.e. the ones with are included in other paths
  path.matrix <- removeDominatedPaths(path.matrix, path.num)
  
  # assign activity mean and variance
  activity.means <- runif(path.num, 0, 1)
  activity.variances <- rexp(path.num, 4)
  
  return(list(path.matrix = path.matrix, means = activity.means, var = activity.variances))
      
}


RemoveDominatedPaths <- function(path.matrix, path.num){
  # Remove dominated or redundant paths i.e. those which are already included in other paths
  #
  
  # order rows
  sum.row <- rowSums(path.matrix)
  path.matrix <- path.matrix[sort(sum.row, decreasing = TRUE, index = T)$ix, ]
  # order columns
  sum.cols <- colSums(path.matrix)
  path.matrix <- path.matrix[, sort(sum.cols, decreasing = TRUE, index = T)$ix]
  
  # loop through the rows and check to see if the subsequent rows are dominated
  dominated <- vector()
  # loop the first row
  for (path_i in seq(1, path.num-1)) {
    # loop subsequent rows
    for (path_j in seq(path_i + 1, path.num)) {
      
      # take the difference between path activities
      difference <- path.matrix[path_j, ] - path.matrix[path_i, ]
      
      # a path is not dominated if any is true
      is.dominated <- all(difference <= 0)
      
      if (is.dominated) {
        dominated <- c(dominated, path_j)
      }
    }
  }
  
  # remove dominated paths
  path.matrix <- path.matrix[-dominated, ]
  
  return(path.matrix)
}


# debugging
activity.num <- 5
path.num <- 5
link.prob <- 0.5


