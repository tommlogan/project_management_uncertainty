 GenerateProjects_B <- function(activity.num, link.prob){
  # generate a project using an activity dependence network, modeled as a random upper triangular binary network
  # draw activity mean and standard deviation times
  # 
  # Args:
  #   activity.num: number of activities in the project
  #   link.prob: likelihood of an activity being in a path - this is the density of the upper triangular matrix
  #   
  # Returns:
  #   path matrix indicating which activities are within a path
  #   list of activity mean and variances
  
  # import libraries
  library(Matrix)
  library(reshape2)
  library(extraDistr)
  # generate an activity network. A one indicates that activity on column is dependent/proceeds an acitivity in row.
  incidence.matrix <- matrix(rbinom(activity.num * activity.num, 1, link.prob), activity.num, activity.num)
  incidence.matrix[lower.tri(incidence.matrix, T)] <- 0
  
  # determine paths from the incidence matrix
  path.matrix <- FindPaths(incidence.matrix)
  
  # assign activity mean and variance
  activity.means <- runif(activity.num, 0, 1)
  activity.coefvariation <- rkumar(activity.num, 5, 1) #runif(activity.num, 0, 1)
  activity.variances <- (activity.coefvariation * activity.means)^2
  activities <- list(means = activity.means, variances = activity.variances, cov = activity.coefvariation)
  
  return(list(path.matrix = path.matrix, activities = activities))
  
}


FindPaths <- function(incidence.matrix){
  # Using a modified depth first search algorithm, returns a path-activity matrix
  # https://stackoverflow.com/questions/20262712/enumerating-all-paths-in-a-directed-acyclic-graph
  # Input: precedence matrix (binary upper triangular representing an activity's dependence on another)
  #       Note: this precedence matrix is a directed acyclic graph
  # Output: matrix with paths on the rows and activities on the columns
  
  activity.num <- nrow(incidence.matrix)
  # initialise the path-acitivity matrix
  path.matrix <- matrix(0, 1, activity.num)
  
  # Get the adajacency list
  adj.list <- melt(incidence.matrix)
  adj.list <- adj.list[adj.list$value > 0, ]
  adj.list <- lapply(seq(activity.num), function(i) adj.list[adj.list$Var1==i,2])
  
  # Identify the starting nodes
  candidate.nodes <- which(colSums(incidence.matrix)==0)
  
  # do the DFS for each of the candidate nodes
  paths <- c()
  for (node in candidate.nodes){
    paths <- c(paths, DepthFirstSearch(adj.list, path = node, paths = list()))
  }
  
  # put into a path-activity matrix
  paths <- stack(setNames(paths, seq_along(paths)))
  path.matrix <- sparseMatrix(as.numeric(paths[,2]), paths[,1], x=1)
  # as.matrix(path.matrix)
  
  
  # remove reduntant/dominated paths, i.e. the ones with are included in other paths
  if (!is.null(dim(path.matrix)[1])){
    path.matrix <- RemoveDominatedPaths(path.matrix)
  }
  
  return(path.matrix)
  
}


DepthFirstSearch <- function(adj.list, path, paths){
  # Depth First search to catch all of the paths through a directed acyclic graph, starting at node
  # Input: adjacency list; a list containing the node to start from; and initialize an empty list for all possible paths
  # Output: nested list of paths
  
  # get the starting node
  init.node <- tail(path, n=1)
  
  # check if this node has any successors
  if (length(adj.list[[init.node]]) > 0 ){
    # if so, repeat the recursion
    for (node in adj.list[[init.node]]){
      path.new <- c(path, node)
      paths <- DepthFirstSearch(adj.list, path.new, paths)
    }
  } else {
    # if not, end the path
    paths <- c(paths, list(path))
  }
  return(paths)
}



RemoveDominatedPaths <- function(path.matrix){
  # Remove dominated or redundant paths i.e. those which are already included in other paths
  #
  
  path.num <- nrow(path.matrix)
  
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
  if (length(dominated) > 0){
    path.matrix <- path.matrix[-dominated, ]
  }
  
  return(path.matrix)
}

