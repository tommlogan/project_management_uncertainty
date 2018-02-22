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
  set.seed(15673)
  
  # User defined variables
  # activity.num <- 8
  # link.prob <- 1/3
  deadline <- list(critical.path.percentile = 0.9, coeff.variation = seq(0,1, 0.1))
  crash.number <- 4
  number.cores <- 11
  
  # import libraries
  library(data.table)
  library(mvtnorm)
  library(pbapply)
  library(ggplot2)
  library(ggthemes)
  
  # import functions
  source('code/generate_projects_B.R')
  
  # init results dataframe
  df <- NULL
  # loop activity and precedence density
  for (i in seq(1,3)){
    activity.num <- sample(5:15,1)
    link.prob <- runif(1, 0.2, 0.4)
      
    
    # Generate random project 
    projects <- GenerateProjects_B(activity.num, link.prob)
    path.matrix <- projects$path.matrix
    activities <- projects$activities
    
    # Calculate path properties
    path.num <- nrow(path.matrix)
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # Calculate deadline mean
    critical.path.index <- which.max(path.properties$mean)
    deadline$mean <- qnorm(deadline$critical.path.percentile, mean = path.properties$mean[[critical.path.index]], sd = sqrt(path.properties$variance[[critical.path.index]]))
    
    # Crash strategies 
    strategy.names <- c('critical.path', 'probability.lateness')
    
    # loop through the deadline's coefficient of variation, calculate probability of on-time completion
    results <- pblapply(deadline$coeff.variation, function(CoV) VaryCoefVariation(CoV, deadline, strategy.names, path.matrix, activities, link.prob), cl = number.cores)
    df <- rbind(df, rbindlist(results, use.names=TRUE, fill=FALSE))
  }

  # save results
  variate.str <- 'CoV'
  fwrite(df, file = paste0('data/results/vary_',variate.str,'.csv'))
  
  # plot the results
  PlotRibbon(df, variate.str, strategy.names, log_axis = F)
  
}


VaryCoefVariation <- function(CoV, deadline, strategy.names, path.matrix, activities, link.prob){
  # calculate the probability of on-time completion with different crashing strategies for a given coefficient of variation in the deadline
  # 
  # Return:
  #   Coefficient of Variation for the deadline
  #   Probability of on-time completion for each of the crashing strategies
  
  # calculate the deadline variance
  deadline$variance <- (deadline$mean * CoV)**2
  
  # loop through the crash strategies and return probability of on-time completion
  prob.ontime <- lapply(strategy.names, function(crash.strategy) CalcCrashCompletionProb(crash.strategy, activities, crash.number, path.matrix, deadline))
  
  # named results
  prob.ontime <- unlist(prob.ontime)
  prob.ontime <- setNames(as.list(prob.ontime), strategy.names)
  
  # return the following:
  # coefficient of variation | number of activities | number of paths | link density | probability of completion time for each crash strategy
  results <- list(CoV = CoV, num_acts = ncol(path.matrix), num_paths = nrow(path.matrix), precedence_density = link.prob)
  results <- c(results, prob.ontime)
  
  return(results)
}


CalcCrashCompletionProb <- function(crash.strategy, activities, crash.number, path.matrix, deadline){
  # calculate the probability of on-time completion using a crashing strategy
  
  # crash the activity and update the activities list
  if (crash.strategy == 'critical.path'){
    activities <- Crash.CriticalPath(activities, crash.number, path.matrix, deadline)
  } else {
    activities <- Crash.LatenessProbability(activities, crash.number, path.matrix, deadline)
  }
  
  # update the path properties
  project <- CalculatePathProperties(path.matrix, activities)
  
  # Determine probability of on-time completion
  prob.ontime <- CalculateCompletionProbability(project, deadline) 
}


CalculatePathProperties <- function(path.matrix, activities){
  # calculate the properties of each path:
  #   mean
  #   variance
  #   covariance between paths
  #
  # Return:
  #   data.table of path properties: mean and variance
  #   cov.matrix is the covariance matrix between paths
  
  # number of paths in project
  path.num <- nrow(path.matrix)
  
  # project properties
  path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
  
  # covariance
  cov.matrix <- diag(path.properties$variance)
  for (i in seq(1,path.num-1)){
    covs <- lapply(seq(i+1,path.num), function(j) {
      activities.shared <- which(colSums(path.matrix[c(i,j),])==2)
      sum(activities$variances[activities.shared])
    })
    cov.matrix[i,(i+1):path.num] <- unlist(covs)
  }
  
  # reflect the covariance matrix
  cov.matrix[lower.tri(cov.matrix)] <- t(cov.matrix)[lower.tri(cov.matrix)]
  
  # return
  project = list(properties = path.properties, covariance = cov.matrix)
  return(project)
}


CalculatePathMeanVariance <- function(path.matrix, activities, path.num){
  # calculate the properties of each path:
  #   mean
  #   variance
  #
  # Return:
  #   data.table of path properties: mean and variance

  # project properties
  path.properties <- data.table(
    mean = lapply(seq(path.num), function(i) sum(activities$means[which(path.matrix[i,]==1)])),
    variance = lapply(seq(path.num), function(i) sum(activities$variances[which(path.matrix[i,]==1)]))
  )
  return(path.properties)
}


CalculateCompletionProbability <- function(project, deadline){
  # Calculate the joint multivariate normal probability of the project finishing on-time
  #
  # Return:
  #   probability of on-time completion
  
  # calculate the difference between the paths and the deadlien
  
  # covariance
  difference.covariance <- project$covariance + deadline$variance
  
  # path means
  difference.mean <- unlist(project$properties[,mean]) - deadline$mean
  
  # multivariate normal probability
  prob.ontime <- pmvnorm(upper = 0, mean = difference.mean,   sigma = difference.covariance)
  
  return(prob.ontime)
}


PlotRibbon <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  
  # colors
  cols <- c("#f45b5b", "#7cb5ec")
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- c(1, 2)
  names(ltype) <- strategy.names

  # plot  
  plt <- ggplot() + theme_hc()
  
  i <- 1
  for (strategy in strategy.names){
    # ribbon
    plt <- plt + stat_summary(data = df, aes_string(x = variate.str, y = strategy), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.25), fun.ymax = function(x) quantile(x, 0.75)) + 
    # mean
      stat_summary(data = df, aes_string(x = variate.str, y = strategy, color = shQuote(strategy.names[[i]]), linetype = shQuote(strategy.names[[i]])), geom="line", fun.y=median, size = 1)
    # iterate
    i <- i + 1
  }

      # plotting
  plt <- plt +
    xlab(variate.str) +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Probability of Completion Before the Deadline') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) +
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )
  if (log_axis){
    plt <- plt + scale_x_continuous(trans='log10', breaks = unique(df$x.list))
  }
  
  # guides(colour = guide_legend(override.aes = list(size=1)))
  # plot
  print(plt)
  
  # save figure
  fig.width <- 3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', variate.str, '_vary.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}

###
# Crash Strategies
###

Crash.CriticalPath <- function(activities, crash.number, path.matrix, deadline){
  # Crash strategy:
  # Crash the activities on the critical path
  
  for (i in seq(1, crash.number)){
    # Calculate path properties
    path.num <- nrow(path.matrix)
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # identify the critical path
    critical.path.index <- which.max(path.properties$mean)
  
    # identify activities on the critical path
    critical.path.activities <- which(path.matrix[critical.path.index,]==1)
    
    # activity to crash is the one on the critical path with the largest mean
    crash.activity <- critical.path.activities[which.max(activities$means[critical.path.activities])]
    
    # crash the activity (reduce it's mean completion time by 50%)
    activities$means[crash.activity] <- activities$means[crash.activity]/2
  }
  
  return(activities)
}

Crash.LatenessProbability <- function(activities, crash.number, path.matrix, deadline){
  # Crash strategy:
  # Crash the activities based on their probability of contributing to late project 
  
  for (i in seq(1, crash.number)){
    # Calculate path properties
    path.num <- nrow(path.matrix)
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # calculate the probability that each path exceeds the deadline
    difference.mean = unlist(path.properties[,mean]) - deadline$mean
    difference.sd = sqrt(unlist(path.properties[,variance]) + deadline$variance)
    probability.late.path <- sapply(seq(1,length(difference.sd)), function(i) pnorm(0, mean = difference.mean[i], sd = difference.sd[i], lower.tail = FALSE))
    
    # for each activity, sum the probability of lateness for each path it is on
    probability.late.activity <- sapply(seq(1,length(activities$means)), function(i) sum(probability.late.path[which(path.matrix[,i]==1)]))
    
    # which activity has the highest probability of being late
    crash.activity <- which.max(probability.late.activity)
    # if multiple activities have the highest lateness probability, select the one with the largest mean
    if (length(crash.activity) > 1){
      crash.activity <- crash.activity[which.max(activities$means[crash.activity])]
    }
    
    # crash the activity (reduce it's mean completion time by 50%)
    activities$means[crash.activity] <- activities$means[crash.activity]/2
  }
  
  return(activities)
}
