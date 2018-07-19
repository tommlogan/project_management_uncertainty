library(data.table)
library(mvtnorm)
library(pbapply)
###
# Crash Strategies
###


Crash.None <- function(activities){
  # Crash strategy:
  # Nothing - just as a point of comparison
  
  return(activities)
}


Crash.CriticalPath.2a <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
  # Crash strategy:
  # Crash the activities on the critical path
  
  path.num <- nrow(path.matrix)
  crashed.activities <- c()
  
  # Calculate path properties
  path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
  
  # identify the critical path
  critical.path.index <- which.max(path.properties$mean)
  
  while (crashing.reserve > 0){
    
    # Calculate path properties
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # identify the critical path
    # critical.path.index <- which.max(path.properties$mean)
    
    # identify activities on the critical path
    critical.path.activities <- which(path.matrix[critical.path.index,]==1)
    
    # remove crashed paths from consideration
    critical.path.activities <- critical.path.activities[! critical.path.activities %in% crashed.activities]
    
    # activity to crash is the one on the critical path with the largest mean
    crash.activity <- critical.path.activities[which.max(activities$means[critical.path.activities])]
    
    # crash the activity (reduce it's mean completion time by 100% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    if (crash.effect.var){
      activities$variances[crash.activity] <- activities$variances[crash.activity] + (activities$cov[crash.activity] * crash.reduction)^2
    }
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(activities)
}


Crash.CriticalPath.2b <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
  # Crash strategy:
  # Crash the activities on the critical path
  
  path.num <- nrow(path.matrix)
  crashed.activities <- c()
  
  while (crashing.reserve > 0){
    
    # Calculate path properties
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # identify the critical path
    critical.path.index <- which.max(path.properties$mean)
    
    # identify activities on the critical path
    critical.path.activities <- which(path.matrix[critical.path.index,]==1)
    
    # remove crashed paths from consideration
    critical.path.activities <- critical.path.activities[! critical.path.activities %in% crashed.activities]
    
    # activity to crash is the one on the critical path with the largest mean
    crash.activity <- critical.path.activities[which.max(activities$means[critical.path.activities])]
    
    # crash the activity (reduce it's mean completion time by 100% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    if (crash.effect.var){
      activities$variances[crash.activity] <- activities$variances[crash.activity] + (activities$cov[crash.activity] * crash.reduction)^2
    }
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(activities)
}


Crash.Optimal <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance){
  # Crash strategy:
  # Crash the activities based on their marginal influence on probability of project being late
  # number of paths
  path.num <- nrow(path.matrix)
  activity.num <- length(activities$means)
  
  crashed.activities <- c()
  probability.contributing.late <- 1
  
  while (crashing.reserve > 0 & sum(probability.contributing.late)>0){
    # Calculate path properties
    project <- CalculatePathProperties(path.matrix, activities)
    
    # calculate the first difference for all paths
    epsi <- crashing.reserve
    path.marginal.completion.prob <- unlist(lapply(seq(1,path.num), function(i) CalculateMarginalEffect(i, epsi, deadline, project, deadline.variance)))
    probability.late.activity <- sapply(seq(1,activity.num), function(i) sum(path.marginal.completion.prob[which(path.matrix[,i]==1)]))
    # probability.late.activity <- sapply(seq(1,activity.num), function(i) WeightedSum(i, path.marginal.completion.prob, path.matrix, project, activities))
    
    
    # which N activity has the highest probability of being late
    probability.late.activity[crashed.activities] <- 0
    probability.contributing.late <- probability.late.activity 
    crash.activity.idx <- which(probability.contributing.late==max(probability.contributing.late)) #which.max(probability.contributing.late)
    # if there are multiple activities with the same max probability of being late, crash the one with the largest mean
    crash.activity <- crash.activity.idx[which.max(activities$means[crash.activity.idx])]
    
    # crash the activity (reduce it's mean completion time by 100% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    if (crash.effect.var){
      activities$variances[crash.activity] <- activities$variances[crash.activity] + (activities$cov[crash.activity] * crash.reduction)^2
    }
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(list(acts = activities, crashed = crashed.activities))
}


CalculateMarginalEffect <- function(path.idx, epsi, deadline, project, deadline.variance){
  
  # difference mean and sd
  path.mean.dif <- project$properties$mean[[path.idx]] - deadline$mean
  if (deadline.variance){
    path.sd <- sqrt(project$properties$variance[[path.idx]] + deadline$variance)
  } else {
    path.sd <- sqrt(project$properties$variance[[path.idx]])
  }
  # Determine probability of on-time completion with added epsi
  prob.ontime.added <- pnorm(q = 0, mean = path.mean.dif, sd = path.sd)
  # prob.ontime.added <- pnorm(q = 0, mean = path.mean.dif + epsi/2, sd = path.sd)
  
  # Determine probability of on-time completion with added epsi
  prob.ontime.less <- pnorm(q = 0, mean = path.mean.dif - epsi, sd = path.sd) #CalculateCompletionProbability(project.less, deadline)
  
  # difference
  first.dif <- (prob.ontime.less - prob.ontime.added)/(2*epsi)
  first.dif[first.dif < 0] <- 0
  
  # first difference divided by standard deviation of path
  marginal.complete.prob <- first.dif#/path.sd
  
  return(marginal.complete.prob)
}


plotCrashEffect <- function(project, deadline){
  # how does reducing the mean change the probability of completion
  project <- CalculatePathProperties(path.matrix, activities)
  path.idx <- which.max(project$properties$mean)
  path.mean <- project$properties$mean[[path.idx]]
  reduce.path.mean <- seq(0, path.mean, 0.01)
  CalcPocCrash <- function(reduce, project, deadline, path.idx){
    # copy project
    project.less <- project
    # reduce mean time
    project.less$properties$mean[[path.idx]] <- project.less$properties$mean[[path.idx]] - reduce
    # calc prob completion
    prob <- CalculateCompletionProbability(project.less, deadline)
    return(prob)
  }
  poc <- lapply(reduce.path.mean, function(reduce) CalcPocCrash(reduce, project, deadline, path.idx))
  plot(reduce.path.mean, poc)
}


Crash.BruteForce <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance){
  # Crash strategy:
  # Crash the activities based on their marginal influence on probability of project being late
  
  # number of paths
  path.num <- nrow(path.matrix)
  activity.num <- length(activities$means)
  
  # activity indices
  activities.index <- seq(1,activity.num)
  crashed.activities <- c()
  
  
  while (crashing.reserve > 0){
    # Calculate path properties
    project <- CalculatePathProperties(path.matrix, activities)
    
    # potential crash activities
    potential.crash.activities <- activities.index[! activities.index %in% crashed.activities]
    
    # calculate the effect on completion probability of crashing each activity
    activity.effect <- unlist(lapply(potential.crash.activities, function(i) ActivityCrashEffect(i, deadline, path.matrix, activities, crash.effect.var, crashing.reserve, deadline.variance)))
    
    # crashing which activity will result in highest on-time completion
    crash.activity <- potential.crash.activities[which.max(activity.effect)]
    
    # crash the activity (reduce it's mean completion time by 50% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    if (crash.effect.var){
      activities$variances[crash.activity] <- activities$variances[crash.activity] + (activities$cov[crash.activity] * crash.reduction)^2
    }
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(list(acts = activities, crashed = crashed.activities))
}


ActivityCrashEffect <- function(i, deadline, path.matrix, activities, crash.effect.var, crashing.reserve, deadline.variance){
  
  # crash the activity
  crash.activity <- i
  crash.reduction <- min(activities$means[i], crashing.reserve)
  # calculate the activities variances resulting from the crash
  if (crash.effect.var){
    activities$variances[crash.activity] <- activities$variances[crash.activity] + (activities$cov[crash.activity] * crash.reduction)^2
  }
  # reduce the mean by the crash
  activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
  
  # recalculate the project information
  project <- CalculatePathProperties(path.matrix, activities)
  
  # Determine probability of on-time completion
  prob <- CalculateCompletionProbability(project, deadline, deadline.variance)
  
  return(prob)
}



Crash.Heuristic <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance){
  # Crash strategy:
  # Crash the activities on the critical path
  
  path.num <- nrow(path.matrix)
  crashed.activities <- c()
  
  while (crashing.reserve > 0){
    
    # Calculate path properties
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # calculate probability exceeding deadline
    path.mean.max <- max(unlist(path.properties$mean))
    difference.mean <- unlist(path.properties$mean) - deadline$mean
    if (deadline.variance){
      difference.sqrt <- sqrt(unlist(path.properties$variance) + deadline$variance)
    } else {
      difference.sqrt <- sqrt(unlist(path.properties$variance))
    }
    
    path.properties$late <- pnorm(q = 0, mean = difference.mean, sd = difference.sqrt, lower.tail = F)
    
    # identify the critical path
    # ignore paths with mean 0, because they cannot be crashed further
    critical.path.index <- which.max(path.properties$late * (path.properties$mean>0))
    
    # identify activities on the critical path
    critical.path.activities <- which(path.matrix[critical.path.index,]==1)
    
    # remove crashed paths from consideration
    critical.path.activities <- critical.path.activities[! critical.path.activities %in% crashed.activities]
    
    # activity to crash is the one on the critical path with the largest mean
    activity.mean <- activities$means[critical.path.activities]
    crash.activity <- critical.path.activities[which.max(activity.mean)]
    
    # crash the activity (reduce it's mean completion time by 100% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    if (crash.effect.var){
      activities$variances[crash.activity] <- activities$variances[crash.activity] + (activities$cov[crash.activity] * crash.reduction)^2
    }
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(list(acts = activities, crashed = crashed.activities))
}

WeightedSum <- function(act.idx, path.marginal.completion.prob, path.matrix, project, activities){
  # path ids
  path.idx <- which(path.matrix[,act.idx]==1)
  # path means
  path.means <- unlist(project$properties$mean[path.idx])
  # fraction activity mean
  act.weight <- activities$means[act.idx]/path.means
  # weighted sum
  influence <- sum(path.marginal.completion.prob[path.idx]*act.weight)
  return(influence)
  
}