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
  crashing.reserve <- 1
  crash.effect.var <- -1 #if negative, the variance increases due to crashing
  number.cores <- 11
  
  # import libraries
  library(data.table)
  library(mvtnorm)
  library(pbapply)
  library(ggplot2)
  library(ggthemes)
  library(RColorBrewer)
  library(directlabels)
  # library(reshape)
  
  
  # import functions
  source('code/generate_projects_B.R')
  
  # init results dataframe
  df <- NULL
  # loop activity and precedence density
  sim_num = 1000
  for (i in seq(1,sim_num)){
    activity.num <- sample(5:10,1)
    link.prob <- runif(1, 0.05, 0.4)
    # Generate random project 
    path_dim <- NULL
    while (is.null(path_dim)){
      projects <- GenerateProjects_B(activity.num, link.prob)
      path.matrix <- projects$path.matrix
      activities <- projects$activities
      path_dim <- dim(path.matrix)[1]
    }
    
    # Calculate path properties
    path.num <- nrow(path.matrix)
    path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num)
    
    # Calculate deadline mean
    critical.path.index <- which.max(path.properties$mean)
    deadline$mean <- qnorm(deadline$critical.path.percentile, mean = path.properties$mean[[critical.path.index]], sd = sqrt(path.properties$variance[[critical.path.index]]))
    
    # Crash strategies 
    strategy.names <- c('brute', 'critical.path', 'heuristic', 'limiting', 'none')
    
    # loop through the deadline's coefficient of variation, calculate probability of on-time completion
    results <- pblapply(deadline$coeff.variation, function(CoV) VaryCoefVariation(CoV, deadline, strategy.names, path.matrix, activities, link.prob, crashing.reserve, crash.effect.var), cl = number.cores)
    df <- rbind(df, rbindlist(results, use.names=TRUE, fill=FALSE))
  }
  
  # save results
  variate.str <- 'CoV'
  fwrite(df, file = paste0('data/results/vary_',variate.str,'.csv'))
  # df <- fread(paste0('data/results/vary_',variate.str,'.csv'))
  # plot the results
  PlotMean(df, variate.str, strategy.names, log_axis = F)
  PlotDifference(df, variate.str, strategy.names, log_axis = F)
  
  # plot the forecast bias caused by not considering deadline uncertainty
  PlotForecastBias(df, sim_num, strategy.names, variate.str)
}


VaryCoefVariation <- function(CoV, deadline, strategy.names, path.matrix, activities, link.prob, crashing.reserve, crash.effect.var){
  # calculate the probability of on-time completion with different crashing strategies for a given coefficient of variation in the deadline
  # 
  # Return:
  #   Coefficient of Variation for the deadline
  #   Probability of on-time completion for each of the crashing strategies
  
  # calculate the deadline variance
  deadline$variance <- (deadline$mean * CoV)**2
  
  # loop through the crash strategies and return probability of on-time completion
  prob.ontime <- lapply(strategy.names, function(crash.strategy) CalcCrashCompletionProb(crash.strategy, activities, crashing.reserve, path.matrix, deadline, crash.effect.var))
  
  # named results
  prob.ontime <- unlist(prob.ontime)
  prob.ontime <- setNames(as.list(prob.ontime), strategy.names)
  
  # return the following:
  # coefficient of variation | number of activities | number of paths | link density | probability of completion time for each crash strategy
  results <- list(CoV = CoV, num_acts = ncol(path.matrix), num_paths = nrow(path.matrix), precedence_density = link.prob)
  results <- c(results, prob.ontime)
  
  return(results)
}


CalcCrashCompletionProb <- function(crash.strategy, activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
  # calculate the probability of on-time completion using a crashing strategy
  
  # crash the activity and update the activities list
  if (crash.strategy == 'critical.path'){
    activities <- Crash.CriticalPath(activities, crashing.reserve, path.matrix, deadline, crash.effect.var)
  } else if (crash.strategy == 'heuristic'){
    activities <- Crash.Optimal(activities, crashing.reserve, path.matrix, deadline, crash.effect.var)
  } else if (crash.strategy == 'none'){
    activities <- Crash.None(activities)
  } else if (crash.strategy == 'limiting'){
    activities$means <- activities$means * 0
  } else if (crash.strategy == 'brute'){
    activities <- Crash.BruteForce(activities, crashing.reserve, path.matrix, deadline, crash.effect.var)
  }
  
  # update the path properties
  project <- CalculatePathProperties(path.matrix, activities)
  
  # Determine probability of on-time completion
  prob.ontime <- CalculateCompletionProbability(project, deadline) 
  
  return(prob.ontime)
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
  #   - covariance is sum of variance of shared actvities. the deadline is now considered a "shared activity"
  #     so can add it to the covariance
  difference.covariance <- project$covariance + deadline$variance
  
  # path means
  difference.mean <- unlist(project$properties[,mean]) - deadline$mean
  
  # remove paths that have no time
  difference.mean <- difference.mean[apply(difference.covariance,1,function(difference.covariance) !all(difference.covariance==0))] 
  difference.covariance <- difference.covariance[apply(difference.covariance,1,function(difference.covariance) !all(difference.covariance==0)),apply(difference.covariance,2,function(difference.covariance) !all(difference.covariance==0))] 
  
  # multivariate normal probability
  prob.ontime <- pmvnorm(upper = 0, mean = difference.mean,   sigma = difference.covariance, algorithm = GenzBretz(abseps = 1e6))
  
  return(prob.ontime)
}


PlotRibbon <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  # plot  
  plt <- ggplot() + theme_hc()
  
  i <- 1
  for (strategy in strategy.names){
    plt <- plt + 
      # ribbon
      stat_summary(data = df, aes_string(x = variate.str, y = strategy), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.25), fun.ymax = function(x) quantile(x, 0.75)) + 
      # mean
      stat_summary(data = df, aes_string(x = variate.str, y = strategy, color = shQuote(strategy.names[[i]]), linetype = shQuote(strategy.names[[i]])), geom="line", fun.y=median, size = 1)
    # geom_dl(aes_string(x = variate.str, y = strategy),  method = list(dl.combine("first.points", "last.points"), cex = 0.8))
    # iterate
    i <- i + 1
  }
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Probability of Completion \n Before the Deadline') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) +
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )   
  # scale_colour_discrete(guide = 'none') +
  # scale_x_discrete(expand=c(0, 1)) +
  # geom_dl(aes(label=strategy), method = list(dl.combine("first.points", "last.points"), cex = 0.8))
  
  if (log_axis){
    plt <- plt + scale_x_continuous(trans='log10', breaks = unique(df$x.list))
  }
  
  # guides(colour = guide_legend(override.aes = list(size=1)))
  # plot
  print(plt)
  
  # save figure
  fig.width <- 7  #3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_', variate.str, '-vary.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}


PlotMean <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  # library(reshape)
  # reshape the df 
  df2 <- melt(df, id.vars=c("CoV", "num_acts",'num_paths','precedence_density'), variable.name = "strategy", value.name = "prob")
  # aggregate(df2$prob, by = list(df2[,c('CoV','strategy')]), mean)
  # mydt.mean <- df2[,lapply(.SD,mean,na.rm=TRUE),by=c(Cov, strategy),.SDcols=colstoavg]
  df3 <- df2[, lapply(.SD, mean), by = .(CoV, strategy)]
  
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  # plot  
  plt <- ggplot(df3, aes(x = CoV, y = prob, group = strategy, colour = strategy)) + theme_hc()+
    # stat_summary(aes(group=strategy), fun.y=mean, colour="red", geom="line",group=1)
    geom_line() + 
    # stat_summary(data = df, geom="line", fun.y=median, size = 1) + 
    # scale_colour_discrete(guide = 'none') +
    # scale_x_discrete(expand=c(0, 1)) +
    geom_dl(aes(label = strategy), method = list(dl.combine( "last.points"), cex = 0.8), size = 1) 
  
  
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Probability of Completion \n Before the Deadline') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) + #, expand=c(0, 0.2)) +
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    # expand_limits(0,1.2) + 
    theme(
      legend.position = "none",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )   
  
  # plot
  print(plt)
  
  # save figure
  fig.width <- 7  #3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_', variate.str, '-vary_mean.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}


PlotDifference <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  
  # reshape the df 
  df$heur_crit <- df$heuristic - df$critical.path
  df$heur_brute <- df$brute - df$heuristic
  df$id <- rep(seq(1,dim(df)[1]/length(unique(df$CoV))),each = length(unique(df$CoV)))
  
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  i <- 1
  strategy <- 'heur_crit'
  # plot  
  plt <- ggplot() + theme_hc()
  plt <- plt + 
    # ribbon
    stat_summary(data = df, aes_string(x = variate.str, y = strategy), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95)) +
    # lines
    # geom_line(data = df, aes_string(x = variate.str, y = strategy, group = 'id', color = 'num_paths'), size = 1, alpha = 0.5) +
    # median
    stat_summary(data = df, aes_string(x = variate.str, y = strategy), geom="line", fun.y=median, size = 1)
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Forecast bias when crashing: \n Heuristic vs Critical Path') +
    coord_cartesian(ylim=c(-0.1,0.1)) +
    # scale_y_continuous(breaks = seq(-10,10, by = 0.05)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) +
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    theme(
      legend.position = "none",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )
  # 
  
  # plot
  print(plt)
  
  # save figure
  fig.width <- 7/2  #3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_critical.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
  
  ###
  # Brute
  ###
  # plot 
  plt <- ggplot() + theme_hc()
  plt <- plt + 
    # ribbon
    stat_summary(data = df, aes_string(x = variate.str, y = 'heur_brute'), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.25), fun.ymax = function(x) quantile(x, 0.75)) + 
    # mean
    stat_summary(data = df, aes_string(x = variate.str, y = 'heur_brute', color = shQuote(strategy.names[[i]]), linetype = shQuote(strategy.names[[i]])), geom="line", fun.y=median, size = 1)
  
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Forecast bias when crashing: \n Brute Force vs Heuristic') +
    # ylim(-0.1,0.1) + 
    coord_cartesian(ylim=c(-0.1,0.1)) +
    # scale_y_continuous(breaks = seq(-0.1,0.1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) +
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    theme(
      legend.position = "none",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )
  
  
  # plot
  print(plt)
  
  # save figure
  fig.width <- 7/2  #3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_brute.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
  
}


PlotExplore <- function(){
  
  df$heur_crit <- df$heuristic - df$critical.path
  df$id <- rep(seq(1,dim(df)[1]/length(unique(df$CoV))),each = length(unique(df$CoV)))
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  ggplot(df, aes(CoV, num_paths, z = heur_crit, color = heur_crit)) + stat_density2d()
  
  strategy = 'heur_crit'
  variate.str = 'CoV'
  i = 1
  ggplot(df[df$num_paths == 4]) + 
    # ribbon
    stat_summary(aes_string(x = variate.str, y = strategy), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95)) +
    # lines
    geom_line(aes_string(x = variate.str, y = strategy, group = 'id', color = 'num_paths'), size = 1, alpha = 0.5) +
    # median
    stat_summary(aes_string(x = variate.str, y = strategy), geom="line", fun.y=median, size = 1)
  
}


PlotForecastBias <- function(df, sim_num, strategy.names, variate.str){
  # calculate forecast bias
  # for 
  # get the first prob for each of the sims
  entry <- dim(df)[1]/sim_num
  zero_vals <- df[seq(1,dim(df)[1],entry),'none'][[1]]
  zero_list <- rep(zero_vals, each=entry)
  
  df$bias = zero_list - df$none
  # plot the ribbon results
  
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  # plot  
  plt <- ggplot() + theme_hc()
  
  i <- 1
  # for (strategy in strategy.names){
  plt <- plt + 
    # ribbon
    stat_summary(data = df, aes_string(x = variate.str, y = 'bias'), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95)) + 
    # mean
    stat_summary(data = df, aes_string(x = variate.str, y = 'bias', color = shQuote(strategy.names[[i]]), linetype = shQuote(strategy.names[[i]])), geom="line", fun.y=median, size = 1)
  # iterate
  # i <- i + 1
  # }
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Forecast Bias When \n Ignoring Deadline Uncertainty') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) +
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    theme(
      legend.position = "none",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )
  
  
  # guides(colour = guide_legend(override.aes = list(size=1)))
  # plot
  print(plt)
  
  # save figure
  fig.width <- 7  #3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_forecast-bias.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}


###
# Crash Strategies
###


Crash.None <- function(activities){
  # Crash strategy:
  # Nothing - just as a point of comparison
  
  return(activities)
}

Crash.CriticalPath <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
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
    activities$variances[crash.activity] <- (activities$cov[crash.activity] * (activities$means[crash.activity] - crash.reduction*crash.effect.var))^2
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(activities)
}


Crash.Optimal <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
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
    epsi <- 1
    path.marginal.completion.prob <- unlist(lapply(seq(1,path.num), function(i) CalculateMarginalEffect(i, epsi, deadline, project)))
    # act.marginal.pooc <- unlist(lapply(seq(1,activity.num), function(i) CalculateMarginalEffectActivity(i, epsi, deadline, project, path.matrix, activities)))
    # for each activity, 1 - product of the probability of all paths being early = prob of at least one path being late
    # probability.late.activity <- act.marginal.pooc
    probability.late.activity <- sapply(seq(1,activity.num), function(i) sum(path.marginal.completion.prob[which(path.matrix[,i]==1)]))
    # probability.late.activity <- sapply(seq(1,activity.num), function(i) WeightedSum(i, path.marginal.completion.prob, path.matrix, project, activities))
    
    # weight this by the mean of each activity
    # prob.late.means <- probability.late.activity * activities$means
    
    # which N activity has the highest probability of being late
    probability.late.activity[crashed.activities] <- 0
    probability.contributing.late <- probability.late.activity #* activities$means
    crash.activity.idx <- which(probability.contributing.late==max(probability.contributing.late)) #which.max(probability.contributing.late)
    # if there are multiple activities with the same max probability of being late, crash the one with the largest mean
    crash.activity <- crash.activity.idx[which.max(activities$means[crash.activity.idx])]
    
    # crash the activity (reduce it's mean completion time by 100% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    activities$variances[crash.activity] <- (activities$cov[crash.activity] * (activities$means[crash.activity] - crash.reduction*crash.effect.var))^2
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(activities)
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


CalculateMarginalEffect <- function(path.idx, epsi, deadline, project){
  
  # subtract epsi
  project.less <- project
  project.less$properties$mean[[path.idx]] <- project.less$properties$mean[[path.idx]] - epsi
  
  # add epsi
  project.add <- project
  project.add$properties$mean[[path.idx]] <- project.add$properties$mean[[path.idx]] + epsi
  
  # Determine probability of on-time completion with added epsi
  prob.ontime.added <- CalculateCompletionProbability(project.add, deadline)
  
  # Determine probability of on-time completion with added epsi
  prob.ontime.less <- CalculateCompletionProbability(project.less, deadline)
  
  # difference
  first.dif <- (prob.ontime.less - prob.ontime.added)/(2*epsi)
  first.dif[first.dif < 0] <- 0
  
  # first difference divided by standard deviation of path
  marginal.complete.prob <- first.dif #/sqrt(project$properties$variance[[path.idx]])
  
  return(marginal.complete.prob)
}



CalculateMarginalEffectActivity <- function(activity.idx, epsi, deadline, project, path.matrix, activities){
  
  # subtract epsi
  activities.less <- activities
  activities.less$means[[activity.idx]] <- activities.less$means[[activity.idx]] - epsi
  project.less <- CalculatePathProperties(path.matrix, activities.less)
  # project.less$properties$mean[[path.idx]] <- project.less$properties$mean[[path.idx]] - epsi
  
  # add epsi
  activities.add <- activities
  activities.add$means[[activity.idx]] <- activities.add$means[[activity.idx]] - epsi
  project.add <- CalculatePathProperties(path.matrix, activities.add)
  # project.add$properties$mean[[path.idx]] <- project.add$properties$mean[[path.idx]] + epsi
  
  # Determine probability of on-time completion with added epsi
  prob.ontime.added <- CalculateCompletionProbability(project.add, deadline)
  
  # Determine probability of on-time completion with added epsi
  prob.ontime.less <- CalculateCompletionProbability(project.less, deadline)
  
  # difference
  first.dif <- (prob.ontime.less - prob.ontime.added)/(2*epsi)
  first.dif[first.dif < 0] <- 0
  
  # first difference divided by standard deviation of path
  marginal.complete.prob <- first.dif #/sqrt(project$properties$variance[[path.idx]])
  
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



Crash.BruteForce <- function(activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
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
    activity.effect <- unlist(lapply(potential.crash.activities, function(i) ActivityCrashEffect(i, deadline, path.matrix, activities, crash.effect.var)))
    
    # crashing which activity will result in highest on-time completion
    crash.activity <- potential.crash.activities[which.max(activity.effect)]
    
    # crash the activity (reduce it's mean completion time by 50% or what's left in our crash budget)
    crash.reduction <- min(activities$means[crash.activity], crashing.reserve)
    # calculate the activities variances resulting from the crash
    activities$variances[crash.activity] <- (activities$cov[crash.activity] * (activities$means[crash.activity] - crash.reduction*crash.effect.var))^2
    # reduce the mean by the crash
    activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
    
    # this reduces our crashing reserve
    crashing.reserve <- crashing.reserve - crash.reduction
    
    # add crashed path to list, cannot re-crash a path
    crashed.activities <- c(crashed.activities, crash.activity)
  }
  
  return(activities)
}


ActivityCrashEffect <- function(i, deadline, path.matrix, activities, crash.effect.var){
  
  # crash the activity
  crash.activity <- i
  crash.reduction <- activities$means[i]
  # calculate the activities variances resulting from the crash
  activities$variances[crash.activity] <- (activities$cov[crash.activity] * (activities$means[crash.activity] - crash.reduction*crash.effect.var))^2
  # reduce the mean by the crash
  activities$means[crash.activity] <- activities$means[crash.activity] - crash.reduction
  
  # recalculate the project information
  project <- CalculatePathProperties(path.matrix, activities)
  
  # Determine probability of on-time completion
  prob <- CalculateCompletionProbability(project, deadline)
  
  return(prob)
}
