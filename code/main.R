cost <- 0.02
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
  deadline <- list(critical.path.percentile = 0.9, coeff.variation = seq(0,4, 0.2))
  crashing.reserve <- 1
  crash.effect.var <- F #if T, the variance increases due to crashing
  number.cores <- length(deadline$coeff.variation) # given there are 11 CV's iterated over
  sim_num <- 1000
  
  # import libraries
  library(data.table)
  library(mvtnorm)
  library(pbapply)
  library(ggplot2)
  library(ggthemes)
  library(RColorBrewer)
  library(directlabels)
  # library(matrixcalc)
  
  
  # import functions
  source('code/generate_projects_B.R')
  source('code/crashing_strategies_reward.R')
  
  # init results dataframe
  df <- NULL
  # loop activity and precedence density
  for (i in seq(1,sim_num)){
    print(i)
    activity.num <- sample(5:10,1)
    link.prob <- runif(1, 0.05, 0.4)
    # Generate random project3
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
    strategy.names <- c('PERT','modified.PERT','enumeration') #' 's.2b', 's.4a', 's.4b', 's.5'
    
    # loop through the deadline's coefficient of variation, calculate probability of on-time completion
    results <- pblapply(deadline$coeff.variation, function(CoV) VaryCoefVariation(CoV, deadline, strategy.names, path.matrix, activities, link.prob, crashing.reserve, crash.effect.var), cl = number.cores)
    df <- rbind(df, rbindlist(results, use.names=TRUE, fill=FALSE))
    # add other information
  }
  
  # save results
  variate.str <- 'CoV'
  fwrite(df, file = paste0('data/results/crashstrategies_vary_',variate.str,'.csv'))
  df <- fread(paste0('data/results/crashstrategies_vary_',variate.str,'.csv'))
  # plot the results
  # PlotMeanRatio(df, variate.str, strategy.names, log_axis = F)
  # PlotScheduleRisk(df, variate.str, strategy.names, log_axis = F)
  PlotRibbon(df, variate.str, strategy.names, log_axis = F)
  PlotMean(df, variate.str, strategy.names, log_axis = F)
  # PlotSRRatio(df, variate.str, strategy.names, log_axis = F)
  PlotDifference(df, variate.str, strategy.names, log_axis = F)
  # Plot_byVarianceRatio(df, variate.str, strategy.names, log_axis = F)
  
  # plot the forecast bias caused by not considering deadline uncertainty
  # PlotForecastBias(df, sim_num, strategy.names, variate.str)
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
  
  # Calculate path properties
  path.num <- nrow(path.matrix)
  path.properties <- CalculatePathMeanVariance(path.matrix, activities, path.num, deadline)
  # project.properties <- CalculateProjectProperties(path.matrix, activities)
  
  # return the following:
  # coefficient of variation | number of activities | number of paths | link density | probability of completion time for each crash strategy
  results <- list(CoV = CoV, num_acts = ncol(path.matrix), num_paths = nrow(path.matrix), precedence_density = link.prob, 
                  DV = deadline$variance, PM = mean(unlist(path.properties$mean)), PV = mean(unlist(path.properties$variance)), AV = mean(activities$variances), ACV = mean(activities$cov) )
  results <- c(results, prob.ontime)
  
  return(results)
}


CalcCrashCompletionProb <- function(crash.strategy, activities, crashing.reserve, path.matrix, deadline, crash.effect.var){
  # calculate the probability of on-time completion using a crashing strategy
  # crash the activity and update the activities list
  if (crash.strategy == 's.1'){
    activities <- Crash.None(activities)
  } else if (crash.strategy == 's.2a') {
    activities <- Crash.CriticalPath.2a(activities, crashing.reserve, path.matrix, deadline, crash.effect.var)
  } else if (crash.strategy == 'PERT') {
    activities <- Crash.CriticalPath.2b(activities, crashing.reserve, path.matrix, deadline, crash.effect.var)
  } else if (crash.strategy == 's.3a'){
    ret <- Crash.Heuristic(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance=F)
    # ret <- Crash.Optimal(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance = F)
    crashed.acts = ret$crashed
    activities = ret$acts
  } else if (crash.strategy == 'Heuristic'){
    ret <- Crash.Heuristic(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance=T)
    # ret <- Crash.Optimal(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance = T)
    crashed.acts = ret$crashed
    activities = ret$acts
  } else if (crash.strategy == 's.4a'){
    ret <- Crash.BruteForce(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance=F)
    crashed.acts = ret$crashed
    activities = ret$acts
  } else if (crash.strategy == 'Enumeration'){
    ret <- Crash.BruteForce(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance=T)
    # ret <- Crash.criticalOrder(activities, crashing.reserve, path.matrix, deadline, crash.effect.var, deadline.variance = T)
    crashed.acts = ret$crashed
    activities = ret$acts
  } else if (crash.strategy == 's.5'){
    activities$means <- activities$means * 0
  } 
  
  # update the path properties
  project <- CalculatePathProperties(path.matrix, activities)
  
  # Determine probability of on-time completion
  prob.ontime <- CalculateCompletionProbability(project, deadline, deadline.variance = T) 
  
  return(1-prob.ontime)
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


CalculatePathMeanVariance <- function(path.matrix, activities, path.num, deadline=NULL){
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
  # z score
  if (!missing(deadline)){
    path.properties$z <- (unlist(path.properties$mean) + deadline$mean)/sqrt(unlist(path.properties$variance) + deadline$variance)
  }
  return(path.properties)
}


CalculateCompletionProbability <- function(project, deadline, deadline.variance = TRUE){
  # Calculate the joint multivariate normal probability of the project finishing on-time
  #
  # Return:
  #   probability of on-time completion
  
  # calculate the difference between the paths and the deadlien
  
  # covariance 
  #   - covariance is sum of variance of shared actvities. the deadline is now considered a "shared activity"
  #     so can add it to the covariance
  if (deadline.variance){
    difference.covariance <- project$covariance + deadline$variance
  } else {
    difference.covariance <- project$covariance
  }
  
  # path means
  difference.mean <- unlist(project$properties[,mean]) - deadline$mean
  
  # remove paths that have no time
  difference.mean <- difference.mean[apply(difference.covariance,1,function(difference.covariance) !all(difference.covariance==0))] 
  difference.covariance <- difference.covariance[apply(difference.covariance,1,function(difference.covariance) !all(difference.covariance==0)),apply(difference.covariance,2,function(difference.covariance) !all(difference.covariance==0))] 
  
  # multivariate normal probability
  prob.ontime <- pmvnorm(upper = 0, mean = difference.mean,   sigma = difference.covariance, algorithm = GenzBretz(abseps = 1e6))[1]
  
  return(prob.ontime)
}


CalculateProjectProperties <- function(path.matrix, activities){
  # calculate some properties of paths 
  # criticality
  # expected completion time
  
  project.properties <- data.table()
  path.num <- nrow(path.matrix)
  
  # expected completion time
  path.mean <- lapply(seq(path.num), function(i) sum(activities$means[which(path.matrix[i,]==1)]))
  project.properties$expected.completion <- max(path.mean)
  
  # criticality
  
  
  
  return(project.properties)
}


Derivative <- function(path.idx, epsi, deadline, project, ignore.deadline){
  
  # difference mean and sd
  path.mean.dif <- project$properties$mean[[path.idx]] - deadline$mean
  if (ignore.deadline){
    path.sd <- sqrt(project$properties$variance[[path.idx]])
  } else {
    path.sd <- sqrt(project$properties$variance[[path.idx]] + deadline$variance)
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


###
# Plotting
###

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
      stat_summary(data = df, aes_string(x = variate.str, y = strategy), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95)) + 
      # mean
      stat_summary(data = df, aes_string(x = variate.str, y = strategy, color = shQuote(strategy.names[[i]]), linetype = shQuote(strategy.names[[i]])), geom="line", fun.y=median, size = 1)
    # geom_dl(aes_string(x = variate.str, y = strategy),  method = list(dl.combine("first.points", "last.points"), cex = 0.8))
    # iterate
    i <- i + 1
  }
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Schedule Risk') +
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
  
  # calculate the variance ratio
  # df$VR <- df$DV/df$PV
  # df$CoV <- df$VR # hack
  
  # reshape the df 
  vars.important <- c("CoV", "num_acts",'num_paths','precedence_density')
  df <- df[,c(vars.important, strategy.names),with=FALSE]
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
  plt <- ggplot(df3, aes(x = CoV, y = prob, group = strategy, colour = strategy, linetype = strategy)) + theme_hc()+
    # geom_point() +
    stat_summary(geom="line", fun.y=median, size = 1) 
    # geom_smooth() + 
    # scale_colour_discrete(guide = 'none') +
    # scale_x_discrete(expand=c(0, 1)) +
    # geom_dl(aes(label = strategy), method = list(dl.combine( "last.points"), cex = 0.8), size = 1) 
  
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Schedule Risk') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) + #, expand=c(0, 0.2)) +
    # scale_x_continuous(breaks = seq(0,100,by=1)) + 
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    # expand_limits(0,1.2) + 
    theme(
      legend.position = "bottom",
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


PlotMeanRatio <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  # library(reshape)
  vars.important <- c("CoV", "num_acts",'num_paths','precedence_density')
  df <- df[,c(vars.important, strategy.names),with=FALSE]
  
  # compute fraction of limiting
  df$s.1 <- NULL #df$s.1/df$s.4b
  df$s.2b <- df$s.2b/df$s.4b
  # df$s.2b <- df$s.2b/df$s.4b
  df$s.3a <- df$s.3a/df$s.4b
  df$s.3b <- df$s.3b/df$s.4b
  df$s.4a <- df$s.4a/df$s.4b
  df$s.5 <- NULL# df$s.5/df$s.4b
  # df$s.6a <- df$s.6a/df$s.4b
  # df$s.6b <- df$s.6b/df$s.4b
  df$s.4b <- df$s.4b/df$s.4b
  
  
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
    ylab('Probability of Completion Before the Deadline \n as Fraction of Optimal by Enumeration') +
    scale_y_continuous(breaks = seq(0,1, by = 0.01)) +
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
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_', variate.str, '-vary_mean_ratio.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}


PlotScheduleRisk <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the probability of running late
  
  # calculate the variance ratio
  # df$VR <- df$DV/df$PV
  # df$CoV <- df$VR # hack
  
  # reshape the df 
  vars.important <- c("CoV", "num_acts",'num_paths','precedence_density')
  df <- df[,c(vars.important, strategy.names),with=FALSE]
  df2 <- melt(df, id.vars=c("CoV", "num_acts",'num_paths','precedence_density'), variable.name = "strategy", value.name = "prob")
  # aggregate(df2$prob, by = list(df2[,c('CoV','strategy')]), mean)
  # mydt.mean <- df2[,lapply(.SD,mean,na.rm=TRUE),by=c(Cov, strategy),.SDcols=colstoavg]
  df3 <- df2[, lapply(.SD, mean), by = .(CoV, strategy)]
  
  # schedule risk is 1-p(running on time)
  df3$risk = 1-df3$prob
  
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  # plot  
  plt <- ggplot(df3, aes(x = CoV, y = risk, group = strategy, colour = strategy)) + theme_hc()+
    # geom_point() +
    stat_summary(geom="line", fun.y=median, size = 1) +
    # geom_smooth() + 
    # scale_colour_discrete(guide = 'none') +
    # scale_x_discrete(expand=c(0, 1)) +
    geom_dl(aes(label = strategy), method = list(dl.combine( "last.points"), cex = 0.8), size = 1) 
  
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Schedule Risk') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    scale_x_continuous(breaks = unique(df[[variate.str]])) + #, expand=c(0, 0.2)) +
    # scale_x_continuous(breaks = seq(0,100,by=1)) + 
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
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_', variate.str, '-vary_schedulerisk-mean.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}


PlotSRRatio <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  # library(reshape)
  vars.important <- c("CoV", "num_acts",'num_paths','precedence_density')
  df <- df[,c(vars.important, strategy.names),with=FALSE]
  
  # calculate schedule risk
  df[,(strategy.names) := lapply(.SD, function(x) 1-x), .SDcols=strategy.names]
  # compute fraction of PERT
  # df$s.1 <- NULL #df$s.1/df$s.4b
  
  # df$s.2b <- df$s.2b/df$s.4b
  df$s.3a <- df$s.3a/df$s.2b
  df$s.3b <- df$s.3b/df$s.2b
  df$s.4a <- df$s.4a/df$s.2b
  df$s.5 <- NULL# df$s.5/df$s.4b
  # df$s.6a <- df$s.6a/df$s.4b
  # df$s.6b <- df$s.6b/df$s.4b
  df$s.4b <- df$s.4b/df$s.2b
  df$s.2b <- df$s.2b/df$s.2b
  
  strategy.names <- c('s.2a', 's.3a', 's.3b', 's.4a', 's.4b')
  df[is.na(df)] <- 1
  df[df>1] <-  1#:= lapply(.SD, function(x) if(x>1){1}), .SDcols=strategy.names]
  
  
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
    ylab("Schedule Risk as a Ratio to PERT's Schedule Risk") +
    scale_y_continuous(breaks = seq(0,1, by = 0.05)) +
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
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_', variate.str, '-vary_schedulerisk_ratio.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}


PlotDifference <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  
  # reshape the df 
  df$heur_crit <- df$PERT - df$modified.PERT
  df$heur_brute <- df$modified.PERT - df$enumeration
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
    stat_summary(data = df, aes_string(x = variate.str, y = strategy), geom="line", fun.y=mean, size = 1)
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab(bquote(atop(Delta ~ plain("Schedule Risk:"), 'Modified-PERT vs PERT'))) +
    # coord_cartesian(ylim=c(-0.1,0.1)) +
    # scale_y_continuous(breaks = seq(-10,10, by = 0.05)) +
    # scale_x_continuous(breaks = unique(df[[variate.str]])) +
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
    stat_summary(data = df, aes_string(x = variate.str, y = 'heur_brute'), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95)) + 
    # mean
    stat_summary(data = df, aes_string(x = variate.str, y = 'heur_brute', color = shQuote(strategy.names[[i]]), linetype = shQuote(strategy.names[[i]])), geom="line", fun.y=mean, size = 1)
  
  
  # plotting
  plt <- plt +
    xlab("Deadline's Coefficient of Variation") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab(bquote(atop(Delta ~ plain("Schedule Risk:"), 'Enumeration vs Modified-PERT'))) +
    # ylim(-0.1,0.1) + 
    coord_cartesian(ylim=c(-0.1,0.1)) +
    # scale_y_continuous(breaks = seq(-0.1,0.1, by = 0.1)) +
    # scale_x_continuous(breaks = unique(df[[variate.str]])) +
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
  
  df$heur_crit <- df$s.4b - df$s.3a
  # df$heur_brute <- df$enumerate - df$critical.path
  df$id <- rep(seq(1,dim(df)[1]/length(unique(df$CoV))),each = length(unique(df$CoV)))
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  # ggplot(df, aes(CoV, num_paths, z = heur_crit, color = heur_crit)) + stat_density2d()
  
  strategy = 'heur_crit'
  variate.str = 'CoV'
  i = 1
  ggplot(df)+#[df$num_paths == 6]) + 
    # ribbon
    stat_summary(aes_string(x = variate.str, y = strategy), alpha = 0.2, fill = cols[[i]], geom="ribbon", fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95)) +
    # lines
    geom_line(aes_string(x = variate.str, y = strategy, group = 'id', color = 'num_paths'), size = 1, alpha = 0.5) +
    # median
    stat_summary(aes_string(x = variate.str, y = strategy), geom="line", fun.y=median, size = 1)+
    geom_dl(aes(x=CoV, y=heur_crit,label = id), method = list(dl.combine( "last.points"), cex = 0.8), size = 1) 
  
  
}


PlotForecastBias <- function(df, sim_num, strategy.names, variate.str){
  # calculate forecast bias
  # for 
  # get the first prob for each of the sims
  entry <- dim(df)[1]/sim_num
  zero_vals <- df[seq(1,dim(df)[1],entry),'s.1'][[1]]
  zero_list <- rep(zero_vals, each=entry)
  
  df$bias = zero_list - df$s.1
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
    coord_cartesian(ylim=c(0,0.3)) +
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


Plot_byVarianceRatio <- function(df, variate.str, strategy.names, log_axis = F){
  # plot the ribbon results
  
  # calculate the variance ratio
  df$VR <- df$DV/df$PV
  # df$CoV <- df$VR # hack
  
  # reshape the df 
  vars.important <- c("VR", "num_acts",'num_paths','precedence_density')
  df <- df[,c(vars.important, strategy.names),with=FALSE]
  df2 <- melt(df, id.vars=c("VR", "num_acts",'num_paths','precedence_density'), variable.name = "strategy", value.name = "prob")
  # aggregate(df2$prob, by = list(df2[,c('CoV','strategy')]), mean)
  # mydt.mean <- df2[,lapply(.SD,mean,na.rm=TRUE),by=c(Cov, strategy),.SDcols=colstoavg]
  df3 <- df2[, lapply(.SD, mean), by = .(VR, strategy)]
  
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # colors
  cols <- cbbPalette[1:length(strategy.names)]
  names(cols) <- strategy.names
  # setNames(cols) <- strategy.names
  ltype <- seq(1, length(strategy.names))
  names(ltype) <- strategy.names
  
  # plot  
  plt <- ggplot(df3, aes(x = VR, y = prob, group = strategy, colour = strategy)) + theme_hc()+
    # geom_point() +
    # stat_summary(geom="line", fun.y=median, size = 1) +
    geom_smooth() +
    # scale_colour_discrete(guide = 'none') +
    # scale_x_discrete(expand=c(0, 1)) +
    geom_dl(aes(label = strategy), method = list(dl.combine( "last.points"), cex = 0.8), size = 1) 
  
  
  # plotting
  plt <- plt +
    xlab("Ratio of Deadline Variane to \n Weighted Mean of Path Variance") +  #Standard deviation for path completion time') + #Mean path completion time') + #Standard deviation between \n mean path completion times ') + #Inter-path correlation') + #
    ylab('Probability of Completion \n Before the Deadline') +
    scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
    # scale_x_continuous(breaks = unique(df[[variate.str]]), expand=c(0, 0.2)) +
    # scale_x_continuous(breaks = seq(0,100,by=1)) + 
    scale_colour_manual(name="",values=cols) +
    scale_linetype_manual(name = '', values= rep(ltype)) +
    # expand_limits(0,1.2) + 
    theme(
      legend.position = "none",
      legend.key.width = unit(1,"cm"),
      text = element_text(size = 10)
    )    + 
    coord_cartesian(xlim=c(0, 25))
  # plot
  print(plt)
  
  # save figure
  fig.width <- 7  #3.25 # inches
  fig.height <- fig.width * (sqrt(5)-1)/2 + 0.67
  ggsave(paste0('fig/', format(Sys.time(), "%Y-%m-%d"),'_', variate.str, '-variance_ratio.pdf'), plot = plt, device = 'pdf', width =  fig.width, height = fig.height , unit = 'in' , dpi = 500) 
  
}
