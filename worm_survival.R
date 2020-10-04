#!/usr/bin/env Rscript

main <- function() {

  ## load libraries
  packages = c('reshape2','tidyr','survival','survminer','dplyr','MASS','stringr')
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, then load them into the R session.
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  msg.trap <- capture.output( suppressMessages( check.packages(packages)) ) # silent mode
  
  ## parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  ## data file (argument #1)
  file <- args[1]  ## needed for importing data
  ## graph type option (argument #2)
  if (is.na(args[2])==FALSE){
    gextension = args[2]
  } else {
    gextension = 'jpeg'
  }
  ## graph height option (argument #3)
  if (is.na(args[3])==FALSE){
    grheight = args[3]
  } else {
    grheight = 3
  }
  ## graph width option (argument #4)
  if (is.na(args[4])==FALSE){
    grwidth = args[4]
  } else {
    grwidth = 7
  }
  
  
  #****************** MAIN SCRIPT *******************
  
  ## Read data from usual dataset (day, condition1, condition2, etc.)
  worms <- read.csv(file, header = FALSE)
  
  ## store important info as variables
  time_type = worms[1,1]
  time_interval_v = as.numeric(as.vector(worms[[1]])[-1])
  starting_worms_v = as.numeric(worms[2,][-1])
  conditions_v = as.character(worms[1,][-1])
  study_length = as.numeric(nrow(worms))-2
  total_worms = sum(starting_worms_v)
  
  ## format intermediate data frame for downstream use
  worms.t <- data.frame(t(worms[-1])) # transpose original dataframe
  colnames(worms.t) = worms[, 1] # rename columns
  names(worms.t)[1] = "condition"
  worms.melted <- melt(worms.t, id.vars=c('condition')) # melt dataframe
  names(worms.melted)[2] = time_type
  names(worms.melted)[3] = "live_count"
  subtr = rep(starting_worms_v, times = study_length+1)
  worms.melted$subtr = subtr
  worms.melted$dead_count = worms.melted$subtr - as.numeric(worms.melted$live_count)
  worms.melted <- with(worms.melted,  worms.melted[order(worms.melted$condition) , ])
  row.names(worms.melted) <- 1:nrow(worms.melted) # renumber rows after reordering
   
  ## Calculate fustat & futime
  ##   fustat = 1 if the worm died before the end of the study
  ##   fustat = 0 if the worm survived the whole study
  ##   futime = time worm died (if fustat = 1), or last day of study (if fustat = 0)
  
  futime = as.vector(c())
  fustat = as.vector(c())
  worms.melted = na.omit(worms.melted)
  for (k in seq(from=0, to=as.numeric(length(conditions_v)*study_length), 
                by=study_length+1)) {
    daily_death_v = as.vector(c())
    for (i in 1:study_length){
      daily_death = as.vector(worms.melted$dead_count[k+i+1] - worms.melted$dead_count[k+i])
      daily_death_v = append(daily_death_v, daily_death)
    }
    ending_dead_worms = sum(daily_death_v)
    ending_live_worms = as.numeric(worms.melted$live_count[k+1]) - ending_dead_worms
    fustat = append(fustat, rep(1, times=ending_dead_worms))
    fustat = append(fustat, rep(0, times=ending_live_worms))
    futime = append(futime, rep(time_interval_v[-1], times=daily_death_v))
    futime = append(futime, rep(study_length, times=ending_live_worms))
  }
  
  ## Build new data frame for stats
  condition = sort(rep(as.vector(worms.t[[1]]), times = starting_worms_v))
  worm.df = data.frame(condition)
  worm.df$fustat = fustat
  worm.df$futime = futime
  
  ## Create new formatted survival dataset and output csv
  write.csv(worm.df, file="FormattedSurvivalData.csv", row.names = TRUE)
  
  ## Fit survival data using the Kaplan-Meier method
  time_gap = time_interval_v[2] # determines breaks on x-axis
  surv_object <- Surv(time = worm.df$futime, event = worm.df$fustat)
  fit1 <- do.call(survfit, list(surv_object ~ condition, data = worm.df))
  ggsurvplot(fit1, data = worm.df, xlim = c(0,study_length), 
                             break.time.by = time_gap, 
                             xlab = paste("Time (",time_type,")",sep=""), legend="right") 
  ggsave(filename=paste("survival_plot",".",gextension,sep=""), device=gextension, 
         height=as.numeric(grheight), width=as.numeric(grwidth))
  
  # Fit a Cox proportional hazards model (for the hell of it)
  fit.coxph <- coxph(surv_object ~ condition, data = worm.df)
  ggforest(fit.coxph, data = worm.df)
  ggsave(filename=paste("cox_plot",".",gextension,sep=""), device=gextension, 
         height=as.numeric(grheight), width=as.numeric(grwidth))
  
  ## Log Rank tests
  cmbo = data.frame(t(combn(conditions_v,2)))
  num_comb = as.numeric(nrow(cmbo))
  surv_pval_v = as.vector(c())
  for (i in 1:num_comb){
    subset_pair = worm.df[which(worm.df$condition==cmbo[i,1]
                       |worm.df$condition==cmbo[i,2]),]
    surv_object <- Surv(time = subset_pair$futime, event = subset_pair$fustat)
    fit1 <- do.call(survfit, list(surv_object ~ condition, data = subset_pair))
    surv_curve <- ggsurvplot(fit1, pval = TRUE, data = subset_pair, xlim = c(0,study_length), 
                          break.time.by = time_gap, 
                          xlab = paste("Time (",time_type,")",sep=""), legend="right")  
    surv_pval = surv_pvalue(fit1, data=subset_pair)$pval
    surv_pval_v = append(surv_pval_v, surv_pval)
  }
  cmbo$pval = surv_pval_v
  colnames(cmbo) <- c('condition1','condition2','pvalue')
  
  ## Output csv for Log Rank Values
  write.csv(cmbo, file="LogRankPvalues.csv", row.names = FALSE)
  
  ## Calculate LT50
  lt50.df = data.frame(conditions_v)
  lt50_v = as.vector(c())
  day.list = c(1:(study_length+1)) #this makes day 0 = day 1, so subtract 1 from Dose, below
  for (i in 1:as.numeric(length(conditions_v))){
    subset = worms.melted[which(worms.melted$condition==conditions_v[i]),]
    subset$live_count = as.numeric(subset$live_count)
    y = cbind(subset$live_count,subset$dead_count)
    model.results = glm(y~day.list,binomial)
    lt50 = dose.p(model.results,p=0.5)[1]-1 # what time gives 50% survival
    lt50 = as.numeric(round(lt50, digits=1))
    lt50_v = append(lt50_v, lt50)
  }
  lt50.df$lt50 = lt50_v
  colnames(lt50.df) = c('condition','LT50')
  
  ## Output csv for LT50 (Lethal Time)
  write.csv(lt50.df, file="LT50.csv", row.names = FALSE)

  print('Note that the file Rplots.pdf is empty and can be deleted.')
}
main()
