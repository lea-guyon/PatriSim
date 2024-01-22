#!/usr/bin/env Rscript
require(ape)

# Get arguments ie working directory
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

dir=args[1]
setwd(dir)

#' get.intervals
#' 
#' Scan through BEAST tree log line-by-line and parse each 
#' Newick string to retrieve the coalescent intervals
#' 
#' @param treefile: (character) relative or absolute path to trees log file
#' @return List containing objects of class coalescentIntervals (ape) keyed
#'         by state (MCMC step number)
get.intervals <- function(treefile) {
  result <- list()
  
  con <- file(treefile, open='r')
  while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
    if (grepl("^tree", line)) {
      state <- gsub("^.+STATE_([0-9]+) .+$", "\\1", line)
      
      #nwkstr <- strsplit(line, " ")[[1]][4]  # does not support BEAST1
      nwkstr <- gsub("^[^(]+\\((.+);.*$", "(\\1;", line)
      phy <- read.tree(text=nwkstr)
      
      intervals <- coalescent.intervals(phy)
      result[[state]] <- intervals
    }
  }
  close(con)
  result
}


#' get.skyline
#' 
#' Generate coordinates of Bayesian skyline plot by averaging the 
#' effective population size estimates over time.
#' 
#' @param logfile: path to BEAST trace log file
#' @param treefile: path to BEAST tree log file
#' @param p.burnin: proportion of chain sample to discard, starting from 0
#' @param thin: thinning interval, defaults to 1 (retain all post-burnin samples)
#' @param bin.count: number of bins to calculate skyline trends
#' 
#' @return S3 object of class 'skyline'
get.skyline <- function(logfile, treefile, p.burnin=0.1, thin=1, bin.count=100) {  
  # load log file
  log <- read.table(logfile, sep='\t', comment.char='#', header=T)
  
  # parse tree file
  intervals <- get.intervals(treefile)
  
  # determine if this is a BEAST1 or BEAST2 logfile
  if (is.element("treeModel.rootHeight", names(log))) {
    t.heights <- log$treeModel.rootHeight
  } else if (is.element("TreeHeight", names(log))) {
    t.heights <- log$TreeHeight
  } else {
    stop("Failed to locate tree height parameter in log file")
  }
  
  # maximum time is the root height's mean/median
  t.mean <- mean(t.heights)
  q <- quantile(t.heights, c(0.025, 0.5, 0.975))
  t.lower <- q[1]  # <- default
  t.median <- q[2]
  t.upper <- q[3]
  max.height <- t.lower
  
  age.youngest <- 0  # plot backwards in time
  min.time <- 0
  max.time <- max.height - age.youngest
  delta <- (max.time - min.time) / (bin.count-1)
  
  group.size.count <- sum(grepl("groupSize", names(log)))
  
  # store age at the end of each group
  group.times <- matrix(NA, nrow=nrow(log), ncol=group.size.count)
  
  # determine time points of population size transitions for every sample
  n.states <- nrow(log)
  if (is.element("Sample", names(log))) {
    states <- log$Sample
  } else if (is.element("state", names(log))) {
    states <- log$state
  } else {
    stop("Failed to locate chain step number in log")
  }
  
  thinned <- seq(round(n.states*p.burnin), n.states, thin)
  for (i in thinned) {
    state <- as.character(states[i])
    group.sizes <- log[i, grepl("groupSize", names(log))]
    
    # get cumulative times
    cumul.t <- cumsum(intervals[[state]]$interval.length)
    group.times[i, ] <- cumul.t[ cumsum(as.integer(group.sizes)) ]
  }
  
  # assign population sizes to bins
  heights <- seq(0, max.height, length.out=bin.count)
  bins <- matrix(nrow=nrow(log), ncol=bin.count)
  for (i in thinned) {
    state <- as.character(states[i])
    pop.sizes <- log[i, grepl("popSize", names(log))]
    index <- as.integer(cut(heights, c(0, group.times[i,]), right=F))
    bins[i,] <- unlist(sapply(index, function(i) {
      ifelse (is.na(i), NA, pop.sizes[i])
    }))
  }
  
  # remove burnin rows before returning
  bins <- na.omit(bins)
  skyline <- list(
    raw.data=bins, 
    mean=apply(bins, 2, mean),
    median=apply(bins, 2, median),
    lower=apply(bins, 2, function(x) quantile(x, 0.025)),
    upper=apply(bins, 2, function(x) quantile(x, 0.975)),
    time=heights
  )
  class(skyline) <- 'skyline'
  skyline
}


#' plot
#' 
#' Generic plot method for objects of class 'skyline'
#' Mito
logfile <- 'Sim_Mito.log'
treefile <- "Sim_Mito.trees"
sky <- get.skyline(logfile, treefile)
sky_df <- data.frame("time" = sky$time, "mean" = sky$mean, "median" = sky$median, 
                     "upper" = sky$upper, "lower" = sky$lower)
write.table(sky_df, paste0(dir, "/skyline_Mito"), col.names = T, row.names = F)

#' Y
logfile <- 'Sim_Y.log'
treefile <- "Sim_Y.trees"
sky <- get.skyline(logfile, treefile)
sky_df <- data.frame("time" = sky$time, "mean" = sky$mean, "median" = sky$median, 
                     "upper" = sky$upper, "lower" = sky$lower)
write.table(sky_df, paste0(dir, "/skyline_Y"), col.names = T, row.names = F)

#plot(sky)
plot.skyline <- function(obj, ylim=NA, 
                         poly.border=NA, poly.col='lightblue',
                         line.col='black', line.lwd=2,
                         xlab='Time', ylab='Effective population size', 
                         ...) {
  if (is.na(ylim)) {
    ylim <- range(obj$raw.data)  
  }
  
  plot(obj$time, obj$mean, type='n', ylim=ylim, log='y', 
       xlab=xlab, ylab=ylab, ...)
  polygon(x=c(obj$time, rev(obj$time)), 
          y=c(obj$upper, rev(obj$lower)), 
          border=poly.border, col=poly.col)
  lines(obj$time, obj$mean, col=line.col, lwd=line.lwd)
}
#plot.skyline(sky)
