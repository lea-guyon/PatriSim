library(cowplot)
library(ggplot2)
library(scales)
library(boot)

mean_curve <- function(data) {
  m=c()
  for (t in unique(data$time)) {
    data_time = data[data$time == t,]
    m = c(m, mean(data_time$median, na.rm = T))
  }
  new_data = data.frame("time" = unique(data$time), "median"  = m)
  return(new_data)
}

skyline <- function(dir, xmin) {
  setwd(paste0(dir, "1/"))
  
  data_Y <- read.table("skyline_Y", header=T)
  data_Mito <- read.table("skyline_Mito", header=T)
  
  data_Y$time = - data_Y$time / 25
  data_Mito$time = - data_Mito$time / 25
  
  ## Demography
  
  p1 <- ggplot(data_Y, aes(time, median/25)) +
    annotate("rect", xmin = xmin, xmax = 0, ymin = 0, ymax = 8000,
             alpha = .1,fill = "orange") +
    annotate("text", x = -35, y = 4000, label = "mtDNA", col = 'brown2', size = 5, fontface = "bold") +
    annotate("text", x = -5, y = 200, label = "Y", col = 'darkgoldenrod2', size = 5, fontface = "bold") +
    geom_line(col = "darkgoldenrod2", alpha = 1/10) +
    geom_line(data = data_Mito, aes(time, median/25), col = "brown2", alpha = 1/10) +
    scale_y_continuous(trans = 'log10', expand = c(0,0), limits = c(30, NA)) +
    scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 0), n.breaks = 6) +
    labs(x="Time (generations)", y=expression(''*N[e]*'')) +
    annotate("segment", x = xmin, xend = xmin, y = 0, yend = 5000, colour = "salmon2", linetype='longdash', linewidth=0.5) +
    annotate("text", x = xmin, y = 6000, label = bquote(''*t[0]*''), col = 'salmon2', size = 6) +
    theme_cowplot()
  
  interpolated_Y = approx(x=data_Y$time, y=data_Y$median/25, xout = seq(xmin-30, 0, 10))
  interpolated_Mito = approx(x=data_Mito$time, y=data_Mito$median/25, xout = seq(xmin-30, 0, 10))
  
  data_MY <- data.frame("time" = interpolated_Y$x, "median" = interpolated_Mito$y/interpolated_Y$y)
  df_MY = data_MY
  
  df_Y = data.frame("time" = interpolated_Y$x, "median" = interpolated_Y$y, "Replicat" = rep(1, length(interpolated_Y$x)))
  df_Mito = data.frame("time" = interpolated_Mito$x, "median" = interpolated_Mito$y, "Replicat" = rep(1, length(interpolated_Mito$x)))
  
  p2 <- ggplot(data_MY, aes(time, median)) +
    annotate("rect", xmin = xmin, xmax = 0, ymin = 0, ymax = 30,
             alpha = .1,fill = "orange") +
    geom_line(col = "darkblue", alpha = 0.05) +
    labs(x="Time (generations)", y=expression(''*N[e]^(mt/Y)*'')) +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 0), n.breaks = 6) +
    geom_hline(yintercept = 17, col = "darkgray", linetype="dashed") +
    annotate("segment", x = xmin, xend = xmin, y = 0, yend = 28, colour = "salmon2", linetype='longdash', linewidth=0.5) +
    annotate("text", x = xmin, y = 29, label = bquote(''*t[0]*''), col = 'salmon2', size = 6) +
    theme_cowplot()
  
  maxPeak = c(max(data_MY$median, na.rm=T))
  initY = c(df_Y$median[df_Y$time == xmin])
  minY = c(min(df_Y$median, na.rm=T))
  df_Y2 = na.omit(df_Y)
  timeMinY = df_Y2[df_Y2$median == min(df_Y2$median, na.rm = T),]$time
  Ne_N = c(750/min(df_Y$median, na.rm=T))
  
  for (i in seq(1, 200)) {
    dir2 = paste0(dir, i, "/")
    setwd(dir2)
    data_Y <- read.table("skyline_Y", header=T)
    data_Mito <- read.table("skyline_Mito", header=T)
    
    data_Y$time = - data_Y$time/25
    data_Mito$time = - data_Mito$time/25
    
    p1 <- p1 + geom_line(data = data_Y, col = "darkgoldenrod2", alpha = 1/10) +
      geom_line(data = data_Mito, aes(time, median/25), col = "brown2", alpha = 1/10)
    
    new_interpolated_Y = approx(x=data_Y$time, y=data_Y$median/25, xout = seq(xmin-30, 0, 10))
    new_interpolated_Mito = approx(x=data_Mito$time, y=data_Mito$median/25, xout = seq(xmin-30, 0, 10))
    data_MY <- data.frame("time" = new_interpolated_Y$x, "median" = new_interpolated_Mito$y/new_interpolated_Y$y)
    
    new_df_Y = data.frame("time" = new_interpolated_Y$x, "median" = new_interpolated_Y$y, "Replicat" = rep(i, length(interpolated_Y$x)))
    new_df_Mito = data.frame("time" = new_interpolated_Mito$x, "median" = new_interpolated_Mito$y, "Replicat" = rep(i, length(interpolated_Mito$x)))
    df_Y = rbind(df_Y, new_df_Y)
    df_Mito = rbind(df_Mito, new_df_Mito)
    df_MY = rbind(df_MY, data_MY)
    
    p2 <- p2 + geom_line(data = data_MY, aes(time, median), col = "darkblue", alpha = 0.05)
    maxPeak = c(maxPeak, max(data_MY$median, na.rm=T))
    initY = c(initY, new_df_Y$median[new_df_Y$time == xmin])
    minY = c(minY, min(new_df_Y$median, na.rm=T))
    df_ratio = data.frame("initY" = initY, "minY" = minY)
    df_Y2 = na.omit(df_Y)
    timeMinY = df_Y2[df_Y2$median == min(df_Y2$median, na.rm = T),]$time
    Ne_N = c(Ne_N, 750/min(new_df_Y$median, na.rm=T))
  }
  
  mean_Y = mean_curve(df_Y)
  mean_Mito = mean_curve(df_Mito)
  mean_MY = mean_curve(df_MY)
  
  # p1 <- p1 + 
  #   geom_smooth(data = df_Y, aes(time, median), col = "darkgoldenrod2", fill = "darkgoldenrod2", span = 0.2) + 
  #   geom_smooth(data = df_Mito, aes(time, median), col = "brown2", fill = "brown2", span = 0.2)
  
  # p1_bis <- p1 +
  #   geom_line(data=mean_Y, aes(time, median), col = "darkgoldenrod2", linewidth = 2) +
  #   geom_line(data=mean_Mito, aes(time, median), col = "brown2", linewidth = 2)
  # 
  # p2_bis <- p2 +
  #    geom_line(data=mean_MY, aes(time, median), col = "darkblue", linewidth = 2)
   
  # p3 <- plot_grid(p1_bis, p2_bis, labels = c("a", "b"), ncol = 2)
  p3 <- c()
  
  return(list(data_MY, mean_MY, p3, maxPeak, mean_Y, mean_Mito, df_ratio, Ne_N))
}

ratio_ic <- function(data, indices) {
  d <- data[indices,]
  mean_y = mean(d$initY, na.rm = T)/mean(d$minY, na.rm = T)
  return(mean_y)
}

boot_ratio <- function(df_ratio, ratio_ic, model) {
  bootstrap <- boot(data=df_ratio, statistic = ratio_ic, R = 10000)
  ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
  mean <- bootstrap$t0
  ICinf <- ci$bca[4]
  ICsup <- ci$bca[5]
  
  bottleneck <- data.frame("mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup, "model" = model)
  return(bottleneck)
}

mean_maxPeak <- function(data, indices) {
  d <- data[indices]
  meanMP = mean(d, na.rm=T)
  return(meanMP)
}

boot_maxPeak <- function(maxPeak, mean_maxPeak, model) {
  bootstrap <- boot(data=maxPeak, statistic = mean_maxPeak, R = 10000)
  ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
  mean <- bootstrap$t0
  ICinf <- ci$bca[4]
  ICsup <- ci$bca[5]
  
  mp <- data.frame("mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup, "model" = model)
  return(mp)
}

bsp <- function(paths, name, legend, xmin, yminY, ymaxY, yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace) {
  # Build legend
  tab_legend = data.frame("model" = name, "value" = rep(1, length(name)))
  ggp_split_1 <- ggplot(tab_legend,                    # Create ggplot2 plot of first subset
                        aes(model, value,
                            color = model,
                            shape = model,
                            linetype = model)) +
    geom_point(size = 5, stroke = 2) +
    geom_line(linewidth = 1) +
    scale_color_manual(labels = legend,
                       values = col) +
    scale_shape_manual(labels = legend,
                        values = sh) +
    scale_linetype_manual(labels = legend,
                          values = lt) +
    labs(color = "Variance", shape = "Variance", linetype = "Variance") +
    theme(legend.key.size = unit(2, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.title = element_text(size=22, face = "bold"),
          legend.text = element_text(size=19)) + 
    guides(color = guide_legend(ncol=nCol, bycol=TRUE))
  ggp_legend_split_1 <- get_legend(ggp_split_1)
  
  flag = F
  for (i in seq(1, length(paths))) {
    path = paths[i]
    colour = col[i]
    shape = sh[i]
    linetype = lt[i]
    model = name[i]
    
    sp <- skyline(path, xmin)
    data_MY <- sp[[1]]
    mean_MY <- sp[[2]]
    mean_Y <- sp[[5]]
    print(mean_Y)
    mean_mt <- sp[[6]]
    Ne_N <- sp[[8]]
    
    if (flag == F) {
      maxMY <- data.frame("max"= c(max(mean_MY$median)), "model" = c(model))
      reduction_factor <- data.frame("redFactor"= c(750/min(mean_Y$median)), "model" = c(model))
      print(reduction_factor)
      
      bspMY <- ggplot(data_MY, aes(time, median)) +
        annotate("rect", xmin = xmin, xmax = 0, ymin = yminMY, ymax = ymaxMY,
                 alpha = .1,fill = "orange") +
        labs(x="Time (generations)", y=expression(''*N[e]^(mt/Y)*'')) +
        scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(0, NA)) +
        scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 10), n.breaks = 6) +
        geom_hline(yintercept = 17, col = "darkgray", linetype="dashed") +
        annotate("segment", x = xmin, xend = xmin, y = 0, yend = segmentMY, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin, y = textMY, label = bquote(''*t[0]*''), col = 'salmon2', size = 8) +
        theme_cowplot() +
        theme(text=element_text(size=fontsize1),
              axis.text = element_text(size=fontsize2),
              panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line(color="bisque")) +
        geom_point(data=mean_MY, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_MY, aes(time, median), col = colour, linewidth = 1, linetype = linetype)
        
      
      bspY <- ggplot(mean_Y, aes(time, median)) +
        annotate("rect", xmin = xmin, xmax = 0, ymin = yminY, ymax = ymaxY,
                 alpha = .1,fill = "orange") +
        labs(x="Time (generations)", y=expression(''*N[e]^(Y)*'')) +
        scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(0, NA)) +
        scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 10), n.breaks = 6) +
        annotate("segment", x = xmin, xend = xmin, y = 0, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin, y = textY, label = bquote(''*t[0]*''), col = 'salmon2', size = 8) +
        theme_cowplot() +
        theme(text=element_text(size=fontsize1),
              axis.text = element_text(size=fontsize2),
              panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line(color="bisque")) +
        geom_point(data=mean_Y, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_Y, aes(time, median), col = colour, linewidth = 1, linetype = linetype)
        
      
      bspM <- ggplot(mean_mt, aes(time, median)) +
        annotate("rect", xmin = xmin, xmax = 0, ymin = yminM, ymax = ymaxM,
                 alpha = .1,fill = "orange") +
        labs(x="Time (generations)", y=expression(''*N[e]^(mt)*'')) +
        scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(0, NA)) +
        scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 10), n.breaks = 6) +
        annotate("segment", x = xmin, xend = xmin, y = 0, yend = segmentM, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin, y = textM, label = bquote(''*t[0]*''), col = 'salmon2', size = 8) +
        theme_cowplot() +
        theme(text=element_text(size=fontsize1),
              axis.text = element_text(size=fontsize2),
              panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line(color="bisque")) +
        geom_point(data=mean_mt, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_mt, aes(time, median), col = colour, linewidth = 1, linetype = linetype)
      
      bootNeN <- boot_maxPeak(Ne_N, mean_maxPeak, model)
      flag = T
    }
    else {
      maxMY = merge(maxMY, data.frame("max"= c(max(mean_MY$median)), "model" = c(model)), all = T)
      reduction_factor = merge(reduction_factor, data.frame("redFactor"= c(750/min(mean_Y$median)), "model" = c(model)), all=T)
      print(reduction_factor)
      
      bspMY <- bspMY +
        geom_point(data=mean_MY, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_MY, aes(time, median), col = colour, linewidth = 1, linetype = linetype) 
      bspY <- bspY +
        geom_point(data=mean_Y, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_Y, aes(time, median), col = colour, linewidth = 1, linetype = linetype) 
      bspM <- bspM +
        geom_point(data=mean_mt, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_mt, aes(time, median), col = colour, linewidth = 1, linetype = linetype) 
      
      bootNeN_bis <- boot_maxPeak(Ne_N, mean_maxPeak, model)
      bootNeN <- merge(bootNeN, bootNeN_bis, all = TRUE)
    }
  }
  p <- plot_grid(bspY + theme(legend.position="none"), 
                 bspM + theme(legend.position="none"),
                 bspMY + theme(legend.position="none"),
                 ncol=3, labels = "auto", label_size = fontsize1)
  
  p_bis <- plot_grid(p, legend=ggp_legend_split_1, ncol = 1, rel_heights = c(1-legendSpace, legendSpace))
  
  return(list(p_bis, maxMY, reduction_factor, bootNeN))
}

bsp_2_transitions <- function(paths, name, legend, xmin1, xmin2, yminY, ymaxY, yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace) {
  # Build legend
  tab_legend = data.frame("model" = name, "value" = rep(1, length(name)))
  ggp_split_1 <- ggplot(tab_legend,                    # Create ggplot2 plot of first subset
                        aes(model, value,
                            color = model,
                            shape = model,
                            linetype = model)) +
    geom_point(size = 5, stroke = 2) +
    geom_line(linewidth = 1) +
    scale_color_manual(labels = legend,
                       values = col) +
    scale_shape_manual(labels = legend,
                       values = sh) +
    scale_linetype_manual(labels = legend,
                          values = lt) +
    labs(color = "Scenarios", shape = "Scenarios", linetype = "Scenarios") +
    theme(legend.key.size = unit(2, "cm"),
          legend.key.height = unit(2, "cm"),
          legend.title = element_text(size=22, face = "bold"),
          legend.text = element_text(size=19)) + 
    guides(color = guide_legend(ncol=nCol, bycol=TRUE))
  ggp_legend_split_1 <- get_legend(ggp_split_1)
  
  flag = F
  for (i in seq(1, length(paths))) {
    path = paths[i]
    colour = col[i]
    shape = sh[i]
    linetype = lt[i]
    model = name[i]
    
    sp <- skyline(path, xmin)
    data_MY <- sp[[1]]
    mean_MY <- sp[[2]]
    mean_Y <- sp[[5]]
    mean_mt <- sp[[6]]
    Ne_N <- sp[[8]]
    
    if (flag == F) {
      maxMY <- data.frame("max"= c(max(mean_MY$median)), "model" = c(model))
      reduction_factor <- data.frame("redFactor"= c(750/min(mean_Y$median)), "model" = c(model))
      print(reduction_factor)
      
      bspMY <- ggplot(data_MY, aes(time, median)) +
        annotate("rect", xmin = xmin1, xmax = xmin2, ymin = yminMY, ymax = ymaxMY,
                 alpha = .1,fill = "orange") +
        annotate("rect", xmin = xmin2, xmax = 0, ymin = yminMY, ymax = ymaxMY,
                 alpha = .05,fill = "orange") +
        labs(x="Time (generations)", y=expression(''*N[e]^(mt/Y)*'')) +
        scale_y_continuous(expand = c(0,0), limits = c(0, NA), n.breaks = 6) +
        scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 10), n.breaks = 6) +
        geom_hline(yintercept = 17, col = "darkgray", linetype="dashed") +
        annotate("segment", x = xmin, xend = xmin, y = 0, yend = segmentMY, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin, y = textMY, label = bquote(''*t[0]*''), col = 'salmon2', size = 8) +
        annotate("segment", x = xmin+100, xend = xmin+100, y = 0, yend = segmentMY, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin+100, y = textMY, label = bquote(''*t[1]*''), col = 'salmon2', size = 8) +
        theme_cowplot() +
        theme(text=element_text(size=fontsize1),
              axis.text = element_text(size=fontsize2),
              panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line(color="bisque")) +
        geom_point(data=mean_MY, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_MY, aes(time, median), col = colour, linewidth = 1, linetype = linetype)
      
      
      bspY <- ggplot(mean_Y, aes(time, median)) +
        annotate("rect", xmin = xmin1, xmax = xmin2, ymin = yminY, ymax = ymaxY,
                 alpha = .1,fill = "orange") +
        annotate("rect", xmin = xmin2, xmax = 0, ymin = yminY, ymax = ymaxY,
                 alpha = .05,fill = "orange") +
        labs(x="Time (generations)", y=expression(''*N[e]^(Y)*'')) +
        scale_y_continuous(expand = c(0,0), limits = c(0, NA), n.breaks = 6) +
        scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 10), n.breaks = 6) +
        annotate("segment", x = xmin, xend = xmin, y = 0, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin, y = textY, label = bquote(''*t[0]*''), col = 'salmon2', size = 8) +
        annotate("segment", x = xmin+100, xend = xmin+100, y = 0, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin+100, y = textY, label = bquote(''*t[1]*''), col = 'salmon2', size = 8) +
        theme_cowplot() +
        theme(text=element_text(size=fontsize1),
              axis.text = element_text(size=fontsize2),
              panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line(color="bisque")) +
        geom_point(data=mean_Y, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_Y, aes(time, median), col = colour, linewidth = 1, linetype = linetype)
      
      
      bspM <- ggplot(mean_mt, aes(time, median)) +
        annotate("rect", xmin = xmin1, xmax = xmin2, ymin = yminM, ymax = ymaxM,
                 alpha = .1,fill = "orange") +
        annotate("rect", xmin = xmin2, xmax = 0, ymin = yminM, ymax = ymaxM,
                 alpha = .05,fill = "orange") +
        labs(x="Time (generations)", y=expression(''*N[e]^(mt)*'')) +
        scale_y_continuous(expand = c(0,0), limits = c(0, NA), n.breaks = 6) +
        scale_x_continuous(expand = c(0,0), limits = c(xmin-30, 10), n.breaks = 6) +
        annotate("segment", x = xmin, xend = xmin, y = 0, yend = segmentM, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin, y = textM, label = bquote(''*t[0]*''), col = 'salmon2', size = 8) +
        annotate("segment", x = xmin+100, xend = xmin+100, y = 0, yend = segmentM, colour = "salmon2", linetype='longdash', linewidth=1) +
        annotate("text", x = xmin+100, y = textM, label = bquote(''*t[1]*''), col = 'salmon2', size = 8) +
        theme_cowplot() +
        theme(text=element_text(size=fontsize1),
              axis.text = element_text(size=fontsize2),
              panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line(color="bisque")) +
        geom_point(data=mean_mt, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_mt, aes(time, median), col = colour, linewidth = 1, linetype = linetype)
      
      bootNeN <- boot_maxPeak(Ne_N, mean_maxPeak, model)
      flag = T
    }
    else {
      maxMY = merge(maxMY, data.frame("max"= c(max(mean_MY$median)), "model" = c(model)), all = T)
      reduction_factor = merge(reduction_factor, data.frame("redFactor"= c(750/min(mean_Y$median)), "model" = c(model)), all=T)
      print(reduction_factor)
      
      bspMY <- bspMY +
        geom_point(data=mean_MY, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_MY, aes(time, median), col = colour, linewidth = 1, linetype = linetype) 
      bspY <- bspY +
        geom_point(data=mean_Y, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_Y, aes(time, median), col = colour, linewidth = 1, linetype = linetype) 
      bspM <- bspM +
        geom_point(data=mean_mt, aes(time, median), col = colour, shape = shape, size = 5, stroke = 2) +
        geom_line(data=mean_mt, aes(time, median), col = colour, linewidth = 1, linetype = linetype) 
      
      bootNeN_bis <- boot_maxPeak(Ne_N, mean_maxPeak, model)
      bootNeN <- merge(bootNeN, bootNeN_bis, all = TRUE)
    }
  }
  p <- plot_grid(bspY + theme(legend.position="none"), 
                 bspM + theme(legend.position="none"),
                 bspMY + theme(legend.position="none"),
                 ncol=3, labels = "auto", label_size = fontsize1)
  
  p_bis <- plot_grid(p, legend=ggp_legend_split_1, ncol = 1, rel_heights = c(1-legendSpace, legendSpace))
  
  return(list(p_bis, maxMY, reduction_factor, bootNeN))
}

#setwd("/set/your/working/directory")

paths = c("/BEAST/bilateral/regular/r=0.01/patrilocal_villages/",
          "/BEAST/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_a/",
          "/BEAST/unilineal/regular/r=0.01/k=0/FT=150/e=0.15/patrilineal_villages_b/",
          "/BEAST/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_c/",
          "/BEAST/unilineal/regular/r=0.01/k=0/FT=150/e=0.15/patrilineal_villages_d/",
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_e/",
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/e=0.15/patrilineal_villages_f/",
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g/",
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/e=0.15/patrilineal_villages_h/")
name = c("1", "2a", "2b", "2c", "2d", "2e", "2f", "2g", "2h")
legend = c("1: bilateral descent", 
           "2a: patrilineal descent, random fission, no variance, no violence",
           "2b: patrilineal descent, random fission, no variance, violence",
           "2c: patrilineal descent, lineal fission, no variance, no violence",
           "2d: patrilineal descent, lineal fission, no variance, violence",
           "2e: patrilineal descent, random fission, variance, no violence",
           "2f: patrilineal descent, random fission, variance, violence",
           "2g: patrilineal descent, lineal fission, variance, no violence",
           "2h: patrilineal descent, lineal fission, variance, violence")
xmin = -100
yminY = 0
ymaxY = 2500
yminM = 0
ymaxM = 2600
yminMY = 0
ymaxMY = 20
segmentY = 2300
textY = 2400
segmentM = 2450
textM = 2530
segmentMY = 18
textMY = 19
col = c("#7bc84e", "#69aa63", "#55c6b2", "#5a81bd", "#4c69d5", "#bd634a", "#d5732d", "#a15cbe", "#c15787")
sh = c(8, 22, 23, 21, 24, 15, 18, 16, 17)
lt = c(rep("dashed", 5), rep("solid", 4))
fontsize1 = 22
fontsize2 = 18
nCol = 2
legendSpace = 0.4
       
result = bsp(paths, name, legend, xmin, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace)
p = result[[1]]
pdf(file = "/Figures/Figures_art_Y_mt/BEAST_YM.pdf", width = 18, height = 9)
p
dev.off()

maxMY = result[[2]]
reduction_factor = result[[3]]

########## VARIANCE ##########
paths = c("/BEAST/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_c/", 
          "/BEAST/unilineal/regular/r=0.01/k=0.05/FT=150/patrilineal_villages_g/", 
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g/", 
          "/BEAST/unilineal/regular/r=0.01/k=0.2/FT=150/patrilineal_villages_g/")
name = c("0", "0.05", "0.1", "0.2")
legend = c("0", "0.05", "0.1", "0.2")
xmin=-100
yminY = 0
ymaxY = 2500
yminM = 0
ymaxM = 2600
yminMY = 0
ymaxMY = 20
segmentY = 2300
textY = 2400
segmentM = 2400
textM = 2500
segmentMY = 18
textMY = 19
col = c("#F0A4C6", "#D976A3", "#B04576", "#8C2353")
sh=c(15, 16, 17, 18)
lt = rep("solid", 4)
fontsize1 = 22
fontsize2 = 18
nCol = 4
legendSpace = 0.2

result = bsp(paths, name, legend, xmin, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace)
p = result[[1]]
pdf(file = "/Figures/Figures_art_Y_mt/BEAST_YM_var.pdf", width = 18, height = 6.5)
p
dev.off()

maxMY = result[[2]]
reduction_factor = result[[3]]

########## FT ##########
paths = c("/BEAST/unilineal/regular/r=0.01/k=0.1/FT=100/patrilineal_villages_g/", 
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g/", 
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=200/patrilineal_villages_g/")
name = c("100", "150", "200")
legend = c("100", "150", "200")
xmin=-100
yminY = 0
ymaxY = 2500
yminM = 0
ymaxM = 2600
yminMY = 0
ymaxMY = 20
segmentY = 2300
textY = 2400
segmentM = 2400
textM = 2500
segmentMY = 18
textMY = 19
col = c("#F0A4C6", "#B04576", "#8C2353")
sh = c(16, 15, 17)
lt = rep("solid", 3)
fontsize1 = 22
fontsize2 = 18
nCol = 3
legendSpace = 0.2

result = bsp(paths, name, legend, xmin, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace)
p = result[[1]]
pdf(file = "/Figures/Figures_art_Y_mt/BEAST_YM_FT.pdf", width = 18, height = 6.5)
p
dev.off()

maxMY = result[[2]]
reduction_factor = result[[3]]

########## POLYGYNY ##########
paths = c("/BEAST/bilateral/regular/r=0.01/patrilocal_kipsigis/")
name = c("Kipsigis-like")
legend = c("Kipsigis-like")
xmin=-100
yminY = 0
ymaxY = 2500
yminM = 0
ymaxM = 2600
yminMY = 0
ymaxMY = 20
segmentY = 2300
textY = 2400
segmentM = 2400
textM = 2500
segmentMY = 18
textMY = 19
col = c("navyblue")
sh = c(8)
lt = c("dashed")
fontsize1 = 22
fontsize2 = 18
nCol = 1
legendSpace = 0.2

result = bsp(paths, name, legend, xmin, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace)
p = result[[1]]
pdf(file = "/Figures/Figures_art_Y_mt/BEAST_YM_polygyny.pdf", width = 18, height = 6.5)
p
dev.off()

maxMY = result[[2]]
reduction_factor = result[[3]]

########## EXTENDED ##########
paths = c("/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_extended/",
          "/BEAST/unilineal/extended/r=0.01/k=0.1/FT=150/patrilineal_2_patrilocal/",
          "/BEAST/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_strict2relaxed/")
name = c("4a", "4b", "4c")
legend = c("4a: patrilineal descent, lineal fission, variance, no violence", 
           "4b: patrilineal descent, lineal fission, variance, no violence \n then bilateral descent", 
           "4c: patrilineal descent, lineal fission, variance, no violence \n then patrilineal descent, lineal fission, no variance, no violence")
xmin1=-100
xmin2=-200
yminY = 0
ymaxY = 5800
yminM = 0
ymaxM = 5800
yminMY = 0
ymaxMY = 20
segmentY = 5100
textY = 5600
segmentM = 5100
textM = 5600
segmentMY = 18
textMY = 19
col = c("#EDCFF9", "#a15cbe", "#742994")
sh = c(15, 16, 17)
lt = c("solid", "solid", "solid")
fontsize1 = 22
fontsize2 = 18
nCol = 1
legendSpace = 0.4

result = bsp_2_transitions(paths, name, legend, xmin1, xmin2, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace)
p = result[[1]]
pdf(file = "/Figures/Figures_art_Y_mt/BEAST_YM_extended.pdf", width = 18, height = 8)
p
dev.off()

maxMY = result[[2]]
reduction_factor = result[[3]]
