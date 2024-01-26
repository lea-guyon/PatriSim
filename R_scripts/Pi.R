library(cowplot)
library(ggplot2)
library(scales)
library(boot)

mean_ic <- function(tab, metric) {
  mean <- c()
  ICsup = c()
  ICinf = c()
  generations = unique(tab$Gen)
  u = qnorm(0.975)
  
  for (gen in generations) {
    temp_tab = tab[tab$Gen==gen,]
    m = mean(temp_tab[,metric], na.rm = TRUE)
    mean = c(mean, m)
    var = sd(temp_tab[,metric], na.rm = TRUE)
    n = length(temp_tab[,metric])
    ICsup = c(ICsup, m + u*var/sqrt(n))
    ICinf = c(ICinf, m - u*var/sqrt(n))
  }
  new_tab <- data.frame("Gen" = unique(tab$Gen), "mean" = mean, "var" = var, "ICinf" = ICinf, "ICsup" = ICsup)
  return(new_tab)
}

mean_by_gen_XA <- function(data, indices, gen) {
  d <- data[indices,]
  d <- d[d$Gen == gen,]
  mean_ax = mean(d$Pi_X)/mean(d$Pi_A)
  return(mean_ax)
}

mean_by_gen_MY <- function(data, indices, gen) {
  d <- data[indices,]
  d <- d[d$Gen == gen,]
  mean_ax = mean(d$Pi_M)/mean(d$Pi_Y)
  return(mean_ax)
}

fun_mean <- function(data, indices, gen) {
  d <- data[indices,]
  d <- d[d$Gen == gen,]
  mean = mean(d$Pi)
  return(mean)
}

prepare_df <- function(model) {
  "
  Description : transforms result tables into data frames with mean pi + 
  compute mean pi X/A with bootstrap for confidence intervals
  descent : character chain indicating the descent rule
  "
  tab_Y <- read.table("Pi_Y_mean_by_rep.txt", header=TRUE)
  generations = unique(tab_Y$Gen)
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab_Y, statistic = fun_mean, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  pi_tab_Y <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_Y$model <- model

  tab_M <- read.table("Pi_Mito_mean_by_rep.txt", header=TRUE)
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab_M, statistic = fun_mean, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  pi_tab_Mito <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_Mito$model <- model
  
  tab <- data.frame("Gen"=tab_M$Gen, "Pi_M"=tab_M$Pi, "Pi_Y" = tab_Y$Pi)
  tab <- na.omit(tab)
  generations = unique(tab_M$Gen)
  
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab, statistic = mean_by_gen_MY, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  
  pi_tab_MY <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_MY$model <- model
  
  return(list(pi_tab_Y, pi_tab_Mito, pi_tab_MY))
}

prepare_df_global <- function(model) {
  "
  Description : separate the table of results into different data frames with mean global pi for each 
  chromosome and for each ratio with bootstrap for confidence intervals
  tab : metric table generated with python script
  model : character chain indicating the model
  "
  tab <- read.table('Global_Pi.txt', header = T)
  
  tab_Y <- tab[tab$CHR=='Y', ]
  generations = unique(tab_Y$Gen)
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab_Y, statistic = fun_mean, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  pi_tab_Y <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_Y$model <- model
  
  tab_M <- tab[tab$CHR=='Mito', ]
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab_M, statistic = fun_mean, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  pi_tab_Mito <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_Mito$model <- model
  
  tab <- data.frame("Gen"=tab_M$Gen, "Pi_M"=tab_M$Pi, "Pi_Y" = tab_Y$Pi)
  tab <- na.omit(tab)
  
  mean <- c()
  ICinf <- c()
  ICsup <- c()
  for (gen in generations) {
    bootstrap <- boot(data=tab, statistic = mean_by_gen_MY, R = 10000, gen=gen)
    ci <- boot.ci(bootstrap, conf = 0.95, type="bca")
    mean <- c(mean, bootstrap$t0)
    ICinf <- c(ICinf, ci$bca[4])
    ICsup <- c(ICsup, ci$bca[5])
  }
  
  pi_tab_MY <- data.frame("Gen" = generations, "mean" = mean, "ICinf" = ICinf, "ICsup" = ICsup)
  pi_tab_MY$model <- model
  
  return(list(pi_tab_Y, pi_tab_Mito, pi_tab_MY))
}

pairwiseWilcoxTest <- function(dir, values) {
  df_Y = NULL
  df_M = NULL
  df_MY = NULL
  for (i in values) {
    setwd(paste(dir, i, sep = ""))
    tab_Y <- read.table("Pi_Y_mean_by_rep.txt", header=TRUE)
    tab_Y$model = rep(i, length(tab_Y$Pi))
    
    tab_M <- read.table("Pi_Y_mean_by_rep.txt", header=TRUE)
    tab_M$model = rep(i, length(tab_M$Pi))
    
    tab_MY <- read.table("Pi_Y_mean_by_rep.txt", header=TRUE)
    tab_MY$model = rep(i, length(tab_MY$Pi))
    
    if (is.null(df_Y)) {
      df_Y = data.frame("Gen" = tab_Y$Gen, "Pi" = tab_Y$Pi, "model" = tab_Y$model)
      df_M = data.frame("Gen" = tab_M$Gen, "Pi" = tab_M$Pi, "model" = tab_M$model)
      df_MY = data.frame("Gen" = tab_MY$Gen, "Pi" = tab_MY$Pi, "model" = tab_MY$model)
    }
    else {
      df_Y = rbind.fill(df_Y, tab_Y)
      df_M = rbind.fill(df_M, tab_M)
      df_MY = rbind.fill(df_MY, tab_MY)
    }
  }
  df_Y_100 = df_Y[df_Y$Gen == 100,]
  df_M_100 = df_M[df_M$Gen == 100,]
  Wtest = pairwise.wilcox.test(df_M_100$Pi, df_M_100$model, p.adjust.method = "bonferroni")
  return(Wtest)
}

pi_tab <- function(paths, name, scale='local') {
  posd <- position_dodge(5)
  
  flag = F
  for (i in seq(1, length(paths))) {
    path = paths[i]
    model = name[i]
    
    if (flag == F) {
      var = c('pi_tab_Y', 'pi_tab_Mito', 'pi_tab_MY')
    }
    else {
      var = c('pi_tab_bis_Y', 'pi_tab_bis_Mito', 'pi_tab_bis_MY')
    }
    
    setwd(path)
    if (scale == 'global') {
      val = prepare_df_global(model)
    }
    else {
      val = prepare_df(model)
    }
    for (i in seq(var)) assign(var[i], val[i])
    
    if (flag == T) {
      # merge tables
      pi_tab_Y <- merge(pi_tab_Y, pi_tab_bis_Y, all = TRUE)
      pi_tab_Mito <- merge(pi_tab_Mito, pi_tab_bis_Mito, all = TRUE)
      pi_tab_MY <- merge(pi_tab_MY, pi_tab_bis_MY, all = TRUE)
    }
    
    flag = T
  }
  return(list(pi_tab_Y, pi_tab_Mito, pi_tab_MY))
}

figure <- function(pi_tab_Y, pi_tab_Mito, pi_tab_MY, legend, xmin, yminY, ymaxY,  yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, legendSpace) {
  posd <- position_dodge(5)
  # Build legend
  ggp_split_1 <- ggplot(pi_tab_Y,                    # Create ggplot2 plot of first subset
                        aes(Gen,
                            mean,
                            color = model,
                            shape = model,
                            linetype = model)) +
    geom_errorbar(aes(ymin = ICinf,  ymax = ICsup), position = posd) +
    geom_point(size = 5, stroke = 2) +
    scale_color_manual(labels = legend,
                       values = col) +
    scale_shape_manual(labels = legend,
                       values = sh) +
    scale_linetype_manual(labels = legend,
                          values = lt) +
    labs(color = "Scenarios", shape = "Scenarios", linetype = "Scenarios") +
    theme(legend.key.size = unit(2, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.title = element_text(face = 'bold', size = 22),
          legend.text = element_text(size=19),
          legend.spacing.x = unit(1.0, 'cm')) + 
    guides(color = guide_legend(ncol=nCol, bycol=TRUE))
  ggp_legend_split_1 <- get_legend(ggp_split_1)
  
  pi_tab_Y$Gen = pi_tab_Y$Gen + xmin
  
  p1 <- ggplot(pi_tab_Y, aes(Gen, mean/(2*2.5e-8), shape = model, col= model)) +
    theme_cowplot() +
    theme(text=element_text(size=fontsize1),
          axis.text = element_text(size=fontsize2),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(color="bisque")) +
    annotate("rect", xmin = xmin, xmax = 2, ymin = yminY, ymax = ymaxY,
             alpha = .1,fill = "orange") +
    labs(x= 'Generation', y = expression(''*N[e]^(Y)*'')) +
    scale_x_continuous(breaks = seq(xmin, 0, 20)) +
    scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(yminY, ymaxY)) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    geom_line(aes(linetype = model), position = posd, alpha = 0.5) +
    geom_pointrange(aes(ymin=ICinf/(2*2.5e-8), ymax=ICsup/(2*2.5e-8)), position = posd, size = 1.5, stroke = 2) +
    annotate("segment", x = xmin, xend = xmin, y = yminY, yend = segmentY, colour = "salmon2", linetype='longdash', linewidth = 1) +
    annotate("text", x = xmin, y = textY, label = bquote(''*t[0]*''), col = 'salmon2', size = 8)
  
  pi_tab_Mito$Gen = pi_tab_Mito$Gen + xmin
  p2 <- ggplot(pi_tab_Mito, aes(Gen, mean/(2*5.5e-7), shape = model, col= model)) +
    theme_cowplot() +
    theme(text=element_text(size=fontsize1),
          axis.text = element_text(size=fontsize2),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(color="bisque")) +
    annotate("rect", xmin = xmin, xmax = 2, ymin = yminM, ymax = ymaxM,
             alpha = .1,fill = "orange") +
    labs(x= 'Generation', y = expression(''*N[e]^(mt)*'')) +
    scale_x_continuous(breaks = seq(xmin,0,20)) +
    scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(yminM, ymaxM)) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    geom_line(aes(linetype = model), position = posd, alpha = 0.5) +
    geom_pointrange(aes(ymin=ICinf/(2*5.5e-7), ymax=ICsup/(2*5.5e-7)), position = posd, size = 1.5, stroke = 2) +
    annotate("segment", x = xmin, xend = xmin, y = yminM, yend = segmentM, colour = "salmon2", linetype='longdash', linewidth = 1) +
    annotate("text", x = xmin, y = textM, label = bquote(''*t[0]*''), col = 'salmon2', size = 8)
  
  pi_tab_MY$Gen = pi_tab_MY$Gen + xmin
  p3 <- ggplot(pi_tab_MY, aes(Gen, mean/22, col = model, shape = model)) +
    theme_cowplot() +
    theme(text=element_text(size=fontsize1),
          axis.text = element_text(size=fontsize2),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(color="bisque")) +
    annotate("rect", xmin = xmin, xmax = 2, ymin = yminMY, ymax = ymaxMY,
             alpha = .1,fill = "orange") +
    labs(x= 'Generation', y = expression(''*N[e]^(mt/Y)*'')) +
    scale_x_continuous(breaks = seq(xmin, 0, 20)) +
    scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(yminMY, ymaxMY)) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    geom_line(aes(Gen, mean/22, col = model, linetype = model), position = posd, alpha = 0.5) +
    geom_abline(slope = 0, intercept=1, color='orchid4', linetype='longdash', linewidth = 0.75) +
    geom_pointrange(aes(ymin=ICinf/22, ymax=ICsup/22), position = posd, size = 1.5, stroke = 2) +
    annotate("segment", x = xmin, xend = xmin, y = yminMY, yend = segmentMY, colour = "salmon2", linetype='longdash', linewidth = 1) +
    annotate("text", x = xmin, y = textMY, label = bquote(''*t[0]*''), col = 'salmon2', size = 8)
  
  
  p <- plot_grid(p1 + theme(legend.position="none"), 
                 p2 + theme(legend.position="none"),
                 p3 + theme(legend.position="none"),
                 ncol=3, labels = "auto", label_size = fontsize1)
  
  p_bis <- plot_grid(p, legend=ggp_legend_split_1, ncol = 1, rel_heights = c(1-legendSpace, legendSpace))
  return(p_bis)
}

#setwd("/set/your/working/directory")

paths = c("/Tables/Pi/bilateral/regular/r=0.01/patrilocal_villages", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_a",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/e=0.15/patrilineal_villages_b",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_c",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/e=0.15/patrilineal_villages_d",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_e", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/e=0.15/patrilineal_villages_f",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g",  
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/e=0.15/patrilineal_villages_h")
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
xmin=-100
yminY = 0
ymaxY = 990
yminM = 650
ymaxM = 900
yminMY = 0
ymaxMY = 16
segmentY = 900
textY = 940
segmentM = 880
textM = 890
segmentMY = 14
textMY = 15
col = c("#7bc84e", "#69aa63", "#55c6b2", "#5a81bd", "#4c69d5", 
        "#bd634a", "#d5732d", "#a15cbe", "#c15787")
sh = c(8, 22, 23, 21, 24, 15, 18, 16, 17)
lt = c(rep("dashed", 5), rep("solid", 4))
fontsize1 = 22
fontsize2 = 18
nCol=2
legendSpace=0.4

pi_tabs = pi_tab(paths, name)
pi_tab_Y = pi_tabs[[1]]
pi_tab_M = pi_tabs[[2]]
pi_tab_MY = pi_tabs[[3]]
result = figure(pi_tab_Y, pi_tab_M, pi_tab_MY, legend, xmin, yminY, ymaxY, 
                yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, 
                segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, 
                legendSpace)

pdf(file = "/Figures/Ne_Y_Mito.pdf", width = 18, height = 9)
result
dev.off()

# Y
csv_Y <- pi_tab_Y
csv_Y$mean <- pi_tab_Y$mean/(2*2.5e-8)
csv_Y$ICinf <- pi_tab_Y$ICinf/(2*2.5e-8)
csv_Y$ICsup <- pi_tab_Y$ICsup/(2*2.5e-8)
csv_Y$Gen <- pi_tab_Y$Gen - 100
csv_Y = csv_Y[, c(5, 1, 2, 3, 4)]
csv_Y = csv_Y[order(csv_Y$Gen, csv_Y$model),]
write.csv(csv_Y, "/Tables/Pi/Ne_Y.csv", row.names = F)

# Mito
csv_Mito = pi_tab_Mito
csv_Mito$mean <- pi_tab_Mito$mean/(2*5.5e-7)
csv_Mito$ICinf <- pi_tab_Mito$ICinf/(2*5.5e-7)
csv_Mito$ICsup <- pi_tab_Mito$ICsup/(2*5.5e-7)
csv_Mito$Gen <- pi_tab_Mito$Gen - 100
csv_Mito = csv_Mito[, c(5, 1, 2, 3, 4)]
csv_Mito = csv_Mito[order(csv_Mito$Gen, csv_Mito$model),]
write.csv(csv_Mito, "/Tables/Pi/Ne_Mito.csv", row.names = F)

# M/Y
csv_MY = pi_tab_MY
csv_MY$mean <- pi_tab_MY$mean/(22)
csv_MY$ICinf <- pi_tab_MY$ICinf/(22)
csv_MY$ICsup <- pi_tab_MY$ICsup/(22)
csv_MY$Gen <- pi_tab_MY$Gen - 100
csv_MY = csv_MY[, c(5, 1, 2, 3, 4)]
csv_MY = csv_MY[order(csv_MY$Gen, csv_MY$model),]
write.csv(csv_MY, "/Tables/Pi/Ne_MY.csv", row.names = F)

######## VARIANCE #########

paths = c("/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_c", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.05/FT=150/patrilineal_villages_g", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.2/FT=150/patrilineal_villages_g")
name = c("0", "0.05", "0.1", "0.2")
legend = c("0", "0.05", "0.1", "0.2")
xmin=-100
yminY = 0
ymaxY = 980
yminM = 660
ymaxM = 870
yminMY = 0
ymaxMY = 34
segmentY = 880
textY = 930
segmentM = 850
textM = 860
segmentMY = 30
textMY = 32
col = c("#F0A4C6", "#D976A3", "#B04576", "#8C2353")
sh = c(16, 15, 17, 18)
lt = rep("solid", 4)
fontsize1 = 22
fontsize2 = 18
nCol=4
legendSpace=0.2

pi_tabs = pi_tab(paths, name)
pi_tab_Y = pi_tabs[[1]]
pi_tab_M = pi_tabs[[2]]
pi_tab_MY = pi_tabs[[3]]
result = figure(pi_tab_Y, pi_tab_M, pi_tab_MY, legend, xmin, yminY, ymaxY, 
                yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, 
                segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol,
                legendSpace)

pdf(file = "/Figures/Ne_Y_Mito_var.pdf", width = 18, height = 7)
result
dev.off()

# Y
csv_Y <- pi_tab_Y
csv_Y$mean <- pi_tab_Y$mean/(2*2.5e-8)
csv_Y$ICinf <- pi_tab_Y$ICinf/(2*2.5e-8)
csv_Y$ICsup <- pi_tab_Y$ICsup/(2*2.5e-8)
csv_Y$Gen <- pi_tab_Y$Gen - 100
csv_Y = csv_Y[, c(5, 1, 2, 3, 4)]
csv_Y = csv_Y[order(csv_Y$Gen, csv_Y$model),]
write.csv(csv_Y, "/Tables/Pi/Ne_Y_var.csv", row.names = F)

# Mito
csv_Mito = pi_tab_M
csv_Mito$mean <- pi_tab_M$mean/(2*5.5e-7)
csv_Mito$ICinf <- pi_tab_M$ICinf/(2*5.5e-7)
csv_Mito$ICsup <- pi_tab_M$ICsup/(2*5.5e-7)
csv_Mito$Gen <- pi_tab_M$Gen - 100
csv_Mito = csv_Mito[, c(5, 1, 2, 3, 4)]
csv_Mito = csv_Mito[order(csv_Mito$Gen, csv_Mito$model),]
write.csv(csv_Mito, "/Tables/Pi/Ne_Mito_var.csv", row.names = F)

# M/Y
csv_MY = pi_tab_MY
csv_MY$mean <- pi_tab_MY$mean/(22)
csv_MY$ICinf <- pi_tab_MY$ICinf/(22)
csv_MY$ICsup <- pi_tab_MY$ICsup/(22)
csv_MY$Gen <- pi_tab_MY$Gen - 100
csv_MY = csv_MY[, c(5, 1, 2, 3, 4)]
csv_MY = csv_MY[order(csv_MY$Gen, csv_MY$model),]
write.csv(csv_MY, "/Tables/Pi/Ne_MY_var.csv", row.names = F)

########### FT ############

paths = c("/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=100/patrilineal_villages_g", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=200/patrilineal_villages_g")
name = c("100", "150", "200")
legend = c("100", "150", "200")
xmin=-100
yminY = 0
ymaxY = 980
yminM = 660
ymaxM = 870
yminMY = 0
ymaxMY = 27
segmentY = 880
textY = 930
segmentM = 850
textM = 860
segmentMY = 23
textMY = 25
col = c("#F0A4C6", "#B04576", "#8C2353")
sh = c(16, 15, 17)
lt = rep("solid", 3)
fontsize1 = 22
fontsize2 = 18
nCol=3
legendSpace=0.2

pi_tabs = pi_tab(paths, name)
pi_tab_Y = pi_tabs[[1]]
pi_tab_M = pi_tabs[[2]]
pi_tab_MY = pi_tabs[[3]]
result = figure(pi_tab_Y, pi_tab_M, pi_tab_MY, legend, xmin, yminY, ymaxY, 
                yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, 
                segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol,
                legendSpace)

pdf(file = "/Figures/Ne_Y_Mito_FT.pdf", width = 18, height = 7)
result
dev.off()

# Y
csv_Y <- pi_tab_Y
csv_Y$mean <- pi_tab_Y$mean/(2*2.5e-8)
csv_Y$ICinf <- pi_tab_Y$ICinf/(2*2.5e-8)
csv_Y$ICsup <- pi_tab_Y$ICsup/(2*2.5e-8)
csv_Y$Gen <- pi_tab_Y$Gen - 100
csv_Y = csv_Y[, c(5, 1, 2, 3, 4)]
csv_Y = csv_Y[order(csv_Y$Gen, csv_Y$model),]
write.csv(csv_Y, "/Tables/Pi/Ne_Y_FT.csv", row.names = F)

# Mito
csv_Mito = pi_tab_M
csv_Mito$mean <- pi_tab_M$mean/(2*5.5e-7)
csv_Mito$ICinf <- pi_tab_M$ICinf/(2*5.5e-7)
csv_Mito$ICsup <- pi_tab_M$ICsup/(2*5.5e-7)
csv_Mito$Gen <- pi_tab_M$Gen - 100
csv_Mito = csv_Mito[, c(5, 1, 2, 3, 4)]
csv_Mito = csv_Mito[order(csv_Mito$Gen, csv_Mito$model),]
write.csv(csv_Mito, "/Tables/Pi/Ne_Mito_FT.csv", row.names = F)

# M/Y
csv_MY = pi_tab_MY
csv_MY$mean <- pi_tab_MY$mean/(22)
csv_MY$ICinf <- pi_tab_MY$ICinf/(22)
csv_MY$ICsup <- pi_tab_MY$ICsup/(22)
csv_MY$Gen <- pi_tab_MY$Gen - 100
csv_MY = csv_MY[, c(5, 1, 2, 3, 4)]
csv_MY = csv_MY[order(csv_MY$Gen, csv_MY$model),]
write.csv(csv_MY, "/Tables/Pi/Ne_MY_FT.csv", row.names = F)

######## POLYGYNY #########

paths = c("/Tables/Pi/bilateral/regular/r=0.01/patrilocal_kipsigis")
name = c("Kipsigis-like")
legend = c("Kipsigis-like")
xmin=-100
yminY = 0
ymaxY = 980
yminM = 640
ymaxM = 870
yminMY = 0
ymaxMY = 5.8
segmentY = 880
textY = 930
segmentM = 850
textM = 860
segmentMY = 5
textMY = 5.4
col = c("navyblue")
sh = c(8)
lt = rep("dashed", 1)
fontsize1 = 22
fontsize2 = 18
nCol=3
legendSpace=0.2

pi_tabs = pi_tab(paths, name)
pi_tab_Y = pi_tabs[[1]][[1]]
pi_tab_M = pi_tabs[[2]][[1]]
pi_tab_MY = pi_tabs[[3]][[1]]
result = figure(pi_tab_Y, pi_tab_M, pi_tab_MY, legend, xmin, yminY, ymaxY,
                yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, 
                segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol,
                legendSpace)

pdf(file = "/Figures/Ne_Y_Mito_polygyny.pdf", width = 18, height = 7)
result
dev.off()

# Y
csv_Y <- pi_tab_Y
csv_Y$mean <- pi_tab_Y$mean/(2*2.5e-8)
csv_Y$ICinf <- pi_tab_Y$ICinf/(2*2.5e-8)
csv_Y$ICsup <- pi_tab_Y$ICsup/(2*2.5e-8)
csv_Y$Gen <- pi_tab_Y$Gen - 100
csv_Y = csv_Y[, c(5, 1, 2, 3, 4)]
csv_Y = csv_Y[order(csv_Y$Gen, csv_Y$model),]
write.csv(csv_Y, "/Tables/Pi/Ne_Y_poly.csv", row.names = F)

# Mito
csv_Mito = pi_tab_M
csv_Mito$mean <- pi_tab_M$mean/(2*5.5e-7)
csv_Mito$ICinf <- pi_tab_M$ICinf/(2*5.5e-7)
csv_Mito$ICsup <- pi_tab_M$ICsup/(2*5.5e-7)
csv_Mito$Gen <- pi_tab_M$Gen - 100
csv_Mito = csv_Mito[, c(5, 1, 2, 3, 4)]
csv_Mito = csv_Mito[order(csv_Mito$Gen, csv_Mito$model),]
write.csv(csv_Mito, "/Tables/Pi/Ne_Mito_poly.csv", row.names = F)

# M/Y
csv_MY = pi_tab_MY
csv_MY$mean <- pi_tab_MY$mean/(22)
csv_MY$ICinf <- pi_tab_MY$ICinf/(22)
csv_MY$ICsup <- pi_tab_MY$ICsup/(22)
csv_MY$Gen <- pi_tab_MY$Gen - 100
csv_MY = csv_MY[, c(5, 1, 2, 3, 4)]
csv_MY = csv_MY[order(csv_MY$Gen, csv_MY$model),]
write.csv(csv_MY, "/Tables/Pi/Ne_MY_poly.csv", row.names = F)

########## GLOBAL ###########
paths = c("/Tables/Pi/bilateral/regular/r=0.01/patrilocal_villages", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_a",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/e=0.15/patrilineal_villages_b",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/patrilineal_villages_c",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0/FT=150/e=0.15/patrilineal_villages_d",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_e", 
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/e=0.15/patrilineal_villages_f",
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/patrilineal_villages_g",  
          "/Tables/Pi/unilineal/regular/r=0.01/k=0.1/FT=150/e=0.15/patrilineal_villages_h")
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
xmin=-100
yminY = 0
ymaxY = 990
yminM = 650
ymaxM = 890
yminMY = 0
ymaxMY=2
segmentY = 900
textY = 940
segmentM = 870
textM = 880
segmentMY = 1.8
textMY = 1.9
col = c("#7bc84e", "#69aa63", "#55c6b2", "#5a81bd", "#4c69d5", 
        "#bd634a", "#d5732d", "#a15cbe", "#c15787")
sh = c(8, 22, 23, 21, 24, 15, 18, 16, 17)
lt = c(rep("dashed", 5), rep("solid", 4))
fontsize1 = 22
fontsize2 = 18
nCol=2
legendSpace=0.4

pi_tabs = pi_tab(paths, name, "global")
pi_tab_Y = pi_tabs[[1]]
pi_tab_M = pi_tabs[[2]]
pi_tab_MY = pi_tabs[[3]]
result = figure(pi_tab_Y, pi_tab_M, pi_tab_MY, legend, xmin, yminY, ymaxY, 
                yminM, ymaxM, yminMY, ymaxMY, segmentY, textY, segmentM, textM, 
                segmentMY, textMY, col, sh, lt, fontsize1, fontsize2, nCol, 
                legendSpace)

pdf(file = "/Figures/Ne_Y_Mito_global.pdf", width = 18, height = 9)
result
dev.off()

# Y
csv_Y <- pi_tab_Y
csv_Y$mean <- pi_tab_Y$mean/(2*2.5e-8)
csv_Y$ICinf <- pi_tab_Y$ICinf/(2*2.5e-8)
csv_Y$ICsup <- pi_tab_Y$ICsup/(2*2.5e-8)
csv_Y$Gen <- pi_tab_Y$Gen - 100
csv_Y = csv_Y[, c(5, 1, 2, 3, 4)]
csv_Y = csv_Y[order(csv_Y$Gen, csv_Y$model),]
write.csv(csv_Y, "/Tables/Pi/Ne_Y_global.csv", row.names = F)

# Mito
csv_Mito = pi_tab_M
csv_Mito$mean <- pi_tab_M$mean/(2*5.5e-7)
csv_Mito$ICinf <- pi_tab_M$ICinf/(2*5.5e-7)
csv_Mito$ICsup <- pi_tab_M$ICsup/(2*5.5e-7)
csv_Mito$Gen <- pi_tab_M$Gen - 100
csv_Mito = csv_Mito[, c(5, 1, 2, 3, 4)]
csv_Mito = csv_Mito[order(csv_Mito$Gen, csv_Mito$model),]
write.csv(csv_Mito, "/Tables/Pi/Ne_Mito_global.csv", row.names = F)

# M/Y
csv_MY = pi_tab_MY
csv_MY$mean <- pi_tab_MY$mean/(22)
csv_MY$ICinf <- pi_tab_MY$ICinf/(22)
csv_MY$ICsup <- pi_tab_MY$ICsup/(22)
csv_MY$Gen <- pi_tab_MY$Gen - 100
csv_MY = csv_MY[, c(5, 1, 2, 3, 4)]
csv_MY = csv_MY[order(csv_MY$Gen, csv_MY$model),]
write.csv(csv_MY, "/Tables/Pi/Ne_MY_global.csv", row.names = F)
