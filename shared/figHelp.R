####Plot CIs as polygons

polyCI <- function(upper, lower, col){
  len <- length(upper)
  polygon(x = c(rev(1 : len), 1 : len)
          , y = c(rev(lower), upper)
          , col = col
          , border = NA)
}


#####Spike-and-Slab estimate figure#########
##1. posterior z est vs observed effect
##2. Model estimates with CI

mix.fig.check <- function(effect, z, theta, main = c("A.", "B."), int = c(.1, .1)){
  
  #probability plot
  plot(effect, colMeans(z)
       , pch = 19
       , col = "gray40"
       , ylab = "Slab probability"
       , xlab = "Observed effect"
       , axes = F
       , main = main[1])
  abline(h = .5, col = "gray70")
  low <- ceiling(min(effect * 100))/100
  hig <- floor(max(effect * 100))/100
  axis(side = 1, seq(low, hig, int[1]))
  axis(side = 2, seq(round(min(colMeans(z)), 1), 1, int[2]))
  
  #Effect and estimate plot
  ind <- order(effect)
  est <- z * theta
  upper <- apply(est, 2, quantile, probs = .975)
  lower <- apply(est, 2, quantile, probs = .025)
  thneg <- effect[effect < 0]
  
  plot(1:ncol(z), effect[ind]
       , pch = 19
       , ylab = expression(hat(theta))
       , xlab = "Participants"
       , col = "darkgray"
       , axes = F
       , main = main[2])
  axis(side = 1, c(1, ncol(z)))
  axis(side = 2, seq(low, hig, int[1]))
  
  polyCI(upper[ind], lower[ind], "slategray1")
  points(1:ncol(z)
         , effect[ind]
         , pch = 19
         , col = "darkgray")
  points(1:length(thneg)
         , sort(thneg)
         , col = "darkred")
  points(1:ncol(z)
         , colMeans(est)[ind]
         , pch = 19
         , col = "skyblue4")
  abline(h = 0)
}


#####Spike-and-Slab estimate figure#########
##1. posterior z est vs observed effect
##2. Model estimates with CI

mix.fig <- function(effect, z, theta, main = c("A.", "B."), int = c(.1, .1), bayesfactors, win, new = F){
  
  cols <- sequential_hcl(100, h = 260, c. = c(80, 80), l = c(17, 90))
  
  #Effect and estimate plot
  ind <- order(effect)
  est <- z * theta
  upper <- apply(est, 2, quantile, probs = .975)
  lower <- apply(est, 2, quantile, probs = .025)
  thneg <- effect[effect < 0]
  probs <- round(colMeans(z)*100)
  
  low <- ceiling(min(effect * 100))/100
  hig <- floor(max(upper * 100))/100
  
  plot(1:ncol(z), effect[ind]
       , pch = 3
       , ylab = expression(hat(theta))
       , xlab = "Participants"
       , col = "darkgray"
       , axes = F
       # , main = main[1]
       , ylim = c(min(effect, lower), max(effect, upper))
       , cex = 1.2
       , cex.lab = 1.2)
  axis(side = 1, c(1, ncol(z)))
  axis(side = 2, seq(low, hig, int[1]))
  title(main[1], line = -1)
  
  polyCI(upper[ind], lower[ind], "gainsboro")
  points(1:ncol(z)
         , effect[ind]
         , pch = 3
         , col = "gray40"
         , cex = 1.2)
  if(new == F){
  points(1:length(thneg)
         , sort(thneg)
         , col = "darkred"
         , pch = 3
         , cex = 1.2)}
  points(1:ncol(z)
         , colMeans(est)[ind]
         , pch = 19
         , col = cols[probs[ind]]
         , cex = 1.3)
  points(1:ncol(z)
         , colMeans(est)[ind]
         , col = "gray40"
         , cex = 1.3)
  abline(h = 0)
  
  ##BF
  safe <- par('mgp', 'mar', 'xpd')
  par(mgp = c(2, 1, 0), mar = c(1, 1, 0, 1), xpd = T)
  
  names <- c(
    paste("Spike-and-slab")
    , paste("Unconstrained")
    , paste("Positive-effects")
    , paste("Common-effect")
    , paste("Null")
  )
  
  myCol <- c("winner" = "firebrick"
             , "regular" = "darkgray"
             , "encompassing" = "slategray4"
             , "analytical" = "slategray4")
  
  pos <- matrix(ncol = 2, nrow = length(names))
  x <- c(1/6, .5, 5/6, .25)
  pos[,1] <- c(x[3], x[2], x[1], x[4], x[2])
  pos[,2] <- rev(seq(.1, .9, .2))
  pos[3, 2] <- .55
  pos[4, 2] <- .25
  
  bftext <- bayesfactors
  
  connect <- matrix(0, nrow = length(names), byrow = F, ncol = length(names))
  connect[1, 2] <- ''
  connect[1, 5] <- ''
  connect[2, 3] <- ''
  connect[2, 4] <- ''
  connect[2, 5] <- ''
  
  lcol <- matrix(0, nrow = length(names), byrow = F, ncol = length(names))
  lcol[1, 2] <- myCol["encompassing"]
  lcol[1, 5] <- myCol["encompassing"]
  lcol[2, 3] <- myCol["analytical"]
  lcol[2, 4] <- myCol["analytical"]
  lcol[2, 5] <- myCol["analytical"]
  
  col.order <- rep("regular", 5)
  col.order[win] <- "winner"
  
  G <- plotmat(connect
               , pos
               , name = names
               , curve = 0
               , box.type="rect"
               , box.size=.22
               , box.prop=.15
               , box.lcol = myCol[col.order]
               , arr.length = 0
               , arr.lcol = lcol
               , shadow.size = 0
               , box.cex = 1.1
               , dtext = .1
               # , main = main[2]
               , cex.main = 1.2
  )
  title(main[2], line = -1)
  
  x <- c(.625, .695, .286, .415, .53)
  y <- c(.8, .5, .625, .49, .4)
  srt <- c(36, 70, 28, 65, 89)
  # offset <- c(.3, )
  for(i in 1:length(bftext)){
    text(x = x[i]
         , y = y[i]
         , labels = bftext[i]
         , cex = 1
         , srt = srt[i]
    )
  }
  par(safe)
  
}


mix.fig.simple <- function(effect, z, theta, main = "A.", int = c(.1, .1)){
  
  cols <- sequential_hcl(100, h = 260, c. = c(80, 80), l = c(17, 90))
  
  #Effect and estimate plot
  ind <- order(effect)
  est <- z * theta
  upper <- apply(est, 2, quantile, probs = .975)
  lower <- apply(est, 2, quantile, probs = .025)
  thneg <- effect[effect < 0]
  probs <- round(colMeans(z)*100)
  
  low <- ceiling(min(effect * 100))/100
  hig <- floor(max(upper * 100))/100
  
  plot(1:ncol(z), effect[ind]
       , pch = 3
       , ylab = expression(hat(theta))
       , xlab = "Participants"
       , col = "darkgray"
       , axes = F
       , main = main[1]
       , ylim = c(min(effect, lower), max(effect, upper))
       , cex = 1.2
       , cex.lab = 1.2)
  axis(side = 1, c(1, ncol(z)))
  axis(side = 2, seq(low, hig, int[1]))
  
  polyCI(upper[ind], lower[ind], "gainsboro")
  points(1:ncol(z)
         , effect[ind]
         , pch = 3
         , col = "gray40"
         , cex = 1.2)
  points(1:length(thneg)
         , sort(thneg)
         , col = "darkred"
         , pch = 3
         , cex = 1.2)
  points(1:ncol(z)
         , colMeans(est)[ind]
         , pch = 19
         , col = cols[probs[ind]]
         , cex = 1.3)
  points(1:ncol(z)
         , colMeans(est)[ind]
         , col = "gray40"
         , cex = 1.3)
  abline(h = 0)
}




mix.fig.poster <- function(effect, z, theta, main = c("A.", "B."), int = c(.1, .1), bayesfactors, win, new = F){
  
  cols <- sequential_hcl(100, h = 260, c. = c(80, 80), l = c(17, 90))
  
  #Effect and estimate plot
  ind <- order(effect)
  est <- z * theta
  upper <- apply(est, 2, quantile, probs = .975)
  lower <- apply(est, 2, quantile, probs = .025)
  thneg <- effect[effect < 0]
  probs <- round(colMeans(z)*100)
  
  low <- ceiling(min(effect * 100))/100
  hig <- floor(max(upper * 100))/100
  
  plot(1:ncol(z), effect[ind]
       , pch = 3
       , ylab = expression(hat(theta))
       , xlab = "Participants"
       , col = "darkgray"
       , axes = F
       # , main = main[1]
       , ylim = c(min(effect, lower), max(effect, upper))
       , cex = 1.7
       , cex.lab = 1.5)
  axis(side = 1, c(1, ncol(z)), cex.axis = 1.5)
  axis(side = 2, seq(low, hig, int[1]), cex.axis = 1.5)
  mtext(main[1], 2, line = 14.5, cex = 1.2, adj = 0, las = 2)
  
  polyCI(upper[ind], lower[ind], "gainsboro")
  points(1:ncol(z)
         , effect[ind]
         , pch = 3
         , col = "gray40"
         , cex = 1.2)
  if(new == F){
    points(1:length(thneg)
           , sort(thneg)
           , col = "darkred"
           , pch = 3
           , cex = 1.2)}
  points(1:ncol(z)
         , colMeans(est)[ind]
         , pch = 19
         , col = cols[probs[ind]]
         , cex = 1.3)
  points(1:ncol(z)
         , colMeans(est)[ind]
         , col = "gray40"
         , cex = 1.3)
  abline(h = 0)
  
  ##BF
  safe <- par('mgp', 'mar', 'xpd')
  par(mgp = c(2, 1, 0), mar = c(1, 1, 0, 1), xpd = T)
  
  names <- c(
    paste("Spike-and-slab")
    , paste("Unconstrained")
    , paste("Positive-effects")
    , paste("Common-effect")
    , paste("Null")
  )
  
  myCol <- c("winner" = "firebrick"
             , "regular" = "darkgray"
             , "encompassing" = "slategray4"
             , "analytical" = "slategray4")
  
  pos <- matrix(ncol = 2, nrow = length(names))
  x <- c(.19, .5, .8, .45)
  pos[,1] <- c(x[3], x[2], x[1], x[4], x[3])
  pos[1, 2] <- .95
  pos[2, 2] <- .65
  pos[3, 2] <- .4
  pos[4, 2] <- .2
  pos[5, 2] <- .05
  
  bftext <- bayesfactors
  
  connect <- matrix(0, nrow = length(names), byrow = F, ncol = length(names))
  connect[1, 2] <- ''
  connect[1, 5] <- ''
  connect[2, 3] <- ''
  connect[2, 4] <- ''
  connect[2, 5] <- ''
  
  lcol <- matrix(0, nrow = length(names), byrow = F, ncol = length(names))
  lcol[1, 2] <- myCol["encompassing"]
  lcol[1, 5] <- myCol["encompassing"]
  lcol[2, 3] <- myCol["analytical"]
  lcol[2, 4] <- myCol["analytical"]
  lcol[2, 5] <- myCol["analytical"]
  
  col.order <- rep("regular", 5)
  col.order[win] <- "winner"
  
  G <- plotmat(connect
               , pos
               , name = names
               , curve = 0
               , box.type="rect"
               , box.size=.23
               , box.prop=.15
               , box.lcol = myCol[col.order]
               , arr.length = 0
               , arr.lcol = lcol
               , shadow.size = 0
               , box.cex = 1.5
               , dtext = .1
               # , main = main[2]
               , cex.main = 1.5
  )
  title(main[2], line = -1)
  
  x <- c(.61, .83, .286, .505, .68)
  y <- c(.81, .5, .525, .41, .37)
  srt <- c(45, 89, 39, 84, -63)
  # offset <- c(.3, )
  for(i in 1:length(bftext)){
    text(x = x[i]
         , y = y[i]
         , labels = bftext[i]
         , cex = 1.4
         , srt = srt[i]
    )
  }
  par(safe)
  
}
