prep.models <- function(sub, cond){
  
  I <- length(unique(sub))
  R <- length(sub)
  X.full <- matrix(nrow = R, ncol = 2 * I + 2, 0)
  for (r in 1:R){
    X.full[r, 1] <- 1
    X.full[r, sub[r] + 1] <- 1
    if (cond[r] == 2) {
      X.full[r, I + 2] <- 1
      X.full[r, I + 2 + sub[r]] <- 1}}
  
  gMap.full <- c(rep(0, I), 1, rep(2, I))
  
  X.one <- matrix(nrow = R, ncol = I + 2, 0)
  for (r in 1:R){
    X.one[r, 1] <- 1
    X.one[r, sub[r] + 1] <- 1
    if (cond[r] == 2) {
      X.one[r, I + 2] <- 1
    }}
  
  gMap.one <- c(rep(0, I), 1)
  
  X.null <- matrix(nrow = R, ncol = I + 1, 0)
  for(r in 1:R){
    X.null[r, 1] <- 1
    X.null[r, sub[r] + 1] <- 1
  }
  
  gMap.null <- rep(0, I)
  
  if(length(cond) != length(sub)) return(print("oops, your condition vector does not match your subject vector"))
  return(list(X.full = X.full
              , gMap.full = gMap.full
              , X.one = X.one
              , gMap.one = gMap.one
              , X.null = X.null
              , gMap.null = gMap.null
              , R = R
              , I= I))
}

ratio.greater <- function(M, I, a = alpha, b = beta, sd_mu = sigma_mu){
  
  x <- matrix(nrow = M, ncol = I)
  s2 <- rinvgamma(M, a, b)
  mu <- rcauchy(M, 0, sd_mu)
  for (i in 1:I)
    x[, i] <- rnorm(M, mu, sqrt(s2))
  
  all.greater <- function(x) as.integer(mean(x > 0) == 1)
  return(1/mean(apply(x, 1, all.greater)))
}

makeBF <- function(y, meanScale, effectScale, prep = prep.1, keep = 1001:10000)
{
  mcmc.full <- nWayAOV(y
                       , prep$X.full
                       , prep$gMap.full
                       , rscale = c(1, meanScale, effectScale)
                       , posterior = T
                       , method = "auto"
                       , iterations = max(keep))
  bf.full <- nWayAOV(y
                     , prep$X.full
                     , prep$gMap.full
                     , rscale = c(1, meanScale, effectScale)
                     , posterior = F
                     , method = "auto"
                     , iterations = max(keep))
  bf.one <- nWayAOV(y
                    , prep$X.one
                    , prep$gMap.one
                    , rscale = c(1, meanScale)
                    , posterior = F
                    , method = "auto"
                    , iterations = max(keep))
  mcmc.one <- nWayAOV(y
                      , prep$X.one
                      , prep$gMap.one
                      , rscale = c(1, meanScale)
                      , posterior = T
                      , method = "auto"
                      , iterations = max(keep))
  bf.null <- nWayAOV(y
                       , prep$X.null
                       , prep$gMap.null
                       , rscale = 1
                       , posterior = F
                       , method = "auto"
                     , iterations = max(keep))
  
  ##Positive-effects model encompassing approach
  
  i.theta0 <- prep$I + 2
  i.theta <- (prep$I + 3):(2 * prep$I + 2)
  
  myTheta <- mcmc.full[keep, i.theta] + mcmc.full[keep, i.theta0]
  good <- myTheta > 0
  all.good <- apply(good, 1, mean)
  PostCount <- mean(all.good == 1)
  
  #prior settings
  R <- max(keep) * 10
  beta <- .5 * effectScale^2
  alpha <- .5
  mu.theta.sd <- .5 * meanScale
  
  x <- matrix(nrow = R, ncol = prep$I)
  s2 <- rinvgamma(R, alpha, beta)
  mu <- rcauchy(R, 0, mu.theta.sd)
  for (i in 1:prep$I)
    x[,i] <- rnorm(R, mu, sqrt(s2))
  
  all.greater <- function(x) as.integer(mean(x > 0) == 1)
  PriorCount <- mean(apply(x, 1, all.greater))
  
  ##Positive Common-effect model encompassing approach

  myTheta1 <- mcmc.one[keep, i.theta0]
  PostCount.1 <- mean(myTheta1 > 0)
  
  #prior settings
  mu <- rnorm(R, 0, mu.theta.sd)
  PriorCount.1 <- mean(mu > 0)
  
  bf.FP <- PriorCount/PostCount
  bf.F0 <- exp(bf.full$bf - bf.null$bf)
  bf.F1 <- exp(bf.full$bf - bf.one$bf)
  bf.11p <- PriorCount.1/PostCount.1
  bf.F1p <- bf.F1 * bf.11p
  
  m <- apply(myTheta, 2, mean)
  new.sd <- sd(m)
  new.mean <- mean(m)
  return(list(mean = new.mean, sd = new.sd
              , bf.1f = 1/ bf.F1, bf.pf =  1/ bf.FP
              , bf.0f = 1 / bf.F0, bf.1pF = 1/bf.F1p
              , bf.full = bf.full, bf.one = bf.one, bf.null = bf.null
              , est.full = mcmc.full, est.one = mcmc.one
              , prior.c = PriorCount, post.c = PostCount
              , prior.c.1 = PriorCount.1, post.c.1 = PostCount.1
              , myTheta = myTheta
              , bfs = c(bf.1f = 1/ bf.F1, bf.pf =  1/ bf.FP
                        , bf.0f = 1 / bf.F0, bf.1pF = 1/bf.F1p)))
}

ez.mix <- function(y, sub, cond, priors, keep){
  M <- max(keep)
  keep.int <- 10
  I <- length(unique(sub))
  ybar <- tapply(y, list(sub, cond), mean)
  effect <- ybar[, 2] - ybar[, 1]
  
  mu.alpha <- 1:M
  alpha <- matrix(nrow = M, ncol = I)
  mu.theta <- 0:(M - 1)
  theta <- matrix(nrow = M, ncol = I)
  
  z <- matrix(nrow = M, ncol = I)
  rho <- 1:M
  
  s2 <- 1:M
  g.alpha <- 1:M
  g.theta <- 1:M
  g.mu.th <- 1:M
  
  rho.a <- 1
  rho.b <- 1
  r.alpha <- priors[1]
  r.theta <- priors[2]
  r.mu.th <- priors[3]
  
  alpha[1, ] <- ybar[, 1] - mean(ybar[, 1])
  theta[1, ] <- effect
  rho[1] <- rho.a/sum(rho.a, rho.b)
  z[1, ] <- ifelse(effect < 0, 0, 1)
  g.theta[1] <- r.theta
  g.mu.th[1] <- r.mu.th
  
  sub.i <- sub[cond ==1]
  
  sd.cand <- .015
  count.decor <- 0
  keep.div <- seq(2,min(keep), keep.int)[-1]
  
  #chain
  for(m in 2:M){
    slab <- as.logical(z[m-1, ])
    
    #mu.alpha
    Y <- y - alpha[m-1, sub] - cond * (theta[m-1, ] * slab)[sub]
    mu.alpha[m] <- rnorm(1, mean(Y), sqrt(s2[m-1] / N)) #flat prior
    
    #alpha
    Y <- y - mu.alpha[m] - cond * (theta[m-1, ] * slab)[sub]
    c <- tapply(Y, sub, sum) / s2[m-1]
    v <- 1/(tapply(Y, sub, length) / s2[m-1] + 1/(g.alpha[m-1] * s2[m-1]))
    alpha[m, ] <- rnorm(I, c * v, sqrt(v))
    
    #Decorrelating step for mu and alpha, see Morey et al., 2008
    cand.add <- rnorm(1, 0, sd.cand)
    u <- exp(-1/2 * 
               sum((alpha[m,] + cand.add)^2 - alpha[m, ]^2)/
               (g.alpha[m-1] * s2[m-1]))
    if(rbinom(1, 1, min(u, 1)) == 1){
      alpha[m, ] <- alpha[m,] + cand.add
      mu.alpha[m] <- mu.alpha[m] - cand.add
      count.decor <- count.decor + 1
    }
    if(m %in% keep.div){
      acc.prob <- count.decor/keep.int
      if(acc.prob > .5){
        sd.cand <- sd.cand * 1.25
      }
      if(acc.prob < .25){
        sd.cand <- sd.cand * .75
      }
      count.decor <- 0
    }
    
    #theta
    Y <- y - mu.alpha[m] - alpha[m, sub]
    
    c <- tapply(Y, list(sub, cond), sum)[, 2] / s2[m-1] + mu.theta[m-1] / (g.theta[m-1] * s2[m-1])
    v <- 1/(tapply(Y, list(sub, cond), length)[, 2] / s2[m-1] + 1 / (g.theta[m-1] * s2[m-1]))
    theta[m, slab] <- rnorm(sum(slab), c[slab] * v[slab], sqrt(v[slab]))
    theta[m, !slab] <- rnorm(sum(!slab), mu.theta[m-1], sqrt(g.theta[m-1] * s2[m-1]))
    
    #s2
    Error <- y - mu.alpha[m] - alpha[m, sub] - cond * (theta[m, ] * slab)[sub]
    scale <- sum(Error^2)/2 + 
      sum(alpha[m, ]^2)/(2*g.alpha[m-1]) + 
      sum((theta[m, ] - mu.theta[m-1])^2)/(2*g.theta[m-1]) + 
      mu.theta[m-1]^2/(2 * g.mu.th[m-1])
    s2[m] <- rinvgamma(1, (N + I + I + 1)/2, scale)
    
    #z
    dens.spike <- tapply(dnorm(Y[cond == 1], 0, sqrt(s2[m]), log = T)
                         , sub[cond ==1], sum)
    dens.slab <- tapply(dnorm(Y[cond == 1], theta[m, sub[cond ==1]], sqrt(s2[m]), log = T)
                        , sub[cond ==1], sum)
    prob <- 1/(1 + (1- rho[m-1]) / rho[m-1] * exp(dens.spike - dens.slab))
    # prob[is.na(prob)] <- 0
    z[m,] <- rbinom(I, 1, prob)
    
    #rho
    rho[m] <- rbeta(1, rho.a + sum(z[m, ]), rho.b + sum(!z[m, ]))
    
    #mu.theta
    c <- sum(theta[m,])/(g.theta[m-1] * s2[m-1])
    v <- 1/(I/(g.theta[m-1] * s2[m-1]) + 1/(g.mu.th[m-1] * s2[m-1]))
    mu.theta[m] <- rnorm(1, c * v, sqrt(v))
    
    #gs
    gscale <- sum(alpha[m, ]^2) / (2*s2[m])
    g.alpha[m] <- rinvgamma(1, .5 + I/2, gscale + r.alpha^2/2)
    
    gscale <- sum((theta[m, ] - mu.theta[m])^2) / (2*s2[m])
    g.theta[m] <- rinvgamma(1, .5 + I/2, gscale + r.theta^2/2)
    
    gscale <- mu.theta[m]^2 / (2*s2[m])
    g.mu.th[m] <- rinvgamma(1, 1, gscale + r.mu.th^2/2)
    
    # percents <- round(seq(1, M, length.out = 20))
    # if(m %in% percents) print(paste0(which(percents == m) * 5, "%"))
  }
  
  Samp.full <- cbind(mu.alpha, alpha, mu.theta, z * theta - mu.theta, s2, g.alpha, g.mu.th, g.theta)
  pm.full <- colMeans(Samp.full[keep,])
  rho.full <- rho
  z.full <- z
  theta.full <- theta
  
  #prior calculation
  Mprior <- 100000
  
  #Prior all slab or all spike
  prior.rho <- rbeta(Mprior, rho.a, rho.b)
  prior.z <- rbinom(Mprior, I, prior.rho)
  prior.slab <- mean(prior.z == I)
  prior.spike <- mean(prior.z == 0)
  
  #Prior all slab and all positive
  p.g.mu <- rinvgamma(Mprior, .5, .5 * r.mu.th^2)
  p.mu <- rnorm(Mprior, 0, sqrt(p.g.mu))
  g <- rinvgamma(Mprior, .5, .5 * r.theta^2)
  count <- 1:Mprior
  p.theta <- matrix(ncol = I, nrow = Mprior)
  for (m in 1:Mprior) p.theta[m,] <- rnorm(I, p.mu[m], sqrt(g[m]))
  count <- rowSums(p.theta > 0) == I
  prior.pos.slab <- mean(prior.z == I & count)
  prior.pos <- mean(count)
  
  #Bayes factors
  all.pos <- rowSums(theta.full[keep,] > 0)
  all.slab <- mean(rowSums(z.full[keep,]) == I)
  all.spike <- mean(rowSums(z.full[keep,]) == 0)
  all.pos.slab <- mean(all.pos == I & rowSums(z.full[keep,]) == I)
  all.pos.ss <- mean(all.pos == I)
  
  #BF
  bf.spike <- all.spike / prior.spike
  bf.slab <- all.slab / prior.slab
  bf.pos.slab <- all.pos.slab / prior.pos.slab
  bf.pos.ss <- all.pos.ss/prior.pos
  
  return(list(
    keep = keep
    , sample = Samp.full
    # , theta = theta.full
    , z = colMeans(z.full[keep,])
    , rho = rho
    , counter = count.decor
    , bf0s = bf.spike
    , bffs = bf.slab
    , bfps = bf.pos.slab
    , bfpss = bf.pos.ss
    , prior.prob.spike = prior.spike
    , prior.prob.slab = prior.slab
    , prior.prob.pos.slab = prior.pos.slab
    , post.prob.spike = all.spike
    , post.prob.slab = all.slab
    , post.prob.pos.slab = all.pos.slab
  ))
  
}


all.comp <- function(y, sub, cond, priors, keep){
  M <- max(keep)
  keep.int <- 10
  I <- length(unique(sub))
  ybar <- tapply(y, list(sub, cond), mean)
  effect <- ybar[, 2] - ybar[, 1]
  
  cond.pack <- ifelse(cond == 0, 1, 2)
  prep <- prep.models(sub, cond.pack)
  bf <- makeBF(y, priors[3], priors[2], prep = prep, keep = keep)
  
  mu.alpha <- 1:M
  alpha <- matrix(nrow = M, ncol = I)
  mu.theta <- 0:(M - 1)
  theta <- matrix(nrow = M, ncol = I)
  
  z <- matrix(nrow = M, ncol = I)
  rho <- 1:M
  
  s2 <- 1:M
  g.alpha <- 1:M
  g.theta <- 1:M
  g.mu.th <- 1:M
  
  rho.a <- 1
  rho.b <- 1
  r.alpha <- priors[1]
  r.theta <- priors[3]
  r.mu.th <- priors[2]
  
  alpha[1, ] <- ybar[, 1] - mean(ybar[, 1])
  theta[1, ] <- effect
  rho[1] <- rho.a/sum(rho.a, rho.b)
  z[1, ] <- ifelse(effect < 0, 0, 1)
  g.theta[1] <- r.theta
  g.mu.th[1] <- r.mu.th
  
  sub.i <- sub[cond ==1]
  
  sd.cand <- .015
  count.decor <- 0
  keep.div <- seq(2,min(keep), keep.int)[-1]
  
  #chain
  for(m in 2:M){
    slab <- as.logical(z[m-1, ])
    
    #mu.alpha
    Y <- y - alpha[m-1, sub] - cond * (theta[m-1, ] * slab)[sub]
    mu.alpha[m] <- rnorm(1, mean(Y), sqrt(s2[m-1] / N)) #flat prior
    
    #alpha
    Y <- y - mu.alpha[m] - cond * (theta[m-1, ] * slab)[sub]
    c <- tapply(Y, sub, sum) / s2[m-1]
    v <- 1/(tapply(Y, sub, length) / s2[m-1] + 1/(g.alpha[m-1] * s2[m-1]))
    alpha[m, ] <- rnorm(I, c * v, sqrt(v))
    
    #Decorrelating step for mu and alpha, see Morey et al., 2008
    cand.add <- rnorm(1, 0, sd.cand)
    u <- exp(-1/2 * 
               sum((alpha[m,] + cand.add)^2 - alpha[m, ]^2)/
               (g.alpha[m-1] * s2[m-1]))
    if(rbinom(1, 1, min(u, 1)) == 1){
      alpha[m, ] <- alpha[m,] + cand.add
      mu.alpha[m] <- mu.alpha[m] - cand.add
      count.decor <- count.decor + 1
    }
    if(m %in% keep.div){
      acc.prob <- count.decor/keep.int
      if(acc.prob > .5){
        sd.cand <- sd.cand * 1.25
      }
      if(acc.prob < .25){
        sd.cand <- sd.cand * .75
      }
      count.decor <- 0
    }
    
    #theta
    Y <- y - mu.alpha[m] - alpha[m, sub]
    
    c <- tapply(Y, list(sub, cond), sum)[, 2] / s2[m-1] + mu.theta[m-1] / (g.theta[m-1] * s2[m-1])
    v <- 1/(tapply(Y, list(sub, cond), length)[, 2] / s2[m-1] + 1 / (g.theta[m-1] * s2[m-1]))
    theta[m, slab] <- rnorm(sum(slab), c[slab] * v[slab], sqrt(v[slab]))
    theta[m, !slab] <- rnorm(sum(!slab), mu.theta[m-1], sqrt(g.theta[m-1] * s2[m-1]))
    
    #s2
    Error <- y - mu.alpha[m] - alpha[m, sub] - cond * (theta[m, ] * slab)[sub]
    scale <- sum(Error^2)/2 + 
      sum(alpha[m, ]^2)/(2*g.alpha[m-1]) + 
      sum((theta[m, ] - mu.theta[m-1])^2)/(2*g.theta[m-1]) + 
      mu.theta[m-1]^2/(2 * g.mu.th[m-1])
    s2[m] <- rinvgamma(1, (N + I + I + 1)/2, scale)
    
    #z
    dens.spike <- tapply(dnorm(Y[cond == 1], 0, sqrt(s2[m]), log = T)
                         , sub[cond ==1], sum)
    dens.slab <- tapply(dnorm(Y[cond == 1], theta[m, sub[cond ==1]], sqrt(s2[m]), log = T)
                        , sub[cond ==1], sum)
    prob <- 1/(1 + (1- rho[m-1]) / rho[m-1] * exp(dens.spike - dens.slab))
    # prob[is.na(prob)] <- 0
    z[m,] <- rbinom(I, 1, prob)
    
    #rho
    rho[m] <- rbeta(1, rho.a + sum(z[m, ]), rho.b + sum(!z[m, ]))
    
    #mu.theta
    c <- sum(theta[m,])/(g.theta[m-1] * s2[m-1])
    v <- 1/(I/(g.theta[m-1] * s2[m-1]) + 1/(g.mu.th[m-1] * s2[m-1]))
    mu.theta[m] <- rnorm(1, c * v, sqrt(v))
    
    #gs
    gscale <- sum(alpha[m, ]^2) / (2*s2[m])
    g.alpha[m] <- rinvgamma(1, .5 + I/2, gscale + r.alpha^2/2)
    
    gscale <- sum((theta[m, ] - mu.theta[m])^2) / (2*s2[m])
    g.theta[m] <- rinvgamma(1, .5 + I/2, gscale + r.theta^2/2)
    
    gscale <- mu.theta[m]^2 / (2*s2[m])
    g.mu.th[m] <- rinvgamma(1, 1, gscale + r.mu.th^2/2)
    
    # percents <- round(seq(1, M, length.out = 20))
    # if(m %in% percents) print(paste0(which(percents == m) * 5, "%"))
  }
  
  Samp.full <- cbind(mu.alpha, alpha, mu.theta, z * theta - mu.theta, s2, g.alpha, g.mu.th, g.theta)
  pm.full <- colMeans(Samp.full[keep,])
  rho.full <- rho
  z.full <- z
  theta.full <- theta
  
  #prior calculation
  Mprior <- 100000
  
  #Prior all slab or all spike
  prior.rho <- rbeta(Mprior, rho.a, rho.b)
  prior.z <- rbinom(Mprior, I, prior.rho)
  prior.slab <- mean(prior.z == I)
  prior.spike <- mean(prior.z == 0)
  
  #Prior all slab and all positive
  p.g.mu <- rinvgamma(Mprior, .5, .5 * r.mu.th^2)
  p.mu <- rnorm(Mprior, 0, sqrt(p.g.mu))
  g <- rinvgamma(Mprior, .5, .5 * r.theta^2)
  count <- 1:Mprior
  p.theta <- matrix(ncol = I, nrow = Mprior)
  for (m in 1:Mprior) p.theta[m,] <- rnorm(I, p.mu[m], sqrt(g[m]))
  count <- rowSums(p.theta > 0) == I
  prior.pos.slab <- mean(prior.z == I & count)
  prior.pos <- mean(count)
  
  #Bayes factors
  all.pos <- rowSums(theta.full[keep,] > 0)
  all.slab <- mean(rowSums(z.full[keep,]) == I)
  all.spike <- mean(rowSums(z.full[keep,]) == 0)
  all.pos.slab <- mean(all.pos == I & rowSums(z.full[keep,]) == I)
  all.pos.ss <- mean(all.pos == I)
  
  #BF
  bf.spike <- all.spike / prior.spike
  bf.slab <- all.slab / prior.slab
  bf.pos.slab <- all.pos.slab / prior.pos.slab
  bf.pos.ss <- all.pos.ss/prior.pos
  
  return(list(
      bf0s = bf.spike
      , bffs = bf.slab
      , bfps = bf.pos.slab
      , bfpss = bf.pos.ss
      , bfpf = bf$bf.pf
      , bf1f = bf$bf.1f
      , bf1pf = bf$bf.1pF
      , bf0f = bf$bf.0f
      , prior.prob.spike = prior.spike
      , prior.prob.slab = prior.slab
      , prior.prob.pos.slab = prior.pos.slab
      , post.prob.spike = all.spike
      , post.prob.slab = all.slab
      , post.prob.pos.slab = all.pos.slab
    ))
}
