##Positive Slab help funtions

##Sample mean with metropolis hastings step
sampleNu <- function(nu.0, muSD, mu, m.nu, g.nu, counter.nu, g.0, s2){
  #sample nu
  nu <- nu.0
  cand <- rnorm(1, nu, muSD)
  
  if(cand > 0){
    full.cur <- sum(dnorm(mu, nu, sqrt(g.0 * s2), log = T)) + dnorm(nu, m.nu, sqrt(g.nu * s2), log = T) - I * pnorm(0, -nu, sqrt(g.0 * s2), log = T)
    
    full.cand <- sum(dnorm(mu, cand, sqrt(g.0 * s2), log = T)) + dnorm(cand, m.nu, sqrt(g.nu * s2), log = T) - I * pnorm(0, -cand, sqrt(g.0 * s2), log = T)
    
    prob <- min(exp(full.cand - full.cur), 1)
    
    if (rbinom(1, 1, prob) == 1){
      nu <- cand
      counter.nu <- counter.nu + 1
    }
  }
  
  return(c(nu, counter.nu))
}

##Sample mean with metropolis hastings step with nu not restricted to be greater than 0
sampleNuFree <- function(nu.0, muSD, mu, m.nu, g.nu, counter.nu, g.0, s2){
  #sample nu
  nu <- nu.0
  cand <- rnorm(1, nu, muSD)
  
  # if(cand > 0){
    full.cur <- sum(dnorm(mu, nu, sqrt(g.0 * s2), log = T)) + dnorm(nu, m.nu, sqrt(g.nu * s2), log = T) - I * pnorm(0, -nu, sqrt(g.0 * s2), log = T)
    
    full.cand <- sum(dnorm(mu, cand, sqrt(g.0 * s2), log = T)) + dnorm(cand, m.nu, sqrt(g.nu * s2), log = T) - I * pnorm(0, -cand, sqrt(g.0 * s2), log = T)
    
    prob <- min(exp(full.cand - full.cur), 1)
    
    if (rbinom(1, 1, prob) == 1){
      nu <- cand
      counter.nu <- counter.nu + 1
    }
  # }
  
  return(c(nu, counter.nu))
}

#sample variance with metropolis hastings step
sampleEta <- function(nu.0, mu, g.0, varSD, r, counter.eta, s2){
  #sample eta
  g <- g.0
  cand <- rnorm(1, g, varSD)
  
  if(cand > 0){
    full.cur <- sum(dnorm(mu, nu.0, sqrt(g * s2), log = T)) + log(dinvgamma(g, 1/2, r^2/2)) - I * pnorm(0, -nu.0, sqrt(g * s2), log = T)
    
    full.cand <- sum(dnorm(mu, nu.0, sqrt(cand * s2), log = T)) + log(dinvgamma(cand, 1/2, r^2/2)) - I * pnorm(0, -nu.0, sqrt(cand * s2), log = T)
    
    prob <- min(exp(full.cand - full.cur), 1)
    
    if (rbinom(1, 1, prob) == 1){
      g <- cand
      counter.eta <- counter.eta + 1
    }
  }
  return(c(g, counter.eta))
}

#function to get the log of the posterior distribution of s2 if theta is not truncated
logfS2 <- function(s2, shape, rate){
  - (shape + 1) * log(s2) - rate / s2
}

make.rate <- function(y, mum, alpham, thetam
                      , g.alpham, g.thetam, g.mu.thetam, mu.thetam
                      , subm, condm, slabm){
  Error <- y - mum - alpham[subm] - condm * (slabm * thetam)[subm]
  
  sum(Error^2)/2 + 
    sum(alpham^2) / (2 * g.alpham) + 
    sum((thetam - mu.thetam)^2) / (2 * g.thetam) + 
    mu.thetam^2 / (2 * g.mu.thetam)
}

#sample variance with metropolis hastings step
sampleS2 <- function(curS2, shape, rate, muTheta, gTheta, SD, counter){
  #sample s2
  s2 <- curS2
  cand <- rnorm(1, s2, SD)
  
  if(cand > 0){
    full.cur <- logfS2(s2, shape, rate) - I * pnorm(0, -muTheta, sqrt(gTheta * s2), log = T)
    
    full.cand <- logfS2(cand, shape, rate) - I * pnorm(0, -muTheta, sqrt(gTheta * cand), log = T)
    
    prob <- min(exp(full.cand - full.cur), 1)
    
    if (rbinom(1, 1, prob) == 1){
      s2 <- cand
      counter <- counter + 1
    }
  }
  return(c(s2, counter))
}


#adjust sampling variance during burn-in
check.scale <- function(iteration, check.nu, check.g, check.gnu, keep, counter.nu, counter.g, counter.gnu, muSD, varSD, s2SD){
  if(iteration == check.nu*100 & iteration < keep){
    if(counter.nu < 30){muSD <<- .5*muSD}
    if(counter.nu > 70){muSD <<- 2*muSD}
    counter.nu <- 0
    check.nu <- check.nu + 1
    
    if(counter.g < 40){varSD <<- .5*varSD}
    if(counter.g > 60){varSD <<- 2*varSD}
    counter.g <- 0
    check.g <- check.g + 1
    
    if(counter.gnu < 40){s2SD <<- .5*s2SD}
    if(counter.gnu > 60){s2SD <<- 2*s2SD}
    counter.gnu <- 0
    check.gnu <- check.gnu + 1
  }
  return(c(counter.nu, counter.g, counter.gnu, check.nu, check.g, check.gnu))
}

##addign multiple values to multiple variables - handy in the chain
let <- function(x, value) {
  mapply(
    assign,
    as.character(substitute(x)[-1]),
    value,
    MoreArgs = list(envir = parent.frame()))
  invisible()
}

counter.it <- function(m, M){
  percents <- round(seq(1, M, length.out = 20))
  if(m %in% percents) paste0(which(percents == m), "%")
  if(m == round(M/2)) paste("Half-way there! Wohoo!")
  if(m == M) paste("Done! Success! Watch this video: https://www.youtube.com/watch?v=taOw0intSWk You've earned it.")
}