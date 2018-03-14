#####This script is used for a simulation study. I want to know how well the unconstrained model can be used to detect mixtures
####I am simulating from three truths: Two Mixtures and One graded model.

# library("msm")
# library("gridBase")
# library("ggplot2")
# library("grid")
# library("mvtnorm")
# library("tmvtnorm")
# library('BayesFactor')
# library('MCMCpack')

I <- 30
K <- 50
J <- 2
N <- I*J*K

sub <- rep(1:I, each=J*K)
trial <- rep(1:(J*K), I)
cond <- rep((1:J) - 1, each = K, I) #cond = 1 for incongruent, =0 for congruent

set.seed(124)

t.alpha <- runif(I, .75, .75) #Slightly variable baseline
t.sigma <- .4 #True sample noise
rscale <- c(1, 1/6, 1/10) #scales for alpha, mu_theta, theta
# source("modelcomp.R")

#Function to simulate data from different true values
makeData <- function(t.theta, alpha = t.alpha, sigma = t.sigma
                     , R = N, condition = cond, subject = sub){
	y <- rtnorm(R, alpha[subject] + condition * t.theta[subject], sigma, lower = 0)
	dat<- data.frame(sub = subject, cond = condition, y)
	return(dat)
}

###Function to compute Bayes factors for simulated data
bayesBF <- function(dat, rScale = rscale, M = 10000)
{
  dat$cond <- dat$cond + 1
  prep.bf <- prep.models(dat$sub, dat$cond)
  makeBF(dat$y, meanScale = rScale[2], effectScale = rScale[3]
         , prep = prep.bf, keep = 1001:M)$bfs
}

#Simulation wrapper
make.sim <- function(tData){ #input true data
  dat <- makeData(tData)
  # print(dat$y[1])
  return(bayesBF(dat))
}

#Helper functions
DIM <- function( ... ){
  args <- list(...)
  lapply( args , function(x) { if( is.null( dim(x) ) )
    return( length(x) )
    dim(x) } )
}

sameLength <- function(...){
  x <- unlist(DIM(...))
  var(x) == 0
}

#sample from mixture
rmixture <- function(n, prob, mu, sigma, spike = F, low, up){
  if(spike == T){
    mu <- c(0, mu)
    sigma <- diag(dim(sigma)[1] + 1) * c(10e-20, diag(sigma))
  }
  if(sameLength(prob, mu, sigma, low, up) == F){
    stop("Parameters have wrong dimensions.")}
    return(apply(matrix(1:n), 1, function(x){rtmvnorm(1, mu, sigma, lower = low, upper = up) %*% rmultinom(1, 1, prob)}))
}


# 3 simulations,
#normal 
#mixture of a negative and a positive distribution
#mixtrue of a negative, a spike at zero, and a positive distribution
#Check plot below for true data


#true values:
set.seed(123)
normPar <- c(.06, .08)
mixPar1 <- list(
  prob = c(.22, .78)
  , mu = c(-.06, .08)
  , sigma = c(.03, .06)
  , low = c(-1, -Inf, 0)
  , up = c(1, 0, Inf)
)
mixPar2 <- list(
  prob = c(.15, .2, .65) #prob of spike is first
  , mu = c(-.05, .1)
  , sigma = c(.03, .06)
  , low = c(-1, -Inf, 0)
  , up = c(1, 0, Inf)
)

tNorm <- qnorm((1:I)/(I+1), normPar[1], normPar[2])

r1Mix <- rmixture(n = 10000, prob = mixPar1$prob
                  , mu = mixPar1$mu, sigma = diag(2) * mixPar1$sigma^2
                  , spike = F
                  , low = mixPar1$low[-1]
                  , up = mixPar1$up[-1])
t1Mix <- quantile(r1Mix, (1:I)/(I+1))

r2Mix <- rmixture(n = 10000, prob = mixPar2$prob
                  , mu = mixPar2$mu, sigma = diag(2) * mixPar2$sigma^2
                  , spike = T
                  , low = mixPar2$low
                  , up = mixPar2$up)
t2Mix <- quantile(r2Mix, (1:I)/(I+1))

tDat <- data.frame(tNorm, t1Mix, t2Mix)

# plot(tNorm
#      , ylim = c(-.25, .3)
#      )
# points(t1Mix, col = "darkred")
# points(t2Mix, col = "darkgreen")
# abline(h=0)
meansSim <- apply(tDat, 2, mean)
sdsSim <- apply(tDat, 2, sd)

#sims
R <- 100 #Number of simulations per model
nModels <- 3

sims <- rep(1:R, nModels)
mods <- rep(1:nModels, each = R)

simResults <- sapply(mods, function(x){
  message(x)
  make.sim(tDat[, x])
  })
simRes <- data.frame(cbind(mods, sims, t(simResults)))
simRes$mods <- factor(mods, labels = c("Normal", "Mixture 1", "Mixture 2"))
