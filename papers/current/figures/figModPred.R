sd0 <- .07
eta = .7

gamma <- seq(-.25, .25, .0025)

kern <- convKernel(sigma = 8, k = "gaussian")

nrmlz <- function(mat)
{
  tot <- sum(mat)
  mat/tot
}


#Conditional model specification
norm0 <- function(theta1, theta2, Sigma) dnorm(theta1, 0,Sigma) * dnorm(theta2, 0, Sigma)
norm <- function(theta1, theta2, Sigma) dmvnorm(cbind(theta1, theta2), c(0,0), Sigma)
normT1 <- function(theta1, theta2, Sigma, l, u) dtmvnorm(cbind(theta1, theta2)
                                                   , c(0,0)
                                                   , Sigma
                                                   , lower = rep(l, 2)
                                                   , upper = rep(u, 2))
normT <- function(theta1, theta2, Sigma, l , u){
  dtnorm(theta1, 0, Sigma, lower = l, upper = u) * dtnorm(theta2, 0, Sigma, lower = l, upper = u)
}

Null <- outer(gamma, gamma, norm0, Sigma = .002)
Null <- nrmlz(Null)
One <- outer(gamma
                   , gamma
                   , normT1
                   , Sigma = matrix(c(sd0^2, sd0^2.001, sd0^2.001, sd0^2)
                                    , nrow = 2)
                   , l = 0
                   , u = Inf) 
One <- nrmlz(One)
Pos <- outer(gamma
                   , gamma
                   , normT
                   , sd0
                   , l = 0
                   , u = Inf)
Pos <- nrmlz(Pos)
General <- outer(gamma
                 , gamma
                 , norm
                 , Sigma = matrix(c(sd0^2, 0, 0, sd0^2)
                                  , nrow = 2))
General <- nrmlz(General)

priorPos1 <- outer(gamma
                   , gamma
                   , normT1
                   , Sigma = matrix(ncol = 2, c(sd0^2, 0, 0, .005^2))
                   , l = 0
                   , u = Inf)
priorPos1 <- nrmlz(priorPos1)

priorPos2 <- outer(gamma
                   , gamma
                   , normT1
                   , Sigma = matrix(ncol = 2, c(.005^2, 0, 0, sd0^2))
                   , l = 0
                   , u = Inf)
priorPos2 <- nrmlz(priorPos2)

priorSpike <- outer(gamma
                   , gamma
                   , normT1
                   , Sigma = matrix(ncol = 2, c(.005^2, 0, 0, .005^2))
                   , l = 0
                   , u = Inf)
priorSpike <- nrmlz(priorSpike)

Mix <- .01 * priorSpike + .5 * priorPos1 + .5 * priorPos2 + 10 * Pos
Mix <- nrmlz(Mix)

#Marginal model specification
GeneralH <- outer(gamma
                  , gamma
                  , norm
                  , Sigma = matrix(c(sd0^2, eta*sd0^2, eta*sd0^2, sd0^2)
                                     , nrow = 2))
GeneralH <- nrmlz(GeneralH)

PosH <- 4 * GeneralH
index <- gamma < 0
PosH[index, ] <- 0
PosH[, index] <- 0
PosH <- nrmlz(PosH)

MixH <- .01 * priorSpike + .5 * priorPos1 + .5 * priorPos2 + 7 * PosH
MixH <- nrmlz(MixH)

#Model Predictions
NullP <- nrmlz(applyFilter(Null, kern))
OneP <- nrmlz(applyFilter(One, kern))
PosP <- nrmlz(applyFilter(PosH, kern))
GeneralP <- nrmlz(applyFilter(GeneralH, kern))
MixP <- .25 * priorSpike + .5 * priorPos1 + .5 * priorPos2 + PosH
MixP <- nrmlz(applyFilter(MixP, kern))

#####Figure
top1 <- max(One, PosH)
top2 <- max(Pos)
top3 <- max(NullP)

modFig <- function(mat, par, ylabel, xlabel, main, top, mod, xax = TRUE, yax = TRUE){
  image(par
        , par
        , mat
        , col = grey((256:0)/256)
        , zlim = c(0, top)
        , axes = FALSE
        , ylab = ylabel
        , xlab = xlabel
        , frame.plot=TRUE
        , main = ""
        , cex.lab = 1.2)
  if(xax == TRUE){
  axis(1, at = seq(-.2, .2, .2), cex.axis = 1.2)}
  if(yax == TRUE){
  axis(2, at = seq(-.2, .2, .2), cex.axis = 1.2)}
  abline(h = 0, col = "gray70")
  abline(v = 0, col = "gray70")
  mtext(mod, side = 2, line = 4)
  mtext(main, side = 3, line = 1)
}

# pdf('figModPred.pdf',width=10,height=20)
layout(matrix(1:15, ncol = 3), widths = c(.37, .28, .35), heights = c(.235, rep(.53/3, 3), .235))

par(mar=c(1,6,3.5,1), mgp = c(2.4,.9,0), pty = "s")
#models
modFig(Mix, gamma
       , ylabel = expression(paste(theta[2])), xlabel = ""
       , mod = expression("M"["SS"^"+"]), top = top3 - .0005
       , main = "Conditional Models"
       , xax = FALSE)

par(mar=c(1,6,1,1), mgp = c(2.4,.9,0))
modFig(General, gamma
       , ylabel = expression(paste(theta[2])), xlabel = ""
       , mod = expression("M"["u"]), top = top2, main = ""
       , xax = FALSE)

modFig(Pos, gamma
       , ylabel = expression(paste(theta[2])), xlabel = ""
       , mod = expression("M"["+"]), top = top2, main = ""
       , xax = FALSE)

modFig(One, gamma
       , ylabel = expression(paste(theta[2])), xlabel = ""
       , mod = expression("M"[1]), top = top1, main = ""
       , xax = FALSE)

par(mar=c(3.5,6,1,1), mgp = c(2.4,.9,0))
modFig(Null, gamma
       , ylabel = expression(paste(theta[2])), xlabel = expression(paste(theta[1]))
       , top = top1, mod = expression("M"[0]), main = "")
points(0, 0, pch = 19)

#marginal
par(mar=c(1,1,3.5,0), mgp = c(2.4,.9,0))
modFig(MixH, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top3 -.0005
       , main = "Marginal Models"
       , xax = FALSE
       , yax = FALSE)

par(mar=c(1,1,1,0), mgp = c(2.4,.9,0))
modFig(GeneralH, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top2, main = ""
       , xax = FALSE
       , yax = FALSE)

modFig(PosH, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top2, main = ""
       , xax = FALSE
       , yax = FALSE)

modFig(One, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top1, main = ""
       , xax = FALSE
       , yax = FALSE)

par(mar=c(3.5,1,1,0), mgp = c(2.4,.9,0))

modFig(Null, gamma
       , ylabel = "", xlabel = expression(paste(theta[1]))
       , top = top1, mod = "", main = ""
       , yax = FALSE)
points(0, 0, pch = 19)

#predictions
par(mar=c(1,4.5,3.5,1), mgp = c(2.4,.9,0))
modFig(MixP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(theta)[2])), xlabel = ""
       , mod = "", top = top2, main = "Predictions"
       , xax = FALSE)

par(mar=c(1,4.5,1,1), mgp = c(2.4,.9,0))
modFig(GeneralP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(theta)[2])), xlabel = ""
       , mod = "", top = top2, main = ""
       , xax = FALSE)

modFig(PosP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(theta)[2])), xlabel = ""
       , mod = "", top = top2, main = ""
       , xax = FALSE)

modFig(OneP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(theta)[2])), xlabel = ""
       , mod = "", top = top2 + .001, main = ""
       , xax = FALSE)

par(mar=c(3.5,4.5,1,1), mgp = c(2.4,.9,0))

modFig(NullP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(theta)[2]))
       , xlabel = expression(paste("Observed MIC, ", hat(theta)[2]))
       , top = top3, mod = "", main = "")

# dev.off()

##Prediction for -.1 and -.09
# Parallel1P[61, 63]/Parallel2P[61, 63]


