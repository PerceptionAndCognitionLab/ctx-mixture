par(mar=c(.1,.1,.1,.1))

# library('diagram')

names=c('Spike-and-Slab','Unconstrained','Positive Effects','Common Effect','Null')

myCol=c('slategray1','white')

pos=matrix(nrow=5,ncol=2)
pos[,1]=c(.25, .75, .5, .75, .25)
pos[,2]=c(.85, .85, .5, .3, .10)

M <- matrix(nrow=length(names), byrow=F, ncol=length(names), data=0)
M[1,3]='E'
M[2,4]=M[4,5]='A'
M[2,3]='E'
M[1, 2] <- 'E'
M[1, 5] <- 'E'

plotmat(M,pos,name=names,curve=0,box.type="rect"
        , box.size=.15,box.prop=.28,box.col=myCol[c(1,1,1,1,1,1)],arr.length=0
        , cex = .9
        , box.cex = .9)
legend("bottomright"
       , legend = c("A = Analytical"
                    , "E = Encompassing")
       , col = "white", bty = "n", cex = 1)
