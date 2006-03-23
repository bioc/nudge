"norm1c" <-
function(logratio,logintensity,span=0.6)
{
 r<-require(stats)
 if(!r){require(modreg)}
 lRbar<-apply(logratio,1,mean)
 lIbar<-apply(logintensity,1,mean)
 lRbar-loess(lRbar~lIbar,span=span)$fitted
}

