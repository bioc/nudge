"norm1a" <-
function(logratio,logintensity,span=0.6)
{
 r<-require(stats)
 if(!r){require(modreg)}
 logratio-loess(logratio~logintensity,span=span)$fitted
}

