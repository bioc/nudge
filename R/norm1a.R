"norm1a" <-
function(logratio,logintensity,span=0.6)
{
 logratio-loess(logratio~logintensity,span=span)$fitted
}

