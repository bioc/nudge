"norm1b" <-
function(logratio,logintensity,span1=0.6,span2=0.2,mean.norm=TRUE)
{
 if(mean.norm){
 lRnorm1<-norm1a(logratio,logintensity,span=span1)
 }else{
 lRnorm1<-logratio
 }
 l<-lRnorm1/loess(abs(lRnorm1)~logintensity,span=span2)$fitted
#Get data on same scale as original
 l*median(abs(lRnorm1))/median(abs(l))
}

