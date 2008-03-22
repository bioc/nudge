"norm1d" <-
function(logratio,logintensity,span=0.6,quant=0.99,dye.swap=FALSE)
{
 n<-nrow(logratio)
 d<-ncol(logratio)
 if(d<2||is.null(d)){stop("Single not multiple replicate data entered")}
 ifelse(!(dye.swap),lRnorm<-norm1c(logratio,logintensity,span),lRnorm<-apply(logratio,1,mean))
 sd<-sqrt(apply(logratio,1,var)/d)
 k<-quantile(sd[sd<=abs(lRnorm)],quant)
 sdnew<-sd
 change<-c(1:n)[sd<=abs(lRnorm)]
 for(i in change)
 {
  sdnew[i]<-k
 }
 lRnorm/sdnew
}

