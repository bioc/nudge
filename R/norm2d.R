"norm2d" <-
function(control.logratio,txt.logratio,control.logintensity,txt.logintensity,span=0.6,quant=0.99)
{
 n<-nrow(txt.logratio)
 d1<-ncol(txt.logratio)
 d2<-ncol(control.logratio)
 lRnorm<-norm2c(control.logratio,txt.logratio,control.logintensity,txt.logintensity,span)
 var.txt<-apply(txt.logratio,1,var)/d1
 var.control<-apply(control.logratio,1,var)/d2
 sd<-sqrt(var.txt+var.control)
 k<-quantile(sd[sd<=abs(lRnorm)],quant)
 sdnew<-sd
 change<-c(1:n)[sd<=abs(lRnorm)]
 for(i in change)
 {
  sdnew[i]<-k
 }
 lRnorm/sdnew
}

