"norm2c" <-
function(control.logratio,txt.logratio,control.logintensity,txt.logintensity,span=0.6)
{
 r<-require(stats)
 if(!r){require(modreg)}
 lRbar<-apply(txt.logratio,1,mean)-apply(control.logratio,1,mean)
 lIbar<-apply(cbind(txt.logintensity,control.logintensity),1,mean)
 lRbar-loess(lRbar~lIbar,span=span)$fitted
}

