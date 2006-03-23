"nudge2" <-
function(control.logratio,txt.logratio,control.logintensity,txt.logintensity,span1=0.2,quant=0.99,z=NULL,tol=0.00001,iterlim=500)
{
#nudge2 is the nudge algorithm when the samples (control versus treatment say) are labelled by same dye in different (either technical or biological) replicates and the reference sample for all is labelled by the other dye (e.g. Apo AI data in paper by Sandrine Dudoit, Yee Hwa Yang, Matt Callow and Terry Speed).

#First check that the data for ratio and intensity are of same size
l<-length(control.logratio)
if(l!=length(control.logintensity)) stop("Control log ratios and log intensities not of same size")
if(length(txt.logintensity)!=length(txt.logratio)) stop("Treatment log ratios and log intensities not of same size")

#Check if data is from replicate or single or single experiments
d<-dim(control.logratio)
#if d is NULL then object is a vector and not a matrix
if(is.null(d))
{
 stop("There is only one replicate, recommend using nudge1 function instead")
} else{
n<-ncol(control.logratio)

#data could still be a single replicate if n=1
if(n==1)
{
 stop("There is only one replicate, recommend using nudge1 function instead")
}else{
if((nrow(control.logratio))!=(nrow(txt.logratio)))
{
stop("Control and treatment number of genes/rows are different")
}else{
if((nrow(control.logratio))!=(nrow(control.logintensity)))
{
stop("Control number of genes'/rows' of log ratios and intensities are different")
}else{
if((nrow(txt.logratio))!=(nrow(txt.logintensity)))
{
stop("Treatment number of genes'/rows' of log ratios and intensities are different")
}else{

#for replicate data
 lRnorm<-norm2d(control.logratio,txt.logratio,control.logintensity,txt.logintensity,span1,quant)
}}}}}

 X<-lRnorm

 #EM algorithm for univariate normal and uniform mixture
 n<-length(X)

#z is matrix for probabilities of outliers, column 1 is probability of not being an outlier, column 2 (=1 - column 1) is probability of being an outlier
if(is.null(z))
{
 z<-matrix(0,n,2)

#Getting starting values for z, using usual rule-of-two cutoff on normalised data
 m<-mean(X)
 s<-sqrt(var(X))
 d<-abs((X-m)/s)>2
 z[,2]<-d^2
 z[,1]<-1-z[,2]
}

#Compute first estimates for mixing probabilities from initial z
 p<-sum(z[,1])/n
 muhat<-sum(z[,1]*X)/sum(z[,1])
 sigma2hat<-sum(z[,1]*(X-muhat)^2)/sum(z[,1])
 sigmahat<-sqrt(sigma2hat)
 llike<-c(0,100)
 criterion<-abs(llike[1]-llike[2])
 iter<-0
#tol is the critical value where the algorithm stops if the values of two consqutive log-likelihoods are within tol of each other

 while((criterion>tol)&(iter<iterlim))
#iteration limit to ensure if no convergence for some reason, algorithm doesn't run for ever
 {
  iter<-iter+1
  #E Step-estimating the z's
  z[,1]<-(p*dnorm(X,muhat,sigmahat))/((p*dnorm(X,muhat,sigmahat))+((1-p)*dunif(X,min(X),max(X))))
  z[,2]<-1-z[,1]

  #M Step-estimating parameters
  p<-sum(z[,1])/n
  muhat<-sum(z[,1]*X)/sum(z[,1])
  sigma2hat<-sum(z[,1]*(X-muhat)^2)/sum(z[,1])
  sigmahat<-sqrt(sigma2hat)

  loglike<-sum(log((p*dnorm(X,muhat,sigmahat))+((1-p)*dunif(X,min(X),max(X)))))
  llike[2]<-llike[1]
  llike[1]<-loglike
#calculate absolute difference between log-likelihood for estimated parameters and previous log-likelihood
  criterion<-abs(llike[1]-llike[2])
 }
 colnames(z)<-c("Prob. non-diff. exp.","Prob diff. exp.")
 list(pdiff=z[,2],lRnorm=X,mu=muhat,sigma=sigmahat,mixprob=p,a=min(X),b=min(X),loglike=loglike,iter=iter)
}

