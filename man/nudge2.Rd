\name{nudge2}
\alias{nudge2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for normalizing data, fitting a normal-uniform mixture and estimating probabilities of differential expression in the case where the two samples are being compared indirectly through a common reference sample}
\description{
After a mean and variance normalization a two component mixture model is fitted to the data. The normal component represents the genes that are not differentially expressed and the uniform component represents the genes that are differentially expressed. Posterior probabilities for differential expression are computed from the fitted model.
}
\usage{
nudge2(control.logratio, txt.logratio, control.logintensity, txt.logintensity,
span1 = 0.2, quant = 0.99, z = NULL, tol = 0.00001,iterlim=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{control.logratio}{A multiple-column matrix of replicates of log (base 2) ratios of gene expressions for the control versus reference slides.}
  \item{txt.logratio}{A multiple-column matrix of replicates of log (base 2) ratios of gene expressions for the treatment versus reference slides.}
  \item{control.logintensity}{A multiple-column matrix of replicates of log (base 2) total intensities (defined as the product) of gene expressions for the control versus reference slides.}
  \item{txt.logintensity}{A multiple-column matrix of replicates of log (base 2) total intensities (defined as the product) of gene expressions for the treatment versus reference slides.}
  \item{span1}{Proportion of data used to fit the loess regression of the (average-across-replicates) log ratio differences on the (average-across-replicates) log intensities for the mean normalization.}
  \item{quant}{Quantile to be used from the distribution of standard deviations of log ratio differences across replicates for all genes whose standard deviation was smaller than their absolute (mean normalized) average-across-replicates log ratio difference.}
  \item{z}{An optional 2-column matrix with each row giving a starting estimate for the probability of the gene (in the corresponding row of the log ratio matrix/vector) not being differentially expressed and a starting estimate for the probability of the gene being differentially expressed. Each row should add up to 1.}
  \item{tol}{A scalar tolerance for relative convergence of the loglikelihood.}
  \item{iterlim}{The maximum number of iterations the EM is run for.}
}
\value{
A list including the following components
  \item{pdiff }{A vector with the estimated posterior probabilities of being in the group of differentially expressed genes.}
  \item{lRnorm }{A vector with the normalized (average-across-replicates) log ratio differences.}
  \item{mu }{The estimated mean of the group of genes that are not differentially expressed.}
  \item{sigma }{The estimated variance of the group of genes that are not differentially expressed.}
  \item{mixprob}{The prior/mixing probability of a gene being in the group of genes that are not differentially expressed.}
  \item{a }{The minimum value of the normalized data.}
  \item{b }{The maximum value of the normalized data.}
  \item{loglike }{The log likelihood for the fitted mixture model.}
  \item{iter }{The number of iterations run by the EM algorithm until either convergence or iteration limit was reached.}
}
\references{N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays.  BMC Bioinformatics. 6, 173-186. 

\url{http://www.biomedcentral.com/1471-2105/6/173}

S. Dudoit, Y. H. Yang, M. Callow and T. Speed (2002). Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. Stat. Sin. 12, 111-139. 
}
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{nudge1}},\code{\link{norm2c}},\code{\link{norm2d}},\code{\link{norm1a}},\code{\link{norm1b}},\code{\link{norm1c}},\code{\link{norm1d}}}
\examples{
apo<-read.csv("http://www.stat.berkeley.edu/users/terry/zarray/Data/ApoA1/rg_a1ko_morph.txt",
header=TRUE)
rownames(apo)<-apo[,1]
apo<-apo[,-1]
apo<-apo+1

lRctl<-log(apo[,c(seq(2,16,2))],2)-log(apo[,c(seq(1,15,2))],2)
lRtxt<-log(apo[,c(seq(18,32,2))],2)-log(apo[,c(seq(17,31,2))],2)
lIctl<-log(apo[,c(seq(2,16,2))],2)+log(apo[,c(seq(1,15,2))],2)
lItxt<-log(apo[,c(seq(18,32,2))],2)+log(apo[,c(seq(17,31,2))],2)
 
result<-nudge2(lRctl,lRtxt,lIctl,lItxt)

}

\keyword{cluster}
\keyword{loess}
