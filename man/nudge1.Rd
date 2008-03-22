\name{nudge1}
\alias{nudge1}
\title{Function for normalizing data, fitting a normal-uniform mixture and estimating probabilities of differential expression in the case where the two samples are being compared directly}
\description{
After a mean and variance normalization, a two component mixture model is fitted to the data. The normal component represents the genes that are not differentially expressed and the uniform component represents the genes that are differentially expressed. Posterior probabilities for differential expression are computed from the fitted model.
}
\usage{
nudge1(logratio, logintensity, dye.swap = FALSE, span1 = 0.6, span2 = 0.2,
quant = 0.99, z = NULL, tol = 0.00001,iterlim=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logratio}{A matrix or vector of log (base 2) ratios of intensity expressions in 2 samples, with rows indexing genes and columns (if necessary) indexing replicates.}
  \item{logintensity}{A matrix or vector of total log (base 2) total intensities (defined as the product) of intensity expressions in 2 samples, with rows indexing genes and columns (if necessary) indexing replicates.}
  \item{dye.swap}{A logical value indicating whether or not the data is from a balanced dye-swap. Only used for multiple replicate experiments.}
  \item{span1}{Proportion of data used to fit the loess regression of the (average-across-replicates) log ratios on the (average-across-replicates) log total intensities for the mean normalization.}
  \item{span2}{Proportion of data used to fit the loess regression of the absolute (mean normalized) log ratios on the log total intensities for the variance normalization. Only used for single replicate experiments.}
  \item{quant}{Quantile to be used from the distribution of standard deviations of log ratios across replicates for all genes whose standard deviation was smaller than their absolute (mean normalized) average-across-replicates log ratio. Only used for multiple replicate experiments.}
  \item{z}{An optional 2-column matrix with each row giving a starting estimate for the probability of the gene (in the corresponding row of the log ratio matrix/vector) not being differentially expressed and a starting estimate for the probability of the gene being differentially expressed. Each row should add up to 1.}
  \item{tol}{A scalar tolerance for relative convergence of the loglikelihood.}
  \item{iterlim}{The maximum number of iterations the EM is run for.}
}
\details{
A balanced dye swap is where a certain number of replicates have a particular dye to sample assigment and the same number of other replicates have the reversed assignment. Note in this case log ratios should be taken with numerators being the same sample and denominators the other sample, i.e. ratios should always be sample i/sample j rather than red dye/green dye for all replicates. }
\value{

A list including the following components
  \item{pdiff }{A vector with the estimated posterior probabilities of being in the group of differentially expressed genes.}
  \item{lRnorm }{A vector with the normalized (average-across-replicates) log ratios.}
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
\seealso{\code{\link{nudge2}},\code{\link{norm1a}},\code{\link{norm1b}},\code{\link{norm1c}},\code{\link{norm1d}},\code{\link{norm2c}},\code{\link{norm2d}}}
\examples{
data(like)
lR<-log(like[,1],2)-log(like[,2],2)
lI<-log(like[,1],2)+log(like[,2],2)

result<-nudge1(lR,lI)

data(hiv)
lR<-log(hiv[,1:4],2)-log(hiv[,5:8],2)
lI<-log(hiv[,1:4],2)+log(hiv[,5:8],2)

result<-nudge1(lR,lI,dye.swap=TRUE)

}
\keyword{cluster}
\keyword{loess}