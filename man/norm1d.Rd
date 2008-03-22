\name{norm1d}
\alias{norm1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for normalizing the mean and variance of average-across-replicates log ratios}
\description{
This performs a robust normalization of the variance of the (mean normalized) average-across-replicates log ratios by scaling the (mean normalized) average-across-replicates log ratio for each gene either by the standard deviation of the log ratios for that gene across replicates (if bigger than its absolute (mean normalized) average-across-replicates log ratio) or scaling by a constant (a quantile of the distribution of standard deviations of log ratios across replicates for all genes whose standard deviation was bigger than their absolute (mean normalized) average-across-replicates log ratio.}
\usage{
norm1d(logratio, logintensity, span = 0.6, quant = 0.99, dye.swap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logratio}{A multiple-column matrix of replicates of log (base 2) ratios of gene expressions in two samples.}
  \item{logintensity}{A multiple-column matrix of replicates of log (base 2) total intensities (defined as the product) of gene expressions in two samples.}
  \item{span}{Proportion of data used to fit the loess regression of the average-across-replicates log ratios on the average-across-replicates log total intensities.}
  \item{quant}{Quantile to be used from the distribution of standard deviations of log ratios across replicates for all genes whose standard deviation was smaller than their absolute (mean normalized) average-across-replicates log ratio.}
  \item{dye.swap}{A logical value indicating whether or not the data consists of a balanced dye swap (if FALSE a mean normalization will be performed prior to the variance normalization rather than a simple averaging across replicates).}
}
\details{
A balanced dye swap is a data set were the data are split into two sets of replicates (same number of replicates in each set) where one set has the reverse dye assignment of the other set.
}
\value{A vector of mean and variance normalized average-across-replicates log ratios.

}
\references{N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays.  BMC Bioinformatics. 6, 173-186. 
 
\url{http://www.biomedcentral.com/1471-2105/6/173}
}
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{norm1c}},\code{\link{norm1a}},\code{\link{norm1b}},\code{\link{norm2c}},\code{\link{norm2d}}}
\examples{
data(hiv)
lR<-log(hiv[,1:4],2)-log(hiv[,5:8],2)
lI<-log(hiv[,1:4],2)+log(hiv[,5:8],2)

lRnorm<-norm1d(lR,lI,dye.swap=TRUE)



}
\keyword{loess}
