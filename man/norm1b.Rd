\name{norm1b}
\alias{norm1b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for normalizing the mean and variance (or just the variance) of single replicate log ratios}
\description{
This performs a robust loess normalization of the variance of the log ratios in a single replicate experiment by regressing the absolute (mean normalized) log ratios on the log intensities and using the fitted values to scale the (mean normalized) log ratio for each gene.}
\usage{
norm1b(logratio, logintensity, span1 = 0.6, span2 = 0.2, mean.norm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logratio}{A vector or single-column matrix of log (base 2) ratios of gene expressions in two samples, if mean.norm is FALSE the log ratios should be already mean normalized.}
  \item{logintensity}{A vector or single-column matrix of log (base 2) total intensities (defined as the product) of gene expressions in the two samples.}
  \item{span1}{Proportion of data used to fit the loess regression of the log ratios on the log intensities for the mean normalization.}
  \item{span2}{Proportion of data used to fit the loess regression of the absolute (mean normalized) log ratios on the log intensities for the variance normalization.}
  \item{mean.norm}{A logical value indicating whether or not a mean normalization should be performed prior to the variance normalization.}
}
\value{A vector or single-column matrix of mean and variance normalized log (base 2) ratios of gene expressions in two samples.
}
\references{N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays.  BMC Bioinformatics. 6, 173-186.

\url{http://www.biomedcentral.com/1471-2105/6/173}
}
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{norm1a}},\code{\link{norm1c}},\code{\link{norm1d}},\code{\link{norm2c}},\code{\link{norm2d}}}
\examples{
data(like)
lR<-log(like[,1],2)-log(like[,2],2)
lI<-log(like[,1],2)+log(like[,2],2)

lRnorm<-norm1b(lR,lI,mean.norm=TRUE)

}
\keyword{loess}
