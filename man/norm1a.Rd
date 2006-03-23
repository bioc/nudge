\name{norm1a}
\alias{norm1a}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for normalizing the mean of single replicate log ratios}
\description{
This is essentially the same as the lowess normalization suggested in the paper "Statistical methods for identifying differentially expressed genes in replicated cDNA microarrays" by Dudoit et al (2002), except the loess function is used instead of lowess and the recommended span is between 0.6 and 0.8. The normalization is done for each gene by subtracting from its log ratio the loess estimated mean for the log ratio based on the regression of log ratios on log intensities.
}
\usage{
norm1a(logratio, logintensity, span = 0.6)
}

\arguments{
  \item{logratio}{A vector or single-column matrix of log (base 2) ratios of gene expressions in two samples.}
  \item{logintensity}{A vector or single-column matrix of log (base 2) total intensities (defined as the product) of gene expressions in the two samples.}
  \item{span}{Proportion of data used to fit the loess regression of the log ratios on the log intensities.}
}
\value{
A vector or single-column matrix of mean normalized log ratios.
}
\references{
N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays. BMC Bioinformatics. 6, 173-186.
 
\url{http://www.biomedcentral.com/1471-2105/6/173}

S. Dudoit, Y. H. Yang, M. Callow and T. Speed (2002). Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. Stat. Sin. 12, 111-139. 
}
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{norm1b}},\code{\link{norm1c}},\code{\link{norm1d}},\code{\link{norm2c}},\code{\link{norm2d}}}
\examples{
data(like)
lR<-log(like[,1],2)-log(like[,2],2)
lI<-log(like[,1],2)+log(like[,2],2)

lRnorm<-norm1a(lR,lI)

}
\keyword{loess}
