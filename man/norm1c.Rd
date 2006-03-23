\name{norm1c}
\alias{norm1c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for normalizing the mean of average-across-replicates log ratios}
\description{
This is essentially the same as the lowess normalization suggested in the paper "Statistical methods for identifying differentially expressed genes in replicated cDNA microarrays" by Dudoit et al (2002), except the loess function and average-across-replicates log ratios are used and the recommended span is between 0.6 and 0.8. The normalization is done for each gene by subtracting from its average-across-replicates log ratio the loess estimated mean for average-across-replicates log ratio based on the loess regression of average-across-replicates log ratios on average-across-replicates log total intensities.
}
\usage{
norm1c(logratio, logintensity, span = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logratio}{A multiple-column matrix of replicates of log (base 2) ratios of gene expressions in two samples.}
  \item{logintensity}{A multiple-column matrix of replicates of log (base 2) total intensities (defined as the product) of gene expressions in two samples.}
  \item{span}{Proportion of data used to fit the loess regression of the average-across-replicates log ratios on the average-across-replicates log total intensities}
}

\value{
A vector of mean normalized average-across-replicates log ratios.
}
\references{N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays.  BMC Bioinformatics. 6, 173-186. 

\url{http://www.biomedcentral.com/1471-2105/6/173}

S. Dudoit, Y. H. Yang, M. Callow and T. Speed (2002). Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. Stat. Sin. 12, 111-139. 
}
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{norm1d}},\code{\link{norm1a}},\code{\link{norm1b}},\code{\link{norm2c}},\code{\link{norm2d}} }
\examples{
data(hiv)
lR<-log(hiv[,1:4],2)-log(hiv[,5:8],2)
lI<-log(hiv[,1:4],2)+log(hiv[,5:8],2)

lRnorm<-norm1c(lR,lI)

}
\keyword{loess}
