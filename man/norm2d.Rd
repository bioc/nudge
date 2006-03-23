\name{norm2d}
\alias{norm2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for normalizing the mean and variance of average-across-replicates log ratio differences}
\description{
This normalization is used when the two samples (control and treatment, say) are not being directly compared on the slides but instead are being compared to a common reference sample. The quantity of interest for each gene is thus the average difference between control and treatment log ratios. This function performs a robust normalization of the variance of the (mean normalized) average-across-replicates log ratio differences by scaling the (mean normalized) average-across-replicates log ratio difference for each gene either by the standard deviation of the log ratio differences for that gene across replicates (if bigger than the absolute (mean normalized) average-across-replicates log ratio difference) or scaling by a constant (a quantile of the distribution of standard deviations of (mean normalized) average-across-replicates log ratio differences for all genes whose standard deviation was bigger than their absolute (mean normalized) average-across-replicates log ratio difference.}
\usage{
norm2d(control.logratio, txt.logratio, control.logintensity, txt.logintensity,
span = 0.6, quant = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{control.logratio}{A multiple-column matrix of replicates of log (base 2) ratios of gene expressions for the control versus reference slides.}
  \item{txt.logratio}{A multiple-column matrix of replicates of log (base 2) ratios of gene expressions for the treatment versus reference slides.}
  \item{control.logintensity}{A multiple-column matrix of replicates of log (base 2) total intensities (defined as the product) of gene expressions for the control versus reference slides.}
  \item{txt.logintensity}{A multiple-column matrix of replicates of log (base 2) total intensities (defined as the product) of gene expressions for the treatment versus reference slides.}
  \item{span}{Proportion of data used to fit the loess regression of the average-across-replicates log ratio differences on the average-across-replicates log intensities.}
  \item{quant}{Quantile to be used from the distribution of standard deviations of log ratio differences across replicates for all genes whose standard deviation was smaller than their absolute (mean normalized) average-across-replicates log ratio difference.}
}
\value{
A vector of mean and variance normalized average-across-replicates log ratio differences.}
\references{N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays.  BMC Bioinformatics. 6, 173-186.

\url{http://www.biomedcentral.com/1471-2105/6/173}

S. Dudoit, Y. H. Yang, M. Callow and T. Speed (2002). Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. Stat. Sin. 12, 111-139. 
}
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{norm2c}},\code{\link{norm1a}},\code{\link{norm1b}},\code{\link{norm1c}},\code{\link{norm1d}}}
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

lRnorm<-norm2d(lRctl,lRtxt,lIctl,lItxt)

}
\keyword{loess}
