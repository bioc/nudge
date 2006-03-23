\name{like}
\docType{data}
\alias{like}
\title{Microarray intensity expression levels for like-like data}
\description{
This data set gives the expression levels for the same sample tested against itself, i.e. one sample split into two and labeled with both red and green dyes. Useful for testing false positive rate of methods for detecting differential expression since there are, theoretically no differentially expressed genes. There are 7680 different genes' expression levels measured.}
\usage{data(like)}
\format{A matrix containing 7680 observations on two variables}
\source{Bumgarner Lab}
\references{
N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays. BMC Bioinformatics. 6, 173-186.
 
\url{http://www.biomedcentral.com/1471-2105/6/173}
}
\keyword{datasets}