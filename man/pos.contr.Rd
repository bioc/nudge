\name{pos.contr}
\docType{data}
\alias{pos.contr}
\title{Indicator of known truly differentially expressed genes in the hiv dataset}
\description{
The ``hiv'' section of this dataset consists of cDNA from CD4+ T cell lines at 1 hour after infection with HIV-1BRU. There are 4608 genes' expression levels. There are four replicates with balanced dye sways i.e. two of the replicates have sample 1 labeled with red and sample 2 labeled with green and the other two have the opposite labeling scheme. It is useful in testing the sensitivity of methods identifying differentially expressed genes, since there are 13 genes known to be differentially expressed (HIV-1 genes), called positive controls, identified in the vector of logical values ``pos.contr''}
\usage{data(hiv)}
\format{A logical vector containing 4608 indicators of known differential expression}
\source{Bumgarner Lab}
\references{
N. Dean and A. E. Raftery (2005). Normal uniform mixture differential gene expression detection for cDNA microarrays. BMC Bioinformatics. 6, 173-186.

\url{http://www.biomedcentral.com/1471-2105/6/173}

A.B. van't Wout, G. K. Lehrma, S. A. Mikheeva, G. C. O'Keeffe, M. G. Katze, R. E. Bumgarner, G. K. Geiss and J. I. Mullins (2003). Cellular gene expression upon human immunodeficiency virus type 1 infection of CD4+ T-Cell lines. J. Virol. 77, 1392-1402.
}
\seealso{\code{\link{hiv}},\code{\link{neg.contr}}}
\keyword{datasets}
