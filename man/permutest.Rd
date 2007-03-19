\name{permutest}
\alias{permutest}
\title{Permutation Test P-value for Multivaraite Correlation}
\description{
This function calculates p-values of the multivariate correlation estimator
by enumerating all permutations. We recommend using Likehood Ratio Test implemented
in function cor.LRtest if your data has moderate to large sample size (>5)
The procedure is same as those permutation tests for Pearson correlation coefficient
or other parameters. Since the approximation of null distribution requires enumerating all permutations.
The computational burden increases in $n^2$.   
}
\usage{
permutest(x, y=NULL, m, G)
}
\arguments{
  \item{x}{data matrix, column represents samples (conditions), and row represents variables (genes), see example below for format information}
  \item{y}{optional, used when x and y are vectors}
  \item{m}{number of replicates}
  \item{G}{number of genes}
}
\details{
See manuscript.
}
\value{
  \item{PV}{P-values of permutation tests}
}
\references{Zhu, D and Li Y. 2007. Multivariate Correlation Estimator for Inferring Functional Relationships
from Replicated 'OMICS' data. Submitted.}
\author{Dongxiao Zhu and Youjuan Li}
\seealso{\code{\link{cor.LRtest}}, \code{\link{cor.LRtest.std}}, \code{\link{cor.test}}}
\examples{
library("CORREP")
library("e1071")
d0 <- NULL
## sample size is set to 5, it takes about a min to finish 
for(l in 1:5)
d0 <- rbind(d0, rnorm(100))
## data must have row variance of 1 
d0.std <- apply(d0, 2, function(x) x/sd(x))
M <- cor.balance(t(d0.std), m = 4, G= 25)
M.pv <- permutest(t(d0.std), m = 4, G= 25)
}
\keyword{multivariate}
\keyword{cluster}
\keyword{models}
\keyword{htest}
