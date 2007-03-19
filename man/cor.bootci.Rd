\name{cor.bootci}
\alias{cor.bootci}
\title{Bootstrap Confidence Interval for Multivaraite Correlation}
\description{
This function calcualtes Bootstrap confidence interval for multivaraite correlation.
The procedure is very similair to those used to calcualte Bootstrape CI for
other parameters. See manuscript for detail.
}
\usage{
cor.bootci(x, y = NULL, m, G, alpha)
}
\arguments{
  \item{x}{data matrix, column represents samples (conditions), and row represents variables (genes), see example below for format information}
  \item{y}{optional, used when x and y are vectors}
  \item{m}{number of replicates}
  \item{G}{number of genes}
  \item{alpha}{significant level}
}
\details{
See manuscript.
}
\value{
  \item{upperCI}{Upper bound of CI}
  \item{lowerCI}{Lower bound of CI}
}
\references{Zhu, D and Li Y. 2007. Multivariate Correlation Estimator for Inferring Functional Relationships
from Replicated 'OMICS' data. Submitted.}
\author{Dongxiao Zhu and Youjuan Li}
\seealso{\code{\link{cor.LRtest1}}, \code{\link{cor.LRtest2}}, \code{\link{cor.test}}, \code{\link{permutest}}}
\examples{
library("CORREP")
d0 <- NULL
## sample size is set to 5 
for(l in 1:5)
d0 <- rbind(d0, rnorm(8))
## data must have row variance of 1 
d0.std <- apply(d0, 2, function(x) x/sd(x))
M <- cor.balance(t(d0.std), m = 2, G= 4)
pv.bootci <- cor.bootci(t(d0.std), m = 2, G= 4, alpha = 0.05)
}
\keyword{multivariate}
\keyword{cluster}
\keyword{models}
\keyword{htest}
