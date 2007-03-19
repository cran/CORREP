\name{cor.unbalance}
\alias{cor.unbalance}
\title{Multivariate Correlation Estimator (Unequal Number of Replicates)}
\description{
cor.unbalance estimates correlation from replicated data of unequal number of replicates.
different from \code{\link{cor.balance}}, \code{\link{cor.unbalance}} takes a pair of varaibles 
at a time because of unequal number of replicates. the variance of each row of the data MUST equal to 1 (see example below)    
}
\usage{
cor.unbalance(x, m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data matrix, column represents samples (conditions), and row represents variables (genes), see example below for format information}
  \item{m1}{number of replicates for one variable (gene)}
  \item{m2}{number of replicates for another variable (gene)}
}
\details{
The multivariate correlation estimator assumes replicated omics data are iid samples
from the multivariate normal distribution. It is derived by maximizing the likelihood 
function. Note that the off-diagonal elements in the returned correlation matrix (G by G) is
the average of off-diagonals of MLE of correlation matrix of a pair of variables (m1+m2 by m1+m2). 
}
\value{
  A correlation matrix containing only one distinct correlation coefficient for the pair of variables (genes) 
}
\references{Zhu, D and Li Y. 2007. Multivariate Correlation Estimator for Inferring Functional Relationships
from Replicated 'OMICS' data. Submitted.}
\author{Dongxiao Zhu and Youjuan Li}
\seealso{\code{\link{cor.unbalance}}, \code{\link{cor}}}
\examples{
library("CORREP")
d0 <- NULL
for(l in 1:10)
d0 <- rbind(d0, rnorm(8))
## The simulated data corresponds to the real-world data of 2 genes and 10 conditions, gene expression
## profiles were replicated 3 and 5 times. 
## Note this function can only take calculate correlation matrix between two genes at a time.
d0<- t(d0)
## This step is to make the standard deviation of each replicate equal to 1
## so that we can model the covariance matrix as correlation matrix.
d0.std <- apply(d0, 1, function(x) x/sd(x))
M <- cor.unbalance(t(d0.std), m1=3, m2=5)
}
\keyword{multivariate}
\keyword{cluster}
\keyword{models}
\keyword{htest}
