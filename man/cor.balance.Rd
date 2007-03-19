\name{cor.balance}
\alias{cor.balance}
\title{Multivariate Correlation Estimator (Equal Number of Replicates)}
\description{
cor.balance estimates correlation matrix from replicated data assuming equal number of replicates.
The data must be formated in the right format (rows correspond to replicates, columns correspond
to conditions, see example below) and the variance of each row of the data MUST equal to 1 (see example below).  
}
\usage{
cor.balance(x, m, G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data matrix, column represents samples (conditions), and row represents variables (genes), see example below for format information}
  \item{m}{number of replicates for each variable (gene)}
  \item{G}{number of variables (genes)}
}
\details{
The multivariate correlation estimator assumes replicated omics data are iid samples
from the multivariate normal distribution. It is derived by maximizing the likelihood 
function. Note that each off-diagonal element in the returned correlation matrix (G by G) is
the average of off-diagonals of MLE of correlation matrix of a pair of variables (m by m). 
}
\value{
  A correlation matrix estimated for G variables (genes) 
}
\references{Zhu, D and Li Y. 2007. Multivariate Correlation Estimator for Inferring Functional Relationships
from Replicated 'OMICS' data. Submitted.}
\author{Dongxiao Zhu and Youjuan Li}
\seealso{\code{\link{cor.balance}}, \code{\link{cor}}}
\examples{
library("CORREP")
d0 <- NULL
for(l in 1:10)
d0 <- rbind(d0, rnorm(100))
## The simulated data corresponds to the real-world data of 25 genes and 10 conditions, each gene expression
## profiles was replicated 4 times.
d0<- t(d0)
## This step is to make the standard deviation of each replicate equals to 1
## so that we can model the covariance matrix as correlation matrix.
d0.std <- apply(d0, 1, function(x) x/sd(x))
M <- cor.balance(t(d0.std), m=4, G=25)
}
\keyword{multivariate}
\keyword{cluster}
\keyword{models}
\keyword{htest}
