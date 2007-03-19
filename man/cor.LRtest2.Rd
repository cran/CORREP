\name{cor.LRtest2}
\alias{cor.LRtest2}
\title{Maximum Likehood Ratio Test for Multivariate Correlation Estimator}
\description{
Performs LRT to test if multivariate correlation vanishes. Note
this code `standardizes' estimated correlation matrix to make sure its determinant is
positive.}
\usage{
cor.LRtest2(x, m1, m2)
}
\arguments{
  \item{x}{data matrix, column represents samples (conditions), and row represents variables (genes), see example below for format information}
  \item{m1}{number of replicates for gene X}
  \item{m2}{number of replicates for gene Y}
}
\details{
Under the multivaraite normal distribution assumption, the column vector of the 
data is iid sample. We test the followiing hypothesis: H0: Z ~ N(mu, Sigma0), H1
~ N(mu, Sigma1). Let M = Inverse(Sigma0)*Sigma1, the likelihood ratio test 
statistic $G^2$, is, n[trace(M)-log(det(M))-(m1+m2)]. Under the H0, $G^2$ follows a chi-square
distribution with $2m1*m2$ degree of freedom. In some case, the determinant of M is negative
so that the log(det(M)) return NaN. There are two ways to deal with this problem,
one, those M's whose determinant are negative tend to consist of very small correlations
that are biological irrelevant. Therefore, we can simply ignore those gene pairs that
the determinant correlation matrix M is negative. Second, we can `standardize' the M
to make the determinant positive (implemented in function cor.LRtest2). In most of cases,
we recommend using function cor.LRtest1. Use cor.LRtest2 only you know what you are doing.      
}
\value{
  \item{p.value}{The p-value of the LRT}
}
\references{Zhu, D and Li Y. 2007. Multivariate Correlation Estimator for Inferring Functional Relationships
from Replicated 'OMICS' data. Submitted.}
\author{Youjuan Li and Dongxiao Zhu}
\seealso{\code{\link{cor.LRtest1}}, \code{\link{cor.test}}}
\examples{
library("CORREP")
library("MASS")
Sigma <- matrix(c(1, 0.8, .5,.5, 0.8, 1,
0.5,0.5,0.5,0.5,1,0.6,0.5,0.5,0.6,1),4,4)
dat <- mvrnorm(50, mu=c(0,0,0,0), Sigma)
dat.std <- apply(dat,2,function(x) x/sd(x))
cor.LRtest2 (t(dat.std), m1=2, m2=2)
}
\keyword{multivariate}
\keyword{cluster}
\keyword{models}
\keyword{htest}
