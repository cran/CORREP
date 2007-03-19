\name{true.member}
\alias{true.member}
\docType{data}
\title{The Pre-defined Class Memberships of 205 Genes According to GO Annotation}
\description{
The four classes are: Biosynthesis; Energy pathways; Nucleobase and	Transport
It is used as true clusters (external knowledge) to compare performance of 
different clustering methods.
}
\usage{data(true.member)}
\source{
http://expression.microslu.washington.edu/expression/kayee/medvedovic2003/medvedovic\_bioinf2003.html
}
\references{
Medvedovic M, Yeung KY and Bumgarner RE. 2004. Bayesian Mixture Model Based Clustering of Replicated Microarray Data.
Bioinformatics, 22;20(8):1222-32.
}
\examples{
data(true.member)
## maybe str(true.member) ; plot(true.member) ...
}
\keyword{datasets}
