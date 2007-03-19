\name{gal_all}
\alias{gal_all}
\docType{data}
\title{Galactose Metabolism Data}
\description{
The data is compiled by Mario Medvedovic et al, 2003 based on the original full
data reported in Ideker et al, 2001. There are a total of 205 rows (genes), 
20 experiments, and 4 repeated measurements in the data. There are 4 classes 
(which correspond to functional categories). The data contains approximately 8% 
of missing data. The missing values were filled by applying k-nearest neighbor 
(k = 12) to impute all the missing values.   
}
\usage{data(gal_all)}
\format{
  A data frame with 205 variables on the following 80 replicated observations.
  \describe{
    \item{\code{wtRG1}}{a numeric vector}
    \item{\code{wtRG2}}{a numeric vector}
    \item{\code{wtRG3}}{a numeric vector}
    \item{\code{wtRG4}}{a numeric vector}
    \item{\code{gal1RG1}}{a numeric vector}
    \item{\code{gal1RG2}}{a numeric vector}
    \item{\code{gal1RG3}}{a numeric vector}
    \item{\code{gal1RG4}}{a numeric vector}
    \item{\code{gal2RG1}}{a numeric vector}
    \item{\code{gal2RG2}}{a numeric vector}
    \item{\code{gal2RG3}}{a numeric vector}
    \item{\code{gal2RG4}}{a numeric vector}
    \item{\code{gal3RG1}}{a numeric vector}
    \item{\code{gal3RG2}}{a numeric vector}
    \item{\code{gal3RG3}}{a numeric vector}
    \item{\code{gal3RG4}}{a numeric vector}
    \item{\code{gal4RG1}}{a numeric vector}
    \item{\code{gal4RG2}}{a numeric vector}
    \item{\code{gal4RG3}}{a numeric vector}
    \item{\code{gal4RG4}}{a numeric vector}
    \item{\code{gal5RG1}}{a numeric vector}
    \item{\code{gal5RG2}}{a numeric vector}
    \item{\code{gal5RG3}}{a numeric vector}
    \item{\code{gal5RG4}}{a numeric vector}
    \item{\code{gal6RG1}}{a numeric vector}
    \item{\code{gal6RG2}}{a numeric vector}
    \item{\code{gal6RG3}}{a numeric vector}
    \item{\code{gal6RG4}}{a numeric vector}
    \item{\code{gal7RG1}}{a numeric vector}
    \item{\code{gal7RG2}}{a numeric vector}
    \item{\code{gal7RG3}}{a numeric vector}
    \item{\code{gal7RG4}}{a numeric vector}
    \item{\code{gal10RG1}}{a numeric vector}
    \item{\code{gal10RG2}}{a numeric vector}
    \item{\code{gal10RG3}}{a numeric vector}
    \item{\code{gal10RG4}}{a numeric vector}
    \item{\code{gal80RG1}}{a numeric vector}
    \item{\code{gal80RG2}}{a numeric vector}
    \item{\code{gal80RG3}}{a numeric vector}
    \item{\code{gal80RG4}}{a numeric vector}
    \item{\code{wtR1}}{a numeric vector}
    \item{\code{wtR2}}{a numeric vector}
    \item{\code{wtR3}}{a numeric vector}
    \item{\code{wtR4}}{a numeric vector}
    \item{\code{gal1R1}}{a numeric vector}
    \item{\code{gal1R2}}{a numeric vector}
    \item{\code{gal1R3}}{a numeric vector}
    \item{\code{gal1R4}}{a numeric vector}
    \item{\code{gal2R1}}{a numeric vector}
    \item{\code{gal2R2}}{a numeric vector}
    \item{\code{gal2R3}}{a numeric vector}
    \item{\code{gal2R4}}{a numeric vector}
    \item{\code{gal3R1}}{a numeric vector}
    \item{\code{gal3R2}}{a numeric vector}
    \item{\code{gal3R3}}{a numeric vector}
    \item{\code{gal3R4}}{a numeric vector}
    \item{\code{gal4R1}}{a numeric vector}
    \item{\code{gal4R2}}{a numeric vector}
    \item{\code{gal4R3}}{a numeric vector}
    \item{\code{gal4R4}}{a numeric vector}
    \item{\code{gal5R1}}{a numeric vector}
    \item{\code{gal5R2}}{a numeric vector}
    \item{\code{gal5R3}}{a numeric vector}
    \item{\code{gal5R4}}{a numeric vector}
    \item{\code{gal6R1}}{a numeric vector}
    \item{\code{gal6R2}}{a numeric vector}
    \item{\code{gal6R3}}{a numeric vector}
    \item{\code{gal6R4}}{a numeric vector}
    \item{\code{gal7R1}}{a numeric vector}
    \item{\code{gal7R2}}{a numeric vector}
    \item{\code{gal7R3}}{a numeric vector}
    \item{\code{gal7R4}}{a numeric vector}
    \item{\code{gal10R1}}{a numeric vector}
    \item{\code{gal10R2}}{a numeric vector}
    \item{\code{gal10R3}}{a numeric vector}
    \item{\code{gal10R4}}{a numeric vector}
    \item{\code{gal80R1}}{a numeric vector}
    \item{\code{gal80R2}}{a numeric vector}
    \item{\code{gal80R3}}{a numeric vector}
    \item{\code{gal80R4}}{a numeric vector}
  }
}
\details{
The 205 genes have been classified into four functional classes based on their
GO annotations. In the data examaple provided in the vignette, we assume the 
four classes as true memberships (external knowledge) and use it to evaluate
the performances of different correlation measured based clustering methods.
}
\source{
http://expression.microslu.washington.edu/expression/kayee/medvedovic2003/medvedovic\_bioinf2003.html
}
\references{
Medvedovic M, Yeung KY and Bumgarner RE. 2004. Bayesian Mixture Model Based Clustering of Replicated Microarray Data.
Bioinformatics, 22;20(8):1222-32.
Ideker, T., Thorsson, V., Siegel, A. and Hood, L. Testing for Differentially-Expressed Genes by Maximum-Likelihood 
Analysis of DNA Microarray Data. Journal of Computational Biology 7: 805-817 (2000).
}
\examples{
data(gal_all)
## maybe str(gal_all) ; plot(gal_all) ...
}
\keyword{datasets}
