`cor.unbalance` <-
function(x, m1, m2)
{
    n <- dim(x)[2]
    if(m1+m2 != nrow(x))
    stop("incorrect number of replicates, check m1 and m2")
    muhat <- c(rep(mean(x[1:m1,]), m1), rep(mean(x[(m1+1):(m1+m2),]), m2))
    M <- matrix(0, nrow(x), nrow(x))
    for (k in 1:n)
        M <- M + (x[,k] - muhat)%o%(x[,k] - muhat)
    Sigmahat <- M/n
    cor.unbal <- mean(c(Sigmahat[1:m1,(m1+1):(m1+m2)], Sigmahat[(m1+1):(m1+m2),1:m1]))
    cor.unbal
}

