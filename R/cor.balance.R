`cor.balance` <-
function(x, m, G)
{
    p <- dim(x)[1]
    n <- dim(x)[2]
    D <- p/m
    if(D != G)
    stop("unqueal number of replicates or incorrect number of replicates or incorrect number of genes (random variables)")
    M <- matrix(NA, D, D)

    for(i in 1:(D-1))
    {
        for(j in (i+1):D)
        {
        subdat <- rbind(x[((i-1)*m + 1):((i-1)*m + m),], x[((j-1)*m + 1):((j-1)*m + m),])
        muhat <- c(rep(mean(subdat[1:m,]), m), rep(mean(subdat[(m+1):(2*m),]), m))
        MB <- matrix(0, nrow(subdat), nrow(subdat))
        for (k in 1:n)
            MB <- MB + (subdat[,k] - muhat)%o%(subdat[,k] - muhat)
        Sigmahat <- (MB/n)
        M[i,j] <- mean(c(Sigmahat[1:m,(m+1):(2*m)],Sigmahat[(m+1):(2*m),1:m]))
        M[j,i] <- M[i,j]
        }
    }
diag(M) <- 1
M
}

