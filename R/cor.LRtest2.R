`cor.LRtest2` <-
function(x, m1, m2)
{
    ## m: number of replicates
    p <- dim(x)[1]
    n <- dim(x)[2]
    stopifnot(p == m1+m2)

    ## get grant means muhat_x, muhat_y
    mu.x <- mean(x[1:m1,])
    mu.y <- mean(x[(m1+1):(m1+m2),])
    mux <-  c(rep(mu.x, m1))
    muy <-  c(rep(mu.y, m2))
    mu <- c(rep(mu.x, m1), rep(mu.y, m2))

    Sigma.x <- matrix(0, m1, m1)
    Sigma.y <- matrix(0, m2, m2)
    Sigma.null <- matrix(0, p, p)
    Sigma.alter <- matrix(0, p, p)
    xnew1 <- NULL
    xnew2 <- NULL
    xnew <- NULL
    ## get Sigmahat_x, Sigmahat_y, Sigma_null, Sigma_alternative
    xnew1 <- x[1:m1,]-mux
    xnew2 <- x[(m1+1):(m1+m2),]-muy
    xnew <- x-mu

    for (j in 1:n)
    {
        Sigma.x <- Sigma.x  +  xnew1[,j]%*%t(xnew1 [,j])
        Sigma.y <- Sigma.y  +  xnew2[,j]%*%t(xnew2 [,j])
        Sigma.alter <- Sigma.alter + xnew[,j]%*%t(xnew[,j])
    }
    Sigma.x <- Sigma.x/n
    Sigma.y <- Sigma.y/n
    Sigma.alter <- Sigma.alter/n
    Sigmax <- matrix(0, m1, m1)
    Sigmay <- matrix(0, m2, m2)
    Sigmaalter <- matrix(0, p, p)

    for (i in 1:m1)
    for (j in 1:m1)
        Sigmax[i,j] <- Sigma.x[i,j]/sqrt(Sigma.x[i,i]*Sigma.x[j,j])

    for (i in 1:m2)
    for (j in 1:m2)
        Sigmay[i,j] <- Sigma.y[i,j]/sqrt(Sigma.y[i,i]*Sigma.y[j,j])

    for (i in 1:p)
    for (j in 1:p)
        Sigmaalter[i,j] <- Sigma.alter[i,j]/sqrt(Sigma.alter[i,i]*Sigma.alter[j,j])

    Sigma.null[1:m1,1:m1] <- Sigma.x
    Sigma.null[(m1+1):(m1+m2), (m1+1):(m1+m2)] <- Sigma.y

    M <- solve(Sigma.null) %*% Sigma.alter
    G <- n* (sum(diag(M))-log(det(M))-(m1+m2))

    p.value <- 1- pchisq(G, df=2*m1*m2)
    p.value
}
