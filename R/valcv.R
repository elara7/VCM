.valcv <- function(x,z,y,p,nknots,m=2,lambda,sigma_error=diag(1,length(y)),Dm){
  df1 <- nknots+p
  fold <- 5
  set.seed(1234)
  k = sample(rep(1:fold,length=n))
  fit1 <- rep(NA,n)
  for (j in 1:fold){
    testset = (1:n)[k==j]
    coef <- .beta(x=x[-testset],z=z[-testset],y=y[-testset],
                 p=p,nknots = nknots,m,lambda=lambda,Dm=Dm)
    l <- length(x[testset])
    NPsi <- bs(z[testset],df=df1,degree = p)
    Npsi <- as.matrix(NPsi[1:l,1:df1])
    Ntheta <- Npsi %*% coef
    fit1[testset] <- diag(x[testset]) %*% Ntheta
  }
  mse <- mean((y-fit1)^2)
  return(mse)
}
