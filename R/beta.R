.beta <- function(x,z,y,p,nknots,lambda,sigma_error=diag(1,length(y)),Dm){
  n <- length(y)
  df1 <- nknots+p
  Psi <- bs(z,df=df1,degree = p)
  Psi <- as.matrix(Psi[1:n,1:df1])
  x1 <- diag(x)
  if(lambda!=0){
    solve(t(Psi) %*% x1 %*% solve(sigma_error) %*% x1 %*% Psi +
            (1/lambda)*(nknots^(2*m-1))*t(Dm)%*%(Dm)) %*%t(Psi) %*% x1 %*% solve(sigma_error) %*% y
  }else
    solve(t(Psi) %*% x1 %*% solve(sigma_error) %*% x1 %*% Psi) %*%
    t(Psi) %*% x1 %*% solve(sigma_error) %*% y
}
