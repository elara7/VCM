#' @export

val <- function(x,z,y,p,nknots,m=2,lambda=NULL,cv=FALSE,sigma_error=diag(1,length(y))){
  require(splines)
  if(any(is.na(x)) | any(is.na(z)) | any(is.na(y)))
    stop("NAs in data!")
  n <- length(y)
  df1 <- nknots+p
  Psi <- bs(z,df=df1,degree = p)
  Psi <- as.matrix(Psi[1:n,1:df1])
  D <- diag(1,df1)
  D[which(D==1)[-df1]+1] <- -1
  Dm <- D
  for(i in 1:(m-1)){Dm <- D %*% Dm}
  Dm <- Dm[-(1:m),]
  x1 <- diag(x)
  if(is.null(lambda)){
    lamb <- seq(0,10,0.01)
    mse <-rep(NA,length(lamb))
    for(i in 1:length(lamb)){
      mse[i] <- .valcv(x=x,z=z,y=y,p=p,nknots=nknots,m=m,
                      lambda=lamb[i],sigma_error=diag(1,length(y)),Dm=Dm)
    }
    lambda<- lamb[which(mse==min(mse))]
  }
  if(lambda==0){
    Q<- solve(t(Psi) %*% x1 %*% solve(sigma_error) %*% x1 %*% Psi)
  }else{
    Q<- solve(t(Psi) %*% x1 %*% solve(sigma_error) %*% x1 %*% Psi +
                (1/lambda)*(nknots^(2*m-1))*t(Dm)%*%(Dm))
  }
  BT <-  t(Psi) %*% x1 %*% solve(sigma_error)
  coef <- Q %*% BT %*% y
  theta <- Psi %*% coef
  fit <- x1 %*% theta
  residuals <- y-fit
  H <- t(BT) %*% Q %*% BT
  gcv <- sum((y-fit)^2)/(n-sum(diag(H)))^2
  res <- list(coef=drop(coef),lambda=lambda,GCV=gcv,fit=fit,
              theta=theta,residuals=residuals,degree=p,nknots=nknots,
              m=m,sigma_error=sigma_error)
  class(res) <- "val"
  res
}
