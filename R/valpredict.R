#' @export

valpredict <- function(model,newx,newz){
  if(class(model)!="val")
    stop("the model must be varying-coffecient!")
  n <- length(newx)
  p <- model$degree
  nknots <- model$nknots
  df <- p+nknots
  Psi <- bs(newz,df=df,degree = p)
  Psi <- as.matrix(Psi[1:n,1:df])
  theta <- Psi %*% model$coef
  predict <- diag(newx) %*% theta
  predict
}
