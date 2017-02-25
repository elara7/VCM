#' @export

LRTval=function(X, Z, Y, p=1,nknots=10,df,sigma.output=diag(rep(1, length(Y)))){
  require(lme4)
  require(RLRsim)
  require(nlme)
  n=length(Y)
  x1 <- diag(X)
  z1 <- rep(1,n)
  for (i in 1:p){z1=cbind(z1, z^i)}
  myknots = quantile(unique(z),seq(0,1,length=(nknots+2))[-c(1,(nknots+2))])
  z2 <- outer(z, myknots, FUN="-")
  z3 <- z2*(z2>0)
  if(p>1) {z3=z3^p}
  A1 <- x1 %*% z1
  A2 <- x1 %*% z3
  Xnames = paste("X",1:ncol(A1),sep="")
  Znames = paste("Z",1:ncol(A2),sep="")
  fixed.model = as.formula(paste("DATA.temp ~ -1+",
                                 paste(paste("X",1:ncol(A1),sep=""),collapse="+")))
  fixed.model2 = as.formula(paste("DATA.temp ~ -1+",
                                  paste(paste("X",1:(ncol(A1)-df),sep=""),collapse="+")))
  random.model = as.formula(paste("~-1+",paste(paste("Z",1:ncol(A2),sep=""),collapse="+")))

  DATA.output = as.vector( sigma.output  %*% Y )
  hat1 = sigma.output %*% A1 ;
  hat2 = sigma.output %*% A2

  DATA.temp= DATA.output
  colnames(hat1)=Xnames
  colnames(hat2)=Znames
  subject<-rep(1,n)
  ALLDATA  = data.frame(cbind(subject,DATA.temp, hat1, hat2))
  mA = lme(fixed=fixed.model, data=ALLDATA,
           random=list(subject = pdIdent(random.model)), method="ML")
  m0 = lm(fixed.model2, data=ALLDATA)
  obs.LRT = as.numeric(2*(logLik(mA)-logLik(m0)))
  pvalue = pchisq(obs.LRT,1,lower.tail = FALSE)
  res <- list(LRT =obs.LRT,p.value = pvalue)
  res
}
