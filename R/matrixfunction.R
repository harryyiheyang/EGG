mad=function(a){
  b=1.483*median(abs(a-median(a)))
  return(b)
}

trace=function(A){
  a=sum(diag(A))
  return(a)
}

matrixsqrt=function(A){
  fit=eigen(A)
  d=fit$value
  d1=d*0
  d1[d>0]=1/d[d>0]
  d=sqrt(d)
  d1=sqrt(d1)
  A=fit$vector%*%(t(fit$vector)*d)
  B=fit$vector%*%(t(fit$vector)*d1)
  C=list(w=A,wi=B)
  return(C)
}

entropyloss=function(A,B,eps=1e-8){
  C=A%*%B
  a=trace(C)
  b=eigen(C)$values
  b=Re(b)
  ind=which(b>eps)
  b=sum(log(b[ind]))
  return(a-b-ncol(A))
}

dtraceloss=function(A,B){
  a=0.5*trace(A%*%A%*%B)-trace(A)
  return(a)
}

positiveadj=function(A,min.eps=0.001){
  a=eigen(A)
  d=a$values
  d[d<min.eps]=0
  B=a$vectors%*%diag(d)%*%t(a$vectors)
  return(B)
}

cov2cor1=function(A,kappa=100){
  fit=eigen(A)
  d=fit$values
  eps=max(d)/kappa
  d[d<eps]=eps
  B=fit$vectors%*%(t(fit$vectors)*d)
  colnames(B)=rownames(B)=colnames(A)
  return(cov2cor(B))
}
