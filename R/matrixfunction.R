mad=function(a){
  b=1.483*median(abs(a-median(a)))
  return(b)
}

trace=function(A){
  a=sum(diag(A))
  return(a)
}

matrixsqrt=function(A){
  fit=matrixEigen(A)
  d=fit$value
  d1=d*0
  d1[d>0]=1/d[d>0]
  d=sqrt(d)
  d1=sqrt(d1)
  A=matrixMultiply(fit$vector,(t(fit$vector)*d))
  B=matrixMultiply(fit$vector,(t(fit$vector)*d1))
  C=list(w=A,wi=B)
  return(C)
}

entropyloss=function(A,B,eps=1e-8){
  C=matrixMultiply(A,B)
  a=trace(C)
  b=matrixEigen(C)$values
  b=Re(b)
  ind=which(b>eps)
  b=sum(log(b[ind]))
  return(a-b-ncol(A))
}

dtraceloss=function(A,B){
  a=0.5*trace(matrixListProduct(list(A,A,B)))-trace(A)
  return(a)
}

positiveadj=function(A,min.eps=0.001){
  a=matrixEigen(A)
  d=a$values
  d[d<min.eps]=0
  B=matrixMultiply(a$vector,(t(a$vector)*d))
  return(B)
}

cov2cor1=function(A,kappa=100){
  fit=matrixEigen(A)
  d=fit$values
  eps=max(d)/kappa
  d[d<eps]=eps
  B=matrixMultiply(fit$vector,(t(fit$vector)*d))
  colnames(B)=rownames(B)=colnames(A)
  return(cov2cor(B))
}
