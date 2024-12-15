spearmancov=function(A){
  p=dim(A)[2]
  s=c(1:p)
  for(i in 1:p){
      s[i]=median(abs(median(A[,i])-A[,i]))*1.483
  }
  R=cor(A,method="spearman")
  R=2*sin(R*pi/6)
  S=t(R*s)*s
  return(S)
}
